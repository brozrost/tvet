#include "shadowing.h"
#include "intersect_AB_t.h"
#include "bounding_box.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

void mu(
    const double (*normals)[3], size_t nof_faces, const double s[3], double *mu_i
) {
    #pragma omp parallel for if(nof_faces > 256)
    for(size_t i = 0; i < nof_faces; i++) {
        double d = normals[i][0]*s[0] + normals[i][1]*s[1] + normals[i][2]*s[2];
        mu_i[i] = (d > 0.0) ? d : 0.0;
    }
}

void non(
    const double *mu_i, const double *mu_e, size_t nof_faces, double *nu_i, double *nu_e
) {
    for(size_t i = 0; i < nof_faces; i++) {
        nu_i[i] = 1.0;
        nu_e[i] = 1.0;
    }

    #pragma omp parallel for if(nof_faces > 256)
    for(size_t i = 0; i < nof_faces; i++) {
        if(mu_i[i] == 0.0 || mu_e[i] == 0.0) {
            nu_i[i] = 0.0;
            nu_e[i] = 0.0;
        }
    }
}

void nu(
    const int (*faces)[3], size_t nof_faces,
    const double (*nodes)[3], size_t nof_nodes,
    const double (*normals)[3],
    const double (*centers)[3],
    const double s[3],
    double *nu_i
) {
    double (*boxes)[4] = (double (*)[4]) malloc(nof_faces * sizeof(*boxes));
    if(!boxes) {return;}

    boundingBox(faces, nof_faces, nodes, nof_nodes, s, boxes);

    #pragma omp parallel for if(nof_faces > 64) schedule(static)
    for(size_t i = 0; i < nof_faces; i++) {
        if(nu_i[i] <= 0.0) {continue;}

        double A[3] = {
            centers[i][0],
            centers[i][1],
            centers[i][2]
        };

        for(size_t j = 0; j < nof_faces && nu_i[i] > 0.0; j++) {
            if(i == j) {continue;}
            if(nu_i[j] <= 0.0) {continue;}

            if((boxes[j][1] < boxes[i][0]) || (boxes[j][0] > boxes[i][1]) ||
                (boxes[j][3] < boxes[i][2]) || (boxes[j][2] > boxes[i][3])
            ) {
                continue;
            }

            const double *ni = normals[i];
            double C[3] = {
                centers[j][0] - centers[i][0],
                centers[j][1] - centers[i][1],
                centers[j][2] - centers[i][2],
            };

            double tmp = C[0]*ni[0] + C[1]*ni[1] + C[2]*ni[2];
            if(tmp <= 0.0) {continue;}

            double t[3][3];
            for(int k = 0; k < 3; k++) {
                int vidx = faces[j][k];
                const double *p = nodes[vidx];

                t[k][0] = p[0];
                t[k][1] = p[1];
                t[k][2] = p[2];
            }

            bool has_solution = false;

            intersect_AB_t(A, s, t, C, &has_solution);

            if(has_solution) {nu_i[i] = 0.0;}
        }
    }

    free(boxes);
}
