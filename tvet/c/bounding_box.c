#include "vector_math.h"
#include "bounding_box.h"
#include <math.h>
#include <stdlib.h>

void boundingBox(
    const int (*faces)[3], size_t nof_faces,
    const double (*nodes)[3], size_t nof_nodes,
    const double s[3],
    double (*boxes)[4]
) {
    double *u = (double *) malloc(nof_nodes * sizeof(double));
    double *v = (double *) malloc(nof_nodes * sizeof(double));

    if(nof_nodes != 0 && (u == NULL || v == NULL)) {
        free(u); free(v);
        return;
    }

    const double hatw[3] = {s[0], s[1], s[2]};
    double hatu[3] = {-s[1], s[0], 0.0};

    if(!normalizeVector3D(hatu)) {
        // s parallel with z axis
        hatu[0] = 1.0;
        hatu[1] = 0.0;
        hatu[2] = 0.0;
    }

    double hatv[3];
    crossProduct3D(hatu, hatw, hatv);

    // Flip sign to match original fortran code
    hatv[0] = -hatv[0];
    hatv[1] = -hatv[1];
    hatv[2] = -hatv[2];

    // Project each node into (u, v) coordinates
    for(size_t i = 0; i < nof_nodes; i++) {
        u[i] = dotProduct3D(hatu, nodes[i]);
        v[i] = dotProduct3D(hatv, nodes[i]);
    }

    // For every face compute min/max of its three projected vertices
    for(size_t j = 0; j < nof_faces; j++) {
        int j0 = faces[j][0];
        int j1 = faces[j][1];
        int j2 = faces[j][2];

        double u0 = u[j0], u1 = u[j1], u2 = u[j2];
        double v0 = v[j0], v1 = v[j1], v2 = v[j2];

        double u_min = u0;
        if(u1 < u_min) {u_min = u1;}
        if(u2 < u_min) {u_min = u2;}

        double u_max = u0;
        if(u1 > u_max) {u_max = u1;}
        if(u2 > u_max) {u_max = u2;}

        double v_min = v0;
        if(v1 < v_min) {v_min = v1;}
        if(v2 < v_min) {v_min = v2;}

        double v_max = v0;
        if(v1 > v_max) {v_max = v1;}
        if(v2 > v_max) {v_max = v2;}

        boxes[j][0] = u_min;
        boxes[j][1] = u_max;
        boxes[j][2] = v_min;
        boxes[j][3] = v_max;
    }

    free(u); free(v);
}
