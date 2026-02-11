#include "vector_math.h"
#include "intersect_AB_t.h"
#include <stdbool.h>

// t is a 3x3 matrix flattened to length 9 in column major
#define T(i,j) (t[(i) + 3*(j)])

/// @brief Computes the intersection of a line and a triangle in 3D space. 
/// @param A Point on the line in 3D space.
/// @param B Direction vector of the line.
/// @param t Triangle vertecies. 
/// @param C Output parameter, intersection point.
/// @param has_solution Boolean result. True if intersection exists.
void intersect_AB_t(
    const double A[3], const double B[3], const double t[9], double C[3], bool *has_solution
) {
    double area[3];
    const double EPS = 0.0; 

    // i is the current vertex, i1 and i2 are the remaining two
    for(int i = 0; i < 3; i++) {
        int i1 = (i + 1) % 3;
        int i2 = (i + 2) % 3;

        // Vector from point A to vertex i1
        double v1[3] = {T(i1, 0) - A[0], T(i1, 1) - A[1], T(i1, 2) - A[2]};

        // Vector from point A to vertex i2
        double v2[3] = {T(i2, 0) - A[0], T(i2, 1) - A[1], T(i2, 2) - A[2]};

        double cross_product[3];
        crossProduct3D(v1, v2, cross_product);

        area[i] = 0.5 * dotProduct3D(B, cross_product);
    }

    // Check if signs are consistent
    if((area[0] >= -EPS) && (area[1] >= -EPS) && (area[2] >= -EPS)) {
        *has_solution = true;
    } else if((area[0] <= +EPS) && (area[1] <= +EPS) && (area[2] <= +EPS)) {
        *has_solution = true;
    } else {*has_solution = false;}

    if(!(*has_solution)) {return;}

    // Compute intersection point C via normalized weights
    double sum_area = area[0] + area[1] + area[2];

    area[0] /= sum_area;
    area[1] /= sum_area;
    area[2] /= sum_area;

    C[0] = 0.0; C[1] = 0.0; C[2] = 0.0;
    for(int i = 0; i < 3; i++) {
        C[0] += area[i] * T(i, 0);
        C[1] += area[i] * T(i, 1);
        C[2] += area[i] * T(i, 2);
    }
}
