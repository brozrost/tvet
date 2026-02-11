#include "vector_math.h"
#include <math.h>

/// @brief Computes dot product of two 3D vectors.
/// @param x 3D vector.
/// @param y 3D vector.
/// @return dot product.
double dotProduct3D(const double x[3], const double y[3]) {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2]; 
}

/// @brief Computes cross product of two 3D vectors.
/// @param x 3D vector.
/// @param y 3D vector.
/// @param out Ouput parameter, resulting vector.
void crossProduct3D(const double x[3], const double y[3], double out[3]) {
    out[0] = x[1]*y[2] - x[2]*y[1];
    out[1] = x[2]*y[0] - x[0]*y[2];
    out[2] = x[0]*y[1] - x[1]*y[0];
}

/// @brief Normalizes a 3D vector.
/// @param v 3D vector.
/// @return 0 for vectors of length <= 0, 1 for all others.
int normalizeVector3D(double v[3]) {
    double length = sqrt(dotProduct3D(v, v));
    if(length <= 0.0) {return 0;}

    v[0] /= length;
    v[1] /= length;
    v[2] /= length;

    return 1;
}