#include "vector_math.h"

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
