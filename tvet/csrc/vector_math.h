#pragma once

/// @brief Computes dot product of two 3D vectors.
/// @param a 3D vector.
/// @param b 3D vector.
/// @return dot product.
double dotProduct3D(const double a[3], const double b[3]);

/// @brief Computes cross product of two 3D vectors.
/// @param a 3D vector.
/// @param b 3D vector.
/// @param out Ouput parameter, resulting vector.
void crossProduct3D(const double a[3], const double b[3], double out[3]);

/// @brief Normalizes a 3D vector.
/// @param v 3D vector.
/// @return 0 for vectors of length <= 0, 1 for all others.
int normalizeVector3D(double v[3]);
