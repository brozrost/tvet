#pragma once
#include <stdbool.h>

/// @brief Computes the intersection of a line and a triangle in 3D space. 
/// @param A Point on the line in 3D space.
/// @param B Direction vector of the line.
/// @param t Triangle vertecies. 
/// @param C Output parameter, intersection point.
/// @param has_solution Boolean result. True if intersection exists.
void intersect_AB_t(
    const double A[3], 
    const double B[3], 
    const double t[3][3], 
    double C[3], 
    bool *has_solution
);
