#pragma once
#include <stddef.h>

/// @brief For each triangle face, compute a 2D axis-aligned bounding box of the
///        triangle after projecting the 3D mesh onto a plane perpendicular to s.
/// @param faces Array of triangles (mesh faces).
/// @param nof_faces Number of triangles (mesh faces).
/// @param nodes Array of vertex coordinates. 
/// @param nof_nodes Number of verticies.
/// @param s Direction vector defining the projection plane.
/// @param boxes Ouput array, each face is assigned u_min, u_max, v_min, v_max
void boundingBox(
    const int (*faces)[3], size_t nof_faces,
    const double (*nodes)[3], size_t nof_nodes,
    const double s[3],
    double (*boxes)[4]
);
