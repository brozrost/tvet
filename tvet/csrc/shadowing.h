#pragma once

#include <stddef.h>

/// @brief Computes directional cosine (illumination) for each face.
/// @param normals Face normals.
/// @param nof_faces Number of faces.
/// @param s 3D direction vector (Sun).
/// @param mu_i Output array.
void mu(const double (*normals)[3], size_t nof_faces, const double s[3], double *mu_i);

/// @brief Masks non-illuminated or non-visible faces.
/// @param mu_i Illumination cosine array. 
/// @param mu_e Emission/visibility cosine array.
/// @param nof_faces Number of faces.
/// @param nu_i Output illumination mask.
/// @param nu_e Output visibility mask.
void non(const double *mu_i, const double *mu_e, size_t nof_faces, double *nu_i, double *nu_e);

/// @brief Computes mutual shadowing for a non-convex mesh.
/// @param faces Triangle indicies.
/// @param nof_faces Number of faces.
/// @param nodes Vertex coordinates.
/// @param nof_nodes Number of nodes.
/// @param normals Face normals.
/// @param centers Face centers.
/// @param s 3D direction vector (Sun).
/// @param nu_i In/out illumination mask.
void nu(
    const int (*faces)[3], size_t nof_faces,
    const double (*nodes)[3], size_t nof_nodes,
    const double (*normals)[3],
    const double (*centers)[3],
    const double s[3],
    double *nu_i
);
