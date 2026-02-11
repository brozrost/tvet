#pragma once

#include <stddef.h>

void mu(const double *normals, size_t nof_faces, const double s[3], double *mu_i);
void non(const double *mu_i, const double *mu_e, size_t nof_faces, double *nu_i, double *nu_e);
void nu(
    const int (*faces)[3], size_t nof_faces,
    const double (*nodes)[3], size_t nof_nodes,
    const double (*normals)[3],
    const double (*centers)[3],
    const double s[3],
    double *nu_i
);
