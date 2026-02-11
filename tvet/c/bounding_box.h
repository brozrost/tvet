#pragma once
#include <stddef.h>

void boundingBox(
    const int (*faces)[3], size_t nof_faces,
    const double (*nodes)[3], size_t nof_nodes,
    const double s[3],
    double (*boxes)[4]
);
