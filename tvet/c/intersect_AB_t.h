#include <stdbool.h>

void intersect_AB_t(
    const double A[3], 
    const double B[3], 
    const double t[3][3], 
    double C[3], 
    bool *has_solution
);
