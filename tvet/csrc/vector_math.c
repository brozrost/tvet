#include "vector_math.h"
#include <math.h>

double dotProduct3D(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; 
}

void crossProduct3D(const double a[3], const double b[3], double out[3]) {
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
}

int normalizeVector3D(double v[3]) {
    double length = sqrt(dotProduct3D(v, v));
    if(length <= 0.0) {return 0;}

    v[0] /= length;
    v[1] /= length;
    v[2] /= length;

    return 1;
}