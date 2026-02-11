#include "vector_math.h"
#include <math.h>

double dotProduct3D(const double x[3], const double y[3]) {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2]; 
}

void crossProduct3D(const double x[3], const double y[3], double out[3]) {
    out[0] = x[1]*y[2] - x[2]*y[1];
    out[1] = x[2]*y[0] - x[0]*y[2];
    out[2] = x[0]*y[1] - x[1]*y[0];
}

int normalizeVector3D(double v[3]) {
    double length = sqrt(dotProduct3D(v, v));
    if(length <= 0.0) {return 0;}

    v[0] /= length;
    v[1] /= length;
    v[2] /= length;

    return 1;
}