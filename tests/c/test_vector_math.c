// g++ test_vector_math.c ../../tvet/csrc/vector_math.c -o out/test_vector_math.out

#include "../../tvet/csrc/vector_math.h"
#include <math.h>
#include <assert.h>

#define EPS_TEST 1e-12

static void test_dotProduct3D(void) {
    const double a[3] = {1.0, 0.0, 0.0};
    const double b[3] = {0.0, 1.0, 0.0};
    assert(fabs(dotProduct3D(a, b)) < EPS_TEST);

    const double c[3] = {1.0, 2.0, 3.0};
    const double d[3] = {4.0, 5.0, 6.0};
    assert(fabs(dotProduct3D(c, d) - 32.0) < EPS_TEST);
    assert(fabs(dotProduct3D(c, c) - 14.0) < EPS_TEST);
}

static void test_crossProduct3D(void) {
    const double a[3] = {1.0, 0.0, 0.0};
    const double b[3] = {0.0, 1.0, 0.0};
    double out[3];

    crossProduct3D(a, b, out);

    assert(fabs(out[0] - 0.0) < EPS_TEST);
    assert(fabs(out[1] - 0.0) < EPS_TEST);
    assert(fabs(out[2] - 1.0) < EPS_TEST);

    assert(fabs(dotProduct3D(a, out)) < EPS_TEST);
    assert(fabs(dotProduct3D(b, out)) < EPS_TEST);
}

static void run_tests(void) {
    test_dotProduct3D();
    test_crossProduct3D();
}

int main(void) {
    run_tests();
    return 0;
}
