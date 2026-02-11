// g++ test_intersect_AB_t.c intersect_AB_t.c -o test_intersect_AB_t.out

#include "intersect_AB_t.h"
#include <assert.h>
#include <math.h>
#include <stdbool.h>

#define EPS_TEST 1e-12

static void test_intersect_AB_t_INTERSECT(void) {
    const double t[3][3] = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}
    };

    const double A[3] = {0.25, 0.25, 1.0};
    const double B[3] = {0.0,  0.0, -1.0};
    double C[3] = {0.0, 0.0, 0.0};
    bool has_solution = false;

    intersect_AB_t(A, B, t, C, &has_solution);

    assert(has_solution == true);
    assert(fabs(C[0] - 0.25) < EPS_TEST);
    assert(fabs(C[1] - 0.25) < EPS_TEST);
    assert(fabs(C[2] - 0.0) < EPS_TEST);
}

static void test_intersect_AB_t_NOINTERSECT(void) {
    const double t[3][3] = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}
    };

    const double A[3] = {2.0, 2.0, 1.0};
    const double B[3] = {0.0,  0.0, -1.0};
    double C[3] = {0.0, 0.0, 0.0};
    bool has_solution = false;

    intersect_AB_t(A, B, t, C, &has_solution);

    assert(has_solution == false);
}

static void run_tests(void) {
    test_intersect_AB_t_INTERSECT();
    test_intersect_AB_t_NOINTERSECT();
}

int main(void) {
    run_tests();
    return 0;
}
