// g++ test_shadowing.c ../shadowing.c ../bounding_box.c ../intersect_AB_t.c ../vector_math.c -o out/test_shadowing.out

#include "../shadowing.h"
#include <math.h>
#include <assert.h>

#define EPS_TEST 1e-12

static void computeCenter(const double nodes[][3], const int face[3], double center[3]) {
    const double *p0 = nodes[face[0]];
    const double *p1 = nodes[face[1]];
    const double *p2 = nodes[face[2]];

    center[0] = (p0[0] + p1[0] + p2[0]) / 3.0;
    center[1] = (p0[1] + p1[1] + p2[1]) / 3.0;
    center[2] = (p0[2] + p1[2] + p2[2]) / 3.0;
}

static void test_shadowing(void) {
    const double nodes[6][3] = {
        // Face 0:
        {-1.0, -1.0, 0.0},
        { 1.0, -1.0, 0.0},
        { 0.0,  1.0, 0.0},

        // Face 1:
        {-0.5, -0.5, 0.5},
        { 0.5, -0.5, 0.5},
        { 0.0,  0.5, 0.5}
    };

    const int faces[2][3] = {
        {0, 1, 2},
        {3, 4, 5}
    };

    const double s[3] = {0.0, 0.0, 1.0};

    double centers[2][3];
    double normals[2][3];

    for(int i = 0; i < 2; i++) {
        computeCenter(nodes, faces[i], centers[i]);
        normals[i][0] = 0.0;
        normals[i][1] = 0.0;
        normals[i][2] = 1.0;
    }

    double mu_i[2] = {0};
    double mu_e[2] = {0};

    mu(normals, 2, s, mu_i);
    mu_e[0] = mu_i[0]; mu_e[1] = mu_i[1];

    double nu_i[2] = {0};
    double nu_e[2] = {0};

    non(mu_i, mu_e, 2, nu_i, nu_e);

    // Both initially visible
    assert(fabs(nu_i[0] - 1.0) < EPS_TEST);
    assert(fabs(nu_i[1] - 1.0) < EPS_TEST);

    nu(
        faces, 2,
        nodes, 6,
        normals,
        centers,
        s,
        nu_i
    );

    // Face 0 should be shadowed, face 1 remain visible
    assert(fabs(nu_i[0] - 0.0) < EPS_TEST);
    assert(fabs(nu_i[1] - 1.0) < EPS_TEST);
}

static void run_tests(void) {
    test_shadowing();
}

int main(void) {
    run_tests();
    return 0;
}
