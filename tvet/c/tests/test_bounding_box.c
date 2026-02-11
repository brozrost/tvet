// g++ test_bounding_box.c ../bounding_box.c ../vector_math.c -o out/test_bounding_box.out

#include "../bounding_box.h"
#include <math.h>
#include <assert.h>

#define EPS_TEST 1e-12

static void test_boundingBox() {
    const double nodes[4][3] = {
        {0.0, 0.0, 0.0},
        {2.0, 0.0, 0.0},
        {0.0, 3.0, 0.0},
        {-1.0, -2.0, 0.0}
    };

    const int faces[2][3] = {
        {0, 1, 2},
        {0, 2, 3} 
    };

    const double s[3] = {0.0, 0.0, 1.0};
    double boxes[2][4] = {0};

    boundingBox(faces, 2, nodes, 4, s, boxes);

    // Face 1:
    assert(fabs(boxes[0][0] - 0.0) < EPS_TEST);
    assert(fabs(boxes[0][1] - 2.0) < EPS_TEST);
    assert(fabs(boxes[0][2] - 0.0) < EPS_TEST);
    assert(fabs(boxes[0][3] - 3.0) < EPS_TEST);

    // Face 2:
    assert(fabs(boxes[1][0] + 1.0) < EPS_TEST);
    assert(fabs(boxes[1][1] - 0.0) < EPS_TEST);
    assert(fabs(boxes[1][2] + 2.0) < EPS_TEST);
    assert(fabs(boxes[1][3] - 3.0) < EPS_TEST);
}

static void run_tests(void) {
    test_boundingBox();
}

int main(void) {
    run_tests();
    return 0;
}