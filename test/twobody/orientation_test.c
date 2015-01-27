#include <twobody/orientation.h>

#include "../numtest.h"

void orientation_test(
    double *params,
    int num_params,
    void *extra_args,
    struct numtest_ctx *test_ctx) {
    (void)extra_args;
    ASSERT(num_params == 3, "");

    double i = params[0] * M_PI;
    double an = (ZEROF(i) || ZEROF(i - M_PI))  ? 0.0 :
        (-1.0 + 2.0*params[1]) * M_PI;
    double arg = (-1.0 + 2.0*params[2]) * M_PI;

    vec4d major = orientation_major_axis(i, an, arg);
    vec4d minor = orientation_minor_axis(i, an, arg);
    vec4d normal = orientation_normal_axis(i, an, arg);

    ASSERT_EQF(mag(major), 1.0,
        "Major axis is an unit vector");
    ASSERT_EQF(mag(minor), 1.0,
        "Minor axis is an unit vector");
    ASSERT_EQF(mag(normal), 1.0,
        "Normal axis is an unit vector");

    ASSERT(ZEROF(major[3]) && ZEROF(minor[3]) && ZEROF(normal[3]),
        "w component is zero");

    ASSERT(ZEROF(dot(major, minor)) &&
        ZEROF(dot(major, normal)) &&
        ZEROF(dot(minor, normal)),
        "Axes are orthogonal");

    ASSERT(eqv4d(cross(major, minor), normal),
        "Axis cross product");

    double i2 = orientation_inclination(major, minor, normal);
    double an2 =
        orientation_longitude_of_ascending_node(major, minor, normal);
    double arg2 = orientation_argument_of_periapsis(major, minor, normal);

    ASSERT(isfinite(i2) && isfinite(an2) && isfinite(arg2),
        "Orientation elements not NaN");

    ASSERT_EQF(i, i2, "Inclination identity");
    ASSERT_EQF(an, an2, "Longitude of ascending node identity");
    ASSERT_EQF(arg, arg2, "Argument of periapsis identity");
}
