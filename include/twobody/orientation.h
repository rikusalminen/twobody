#ifndef TWOBODY_ORIENTATION_H
#define TWOBODY_ORIENTATION_H

#ifndef TWOBODY_NO_SIMD
#include <twobody/simd4d.h>
#include <twobody/math_utils.h>
#include <math.h>
#include <float.h>

static inline vec4d orientation_major_axis(double i, double an, double arg)
    __attribute__((always_inline));
static inline vec4d orientation_major_axis(double i, double an, double arg) {
    return (vec4d) {
        (cos(arg) * cos(an)) - (sin(arg) * sin(an) * cos(i)),
        (sin(arg) * cos(an) * cos(i)) + (cos(arg) * sin(an)),
        sin(arg) * sin(i),
        0.0 };
}

static inline vec4d orientation_minor_axis(double i, double an, double arg)
    __attribute__((always_inline));
static inline vec4d orientation_minor_axis(double i, double an, double arg) {
    return (vec4d) {
        -(cos(arg) * sin(an) * cos(i)) - (sin(arg) * cos(an)),
        (cos(arg) * cos(an) * cos(i)) - (sin(arg) * sin(an)),
        cos(arg) * sin(i),
        0.0 };
}

static inline vec4d orientation_normal_axis(double i, double an, double arg)
    __attribute__((always_inline));
static inline vec4d orientation_normal_axis(double i, double an, double arg) {
    (void)arg;
    return (vec4d) {
        sin(an) * sin(i),
        -cos(an) * sin(i),
        cos(i),
        0.0 };
}

static inline double orientation_inclination(
    vec4d major,
    vec4d minor,
    vec4d normal)
    __attribute__((always_inline));
static inline double orientation_inclination(
    vec4d major,
    vec4d minor,
    vec4d normal) {
    (void)major; (void)minor;

    return acos(clamp(-1.0, 1.0, normal[2]));
}

static inline double orientation_longitude_of_ascending_node(
    vec4d major,
    vec4d minor,
    vec4d normal)
    __attribute__((always_inline));
static inline double orientation_longitude_of_ascending_node(
    vec4d major,
    vec4d minor,
    vec4d normal) {
    (void)major; (void)minor;

    // vector pointing to ascending node
    vec4d nodes = { -normal[1], normal[0], 0.0, 0.0 };
    if(dot(nodes, nodes) < DBL_EPSILON) // equatorial orbit
        return 0.0;

    return atan2(nodes[1], nodes[0]);
}

static inline double orientation_argument_of_periapsis(
    vec4d major,
    vec4d minor,
    vec4d normal)
    __attribute__((always_inline));
static inline double orientation_argument_of_periapsis(
    vec4d major,
    vec4d minor,
    vec4d normal) {
    (void)minor;

    vec4d nodes = { -normal[1], normal[0], 0.0, 0.0 };
    double N = dot(nodes, nodes);
    if(N < DBL_EPSILON) // equatorial orbit
        return sign(normal[2]) * atan2(major[1], major[0]);

    return sign(major[2]) *
        acos(clamp(-1.0, 1.0, dot(nodes, major) / sqrt(N)));
}

#endif

void orientation_major_axis_ptr(double *axis, double i, double an, double arg);
void orientation_minor_axis_ptr(double *axis, double i, double an, double arg);
void orientation_normal_axis_ptr(double *axis, double i, double an, double arg);

double orientation_inclination_ptr(
    const double *major,
    const double *minor,
    const double *normal);
double orientation_longitude_of_ascending_node_ptr(
    const double *major,
    const double *minor,
    const double *normal);
double orientation_argument_of_periapsis_ptr(
    const double *major,
    const double *minor,
    const double *normal);


#endif
