#ifndef TWOBODY_MATH_UTILS_H
#define TWOBODY_MATH_UTILS_H

#include <float.h>
#include <math.h>

static inline double clamp(double min, double max, double x) {
    return (x < min ? min : (x > max ? max : x));
}

static inline int zero(double x) {
    return x*x < DBL_EPSILON;
}

static inline double sign(double x) {
    return x < 0 ? -1.0 : 1.0;
}

static inline double square(double x) {
    return x*x;
}

static inline double cube(double x) {
    return x*x*x;
}

static inline double angle_clamp(double x0) {
    double x = (x0+M_PI)/(2.0*M_PI);
    return -M_PI + 2.0*M_PI * (x - floor(x));
}

#endif
