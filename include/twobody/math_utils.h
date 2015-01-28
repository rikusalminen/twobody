#ifndef TWOBODY_MATH_UTILS_H
#define TWOBODY_MATH_UTILS_H

#include <float.h>

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

#endif
