#ifndef TWOBODY_SIMD4D_H
#define TWOBODY_SIMD4D_H
#ifndef TWOBODY_NO_SIMD

#include <math.h>
#include <float.h>

typedef double vec4d __attribute__((vector_size(4 * sizeof(double))));

static inline vec4d splat4d(double x) __attribute__((always_inline));
static inline vec4d splat4d(double x) {
    return (vec4d){ x, x, x, x };
}

#ifdef __clang__
#define shuffle4d(x, y, a, b, c, d) (__builtin_shufflevector((x), (y), (a), (b), (c), (d)))
#else
typedef long v4sl __attribute__((vector_size(4 * sizeof(long))));
#define shuffle4d(x, y, a, b, c, d) \
    (__builtin_shuffle((x), (y), (v4sl){(a), (b), (c), (d)}))
#endif

#if __GNUC__ <= 4 && __GNUC_MINOR__ <= 6

static inline double dot(vec4d a, vec4d b) __attribute__((always_inline));
static inline double dot(vec4d a, vec4d b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
}

static inline vec4d dot4d(vec4d a, vec4d b) __attribute__((always_inline));
static inline vec4d dot4d(vec4d a, vec4d b) {
    return splat4d(dot(a, b));
}

static inline vec4d cross(vec4d a, vec4d b) __attribute__((always_inline));
static inline vec4d cross(vec4d a, vec4d b) {
    return (vec4d){
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
        0.0 };
}

static inline vec4d xyz4d(vec4d a) __attribute__((always_inline));
static inline vec4d xyz4d(vec4d a) {
    return (vec4d){ a[0], a[1], a[2], 0.0 };
}

#else // GCC > 4.7 or Clang required (__builtin_shuffle)

static inline vec4d dot4d(vec4d a, vec4d b) __attribute__((always_inline));
static inline vec4d dot4d(vec4d a, vec4d b) {
    vec4d prod = a * b;
    vec4d sum1 = prod + shuffle4d(prod, prod, 1, 0, 3, 2);
    vec4d sum2 = sum1 + shuffle4d(sum1, sum1, 2, 2, 0, 0);
    return sum2;
}

static inline double dot(vec4d a, vec4d b) __attribute__((always_inline));
static inline double dot(vec4d a, vec4d b) {
    return dot4d(a, b)[0];
}

static inline vec4d cross(vec4d a, vec4d b) __attribute__((always_inline));
static inline vec4d cross(vec4d a, vec4d b) {
    return shuffle4d(a, a, 1, 2, 0, 3) * shuffle4d(b, b, 2, 0, 1, 3)
        - shuffle4d(a, a, 2, 0, 1, 3) * shuffle4d(b, b, 1, 2, 0, 3);
}

static inline vec4d xyz4d(vec4d a) __attribute__((always_inline));
static inline vec4d xyz4d(vec4d a) {
    const vec4d zero = { 0.0, 0.0, 0.0, 0.0 };
    return shuffle4d(a, zero, 0, 1, 2, 7);
}

#endif

static inline double mag(vec4d a) __attribute__((always_inline));
static inline double mag(vec4d a) {
    return sqrt(dot(a, a));
}

static inline vec4d mag4d(vec4d a) __attribute__((always_inline));
static inline vec4d mag4d(vec4d a) {
    return splat4d(mag(a));
}

static inline vec4d unit4d(vec4d a) __attribute__((always_inline));
static inline vec4d unit4d(vec4d a) {
    return a / mag4d(a);
}

static inline int eqv4d(vec4d a, vec4d b) __attribute__((always_inline));
static inline int eqv4d(vec4d a, vec4d b) {
    double threshold = 1.0e-15; // DBL_EPSILON;
    return ((dot(a,a) < threshold && dot(b,b) < threshold) ||
        (dot(a-b, a-b)/(dot(a, a) + dot(b, b))) < threshold);
}

#endif
#endif
