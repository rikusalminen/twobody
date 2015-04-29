#ifndef TWOBODY_STUMPFF_H
#define TWOBODY_STUMPFF_H

double stumpff_c0(double alpha, double s);
double stumpff_c1(double alpha, double s);
double stumpff_c2(double alpha, double s);
double stumpff_c3(double alpha, double s);
double stumpff_series(int k, double z);

void stumpff_quad(double z, double *cs);
void stumpff_fast(double z, double *cs);

#ifndef TWOBODY_NO_SIMD
#include <twobody/simd4d.h>

static inline vec4d stumpff_simd(double z) __attribute__((always_inline));
static inline vec4d stumpff_simd(double z) {
    vec4d numer = { 1.0, 1.0, 1.0, 1.0 };
    vec4d denom = { 1.0, 1.0, 2.0, 6.0 };
    vec4d k = { 0.0, 1.0, 2.0, 3.0 };
    vec4d c = numer / denom;
    vec4d sum = c;
    vec4d is = { 1.0, 1.0, 1.0, 1.0 };

    int max_steps = 25;
    for(int i = 1; i < max_steps; ++i) {
        numer *= splat4d(-z);
        denom *= k + is;
        is += splat4d(1.0);
        denom *= k + is;
        is += splat4d(1.0);
        c = numer / denom;
        sum += c;
    }

    return sum;
}
#endif

#endif
