#include <twobody/stumpff.h>
#include <twobody/math_utils.h>

#include <math.h>
#include <float.h>
#include <stdbool.h>

double stumpff_c0(double alpha, double s) {
    double z = alpha * s*s;
    double sqrtz = sqrt(fabs(z));

    if(zero(z))             // zero
        return 1.0;
    else if(zero(alpha))    // parabolic
        return 1.0;
    else if(alpha < 0.0)    // hyperbolic
        return cosh(sqrtz);
    else                    // elliptic
        return cos(sqrtz);

}

double stumpff_c1(double alpha, double s) {
    double z = alpha * s*s;
    double sqrtz = sqrt(fabs(z));

    if(zero(z))             // zero
        return 1.0; // XXX: 1.0 or 0.0?
    else if(zero(alpha))     // parabolic
        return s;
    else if(alpha < 0.0)    // hyperbolic
        return sinh(sqrtz) / sqrtz;
    else                    // elliptic
        return sin(sqrtz) / sqrtz;
}

double stumpff_c2(double alpha, double s) {
    double z = alpha * s*s;
    double sqrtz = sqrt(fabs(z));

    if(zero(z))             // zero
        return 1.0 / 2.0;
    else if(zero(alpha))    // parabolic
        return s*s / 2.0;
    else if(alpha < 0.0)    // hyperbolic
        return (cosh(sqrtz) - 1.0) / -z;
    else                    // elliptic
        return (1.0 - cos(sqrtz)) / z;
}

double stumpff_c3(double alpha, double s) {
    double z = alpha * s*s;
    double sqrtz = sqrt(fabs(z));

    if(zero(z))             // zero
        return 1.0 / 6.0;
    else if(zero(alpha))    // parabolic
        return s*s*s / 6.0;
    else if(alpha < 0.0)    // hyperbolic
        return (sinh(sqrtz) - sqrtz) / (-z * sqrtz);
    else                    // elliptic
        return (sqrtz - sin(sqrtz)) / (z * sqrtz);
}

double stumpff_series(int k, double z) {
    // c_k(z) = sum (-z)^i / (k + 2i)!

    int fac = 1;
    for(int i = 1; i <= k; ++i)
        fac = fac * i;

    double numer = 1.0, denom = fac;
    double c = numer / denom, sum = c;

    int max_steps = 30;
    for(int i = 1; i < max_steps && fabs(c) >= DBL_EPSILON; ++i) {
        numer *= -z;
        denom *= (k + 2*i) * (k + 2*i - 1);

        c = numer / denom;
        sum += c;
    }

    return sum;
}

void stumpff_quad(double z, double *cs) {
    double z_min = 0.1;
    int n = 0;

    for(n = 0; z*z >= z_min; ++n) // TODO: use pow and trunc
        z = z / 4;

    double c2 = stumpff_series(2, z);
    double c3 = stumpff_series(3, z);

    double c1 = 1.0 - z*c3;
    double c0 = 1.0 - z*c2;

    while(n--) {
        c3 = (c2 + c0*c3) / 4.0;
        c2 = c1*c1 / 2.0;
        c1 = c0*c1;
        c0 = 2.0 * c0*c0 - 1.0;
    }

    cs[0] = c0; cs[1] = c1; cs[2] = c2; cs[3] = c3;
}

void stumpff_fast(double z, double *cs) {
    *(vec4d*)cs = stumpff_simd(z);
}
