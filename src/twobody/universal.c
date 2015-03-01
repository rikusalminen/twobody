#include <twobody/conic.h>
#include <twobody/stumpff.h>
#include <twobody/universal.h>
#include <twobody/math_utils.h>

#include <math.h>

double universal_alpha(double mu, double r, double v2) {
    return 2.0/r - v2/mu;
}

double universal_period(double mu, double alpha) {
    if(alpha <= 0.0)
        return INFINITY;
    return 2.0*M_PI / sqrt(mu * alpha*alpha*alpha);
}

int universal_parabolic(double alpha) {
    return fabs(alpha) < DBL_EPSILON;
}

int universal_hyperbolic(double alpha) {
    return !universal_parabolic(alpha) && alpha < 0.0;
}

int universal_elliptic(double alpha) {
    return !universal_parabolic(alpha) && alpha > 0.0;
}

double universal_from_eccentric(double p, double e, double dE) {
    double a = conic_semi_major_axis(p, e);

    if(conic_parabolic(e))
        return sqrt(p) * dE;
    if(conic_hyperbolic(e))
        return sqrt(-a) * dE;
    else
        return sqrt(a) * dE;
}

double universal_half_sin_true(
    double p,
    double r0,
    double r,
    double s_half,
    const double *cs_half) {
    return sqrt(p / (r0*r)) * s_half * cs_half[1];
}

double universal_half_cos_true(
    double p,
    double r0, double sigma0,
    double r,
    double s_half,
    const double *cs_half) {
    (void)p;
    return 1.0/sqrt(r0*r) *
        (r0 * cs_half[0] + sigma0 * s_half * cs_half[1]);
}

double universal_half_tan_true(
    double p,
    double r0, double sigma0,
    double s_half,
    const double *cs_half) {
    return sqrt(p) * s_half * cs_half[1] /
        (r0 * cs_half[0] + sigma0 * s_half * cs_half[1]);
}

double universal_to_true(
    double p,
    double r0, double sigma0,
    double r,
    double s_half,
    const double *cs_half) {
    return 2.0 * atan2(
        universal_half_sin_true(p, r0, r, s_half, cs_half),
        universal_half_cos_true(p, r0, sigma0, r, s_half, cs_half));
}

double universal_time(
    double mu,
    double r0,
    double sigma0,
    double s,
    const double *cs) {
    return (r0 * s * cs[1] + sigma0 * s*s * cs[2] + s*s*s * cs[3]) / sqrt(mu);
}

double universal_radius(
    double r0,
    double sigma0,
    double s,
    const double *cs) {
    return r0 * cs[0] + sigma0 * s * cs[1] + s*s * cs[2];
}

double universal_sigma(
    double alpha,
    double r0, double sigma0,
    double s,
    const double *cs) {
    return sigma0 * cs[0] + (1.0 - alpha * r0) * s * cs[1];
}

double universal_guess_s(
    double mu,
    double alpha,
    double r0, double sigma0,
    double time) {
    if(universal_parabolic(alpha))
        return sqrt(mu) / r0  * time; // XXX: this sucks
    else if(universal_hyperbolic(alpha))
        return zero(time) ? 0.0 :
            sign(time) * sqrt(1.0/-alpha) *
            log(-2.0 * alpha * mu * time /
                (sqrt(mu)*sigma0 +
                     sign(time) * sqrt(mu/-alpha) * (1 - alpha*r0)));
    else
        return alpha*sqrt(mu) * time;
}

double universal_iterate_s(
    double mu,
    double alpha,
    double r0, double sigma0,
    double s0, double time,
    int max_steps) {
    if(max_steps <= 0)
        max_steps = 30;

    double threshold = DBL_EPSILON;
    double s = s0;

    for(int step = 0; step < max_steps; ++step) {
        double z = alpha * s*s;
        double cs[4];
        stumpff_fast(z, cs);

        double t = universal_time(mu, r0, sigma0, s, cs);
        double r = universal_radius(r0, sigma0, s, cs);
        double ds = sqrt(mu) * (time - t) / r;

        if(ds*ds < threshold)
            break;

        s = s + ds;
    }

    return s;
}
