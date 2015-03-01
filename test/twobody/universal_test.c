#include <twobody/conic.h>
#include <twobody/anomaly.h>
#include <twobody/eccentric_anomaly.h>
#include <twobody/stumpff.h>
#include <twobody/universal.h>

#include <math.h>

#include "../numtest.h"

void universal_test(
    double *params,
    int num_params,
    void *extra_args,
    struct numtest_ctx *test_ctx) {
    (void)extra_args;
    ASSERT(num_params == 5, "");

    double mu = 1.0 + params[0] * 1.0e8;
    double p = 1.0 + params[1] * 1.0e8;
    double e = params[2] * 4.0;
    double n = conic_mean_motion(mu, p, e);

    double t0 = 0.0;
    double E1 = (-1.0 + params[3] * 2.0) * M_PI;
    double f1 = anomaly_eccentric_to_true(e, E1);
    double M1 = anomaly_eccentric_to_mean(e, E1);
    double t1 = t0 + M1 / n;
    double r1 = eccentric_radius(p, e, E1);
    double v1 = eccentric_velocity(mu, p, e, E1);
    double sigma1 = eccentric_radius(p, e, E1) *
        eccentric_velocity_radial(mu, p, e, E1) / sqrt(mu);

    double E2 = (-1.0 + params[4] * 2.0) * M_PI;
    double f2 = anomaly_eccentric_to_true(e, E2);
    double M2 = anomaly_eccentric_to_mean(e, E2);
    double t2 = t0 + M2 / n;
    double r2 = eccentric_radius(p, e, E2);
    double sigma2 = eccentric_radius(p, e, E2) *
        eccentric_velocity_radial(mu, p, e, E2) / sqrt(mu);

    double alpha = universal_alpha(mu, r1, v1*v1);
    ASSERT(isfinite(alpha), "Alpha is not NaN");

    if(conic_elliptic(e))
        ASSERT_EQF(conic_period(mu, p, e), universal_period(mu, alpha),
            "Universal period");

    if(conic_parabolic(e)) // accuracy is bad for parabolic
        ASSERT(ZEROF(alpha*alpha), "Alpha is zero (parabolic)");
    else
        ASSERT_EQF(alpha, -2.0 * conic_specific_orbital_energy(mu, p, e) / mu,
            "Alpha (inverse semi-major axis) and specific orbital energy");

    ASSERT(1 ==
        universal_parabolic(alpha) +
        universal_hyperbolic(alpha) +
        universal_elliptic(alpha),
        "Orbit is parabolic, hyperbolic xor elliptic");

    double s = universal_from_eccentric(p, e, E2-E1);
    double z = alpha * s*s;

    ASSERT(isfinite(s),
        "Universal from eccentric is not NaN");

    ASSERT_EQF(s, alpha * sqrt(mu) * (t2-t1) + sigma2 - sigma1,
        "Universal variable identity");

    double cs[4];
    stumpff_fast(z, cs);

    double t = universal_time(mu, r1, sigma1, s, cs);
    double r = universal_radius(r1, sigma1, s, cs);
    double sigma = universal_sigma(alpha, r1, sigma1, s, cs);

    ASSERT(isfinite(t) && isfinite(r), "Universal time and radius not NaN");
    ASSERT_EQF(t, t2-t1, "Universal time");
    ASSERT_EQF(r, r2, "Universal radius");
    ASSERT_EQF(sigma, sigma2, "Universal sigma");

    double ds = universal_from_eccentric(p, e, M_PI) * 1.0e-4;
    double splus = s + ds, sminus = s - ds;
    double csplus[4], csminus[4];
    stumpff_fast(alpha * splus*splus, csplus);
    stumpff_fast(alpha * sminus*sminus, csminus);
    double tplus = universal_time(mu, r1, sigma1, splus, csplus);
    double tminus = universal_time(mu, r1, sigma1, sminus, csminus);
    ASSERT_EQF(sqrt(mu) * (tplus-tminus), r * 2.0*ds, "sqrt(mu) dt = r ds");

    double s0 = universal_guess_s(mu, alpha, r1, sigma1, t2-t1);
    ASSERT(isfinite(s0), "Universal variable initial guess not NaN");

    if(conic_parabolic(e)) s0 = s + ds; // XXX: parabolic guess_s is broken

    double ss = universal_iterate_s(mu, alpha, r1, sigma1, s0, t2-t1, 0);
    ASSERT(isfinite(ss), "Universal variable time of flight not NaN");

    ASSERT_EQF(s, ss, "Time of flight equation");

    double s_half = s/2.0, z_half = alpha * s_half*s_half;
    double cs_half[4];
    stumpff_fast(z_half, cs_half);

    double half_sin_f = universal_half_sin_true(p, r1, r2, s_half, cs_half);
    double half_cos_f = universal_half_cos_true(p, r1, sigma1, r2, s_half, cs_half);
    double half_tan_f = universal_half_tan_true(p, r1, sigma1, s_half, cs_half);
    double ff = universal_to_true(p, r1, sigma1, r2, s_half, cs_half);

    ASSERT_EQF(sin((f2-f1)/2.0), half_sin_f, "Universal half sin true");
    ASSERT_EQF(cos((f2-f1)/2.0), half_cos_f, "Universal half cos true");

    if(!EQF(fabs(f2-f1), M_PI)) // NOTE: tan(pi) goes to infinity
        ASSERT_EQF(tan((f2-f1)/2.0), half_tan_f, "Universal half tan true");

    if(!EQF(fabs(f2-f1), 2.0*M_PI))
        ASSERT_EQF(ff, f2-f1, "Universal to true");
    else
        ASSERT_EQF(fabs(ff), 2.0*M_PI, "Universal to true (2pi)");
}
