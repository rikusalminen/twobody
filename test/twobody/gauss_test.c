#include <math.h>
#include <twobody/math_utils.h>
#include <twobody/conic.h>
#include <twobody/anomaly.h>
#include <twobody/true_anomaly.h>
#include <twobody/eccentric_anomaly.h>
#include <twobody/gauss.h>

#include "../numtest.h"

void gauss_test(
    double *params,
    int num_params,
    void *extra_args,
    struct numtest_ctx *test_ctx) {
    (void)extra_args;
    ASSERT(num_params == 5, "");

    double mu = 1.0 + params[0] * 1.0e5;
    double p = 1.0 + params[1] * 1.0e5;
    double e = params[2] * 2.0;

    double n = conic_mean_motion(mu, p, e);
    double t0 = 0.0;

    double E1, E2;
    if(e < 1.0) {
        E1 = (-1.0 + params[3] * 2.0) * M_PI;
        E2 = E1 + (0.1 + 0.8 * params[4]) * M_PI;
    } else {
        E1 = (-1.0 + params[3] * 1.5) * M_PI/2.0;
        E2 = E1 + (M_PI/2.0 - E1) * (0.1 + params[4] * 0.9);
    }

    double f1 = anomaly_eccentric_to_true(e, E1);
    double M1 = anomaly_eccentric_to_mean(e, E1);
    double t1 = t0 + M1 / n;
    double r1 = eccentric_radius(p, e, E1);

    double f2 = anomaly_eccentric_to_true(e, E2);
    double M2 = anomaly_eccentric_to_mean(e, E2);
    double t2 = t0 + M2 / n;
    double r2 = eccentric_radius(p, e, E2);

    double zE = conic_parabolic(e) ? 0.0 :
        (E2-E1)*(E2-E1) * (e < 1.0 ? 1.0 : -1.0);
    double z0 = square(2.0*M_PI * trunc((E2-E1) / (2.0*M_PI))) + square(f2-f1);
    double f, g, fdot, gdot;
    double z = gauss_iterate_z(
        mu, r1, r2, f2-f1, t2-t1, z0,
        &f, &g, &fdot, &gdot,
        30);

    ASSERT(isfinite(z), "z not NaN");
    ASSERT_EQF(z, zE,
        "z = dE^2");

    double sigma1 = r1 * eccentric_velocity_radial(mu, p, e, E1) / sqrt(mu);
    double fE = eccentric_f(mu, p, e, r1, E2-E1);
    double gE = eccentric_g(mu, p, e, r1, sigma1, E2-E1);
    //double fdotE = eccentric_fdot(mu, p, e, r1, r2, E2-E1);
    double gdotE = eccentric_gdot(mu, p, e, r2, E2-E1);

    ASSERT_EQF(f, fE, "Lagrangian coefficient f is equal");
    ASSERT_EQF(g, gE, "Lagrangian coefficient g is equal");
    //ASSERT_EQF(fdot, fdotE, "Lagrangian coefficient fdot is equal");
    ASSERT_EQF(gdot, gdotE, "Lagrangian coefficient gdot is equal");
}
