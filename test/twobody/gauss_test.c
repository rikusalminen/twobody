#include <math.h>
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

    double E1 = (-1.0 + params[3] * 1.5) * M_PI;
    double f1 = anomaly_eccentric_to_true(e, E1);
    double M1 = anomaly_eccentric_to_mean(e, E1);
    double t1 = t0 + M1 / n;
    double r1 = eccentric_radius(p, e, E1);

    double E2 = e < 1.0 ?
        E1 + (0.1 + params[4] * 0.9) * M_PI :
        E1 + (M_PI - E1) * (0.1 + params[4] * 0.9);
    double f2 = anomaly_eccentric_to_true(e, E2);
    double M2 = anomaly_eccentric_to_mean(e, E2);
    double t2 = t0 + M2 / n;
    double r2 = eccentric_radius(p, e, E2);

    double z0 = f2-f1;
    double z = gauss_iterate_z(mu, r1, r2, f2-f1, t2-t1, z0, 30);

    ASSERT(isfinite(z), "z not NaN");
    ASSERT_EQF(z, (E2-E1)*(E2-E1),
        "z = dE^2");
}
