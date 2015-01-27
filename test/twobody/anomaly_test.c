#include <twobody/conic.h>
#include <twobody/anomaly.h>
#include <twobody/math_utils.h>

#include <math.h>
#include <float.h>

#include "../numtest.h"

void anomaly_test(
    double *params,
    int num_params,
    void *extra_args,
    struct numtest_ctx *test_ctx) {
    (void)extra_args;

    ASSERT(num_params == 2, "num_params");

    double e = params[0] * 2.0;
    double t = -1.0 + params[1] * 2.0;

    double maxE = M_PI;
    double maxf = anomaly_eccentric_to_true(e, maxE);
    double maxM = anomaly_eccentric_to_mean(e, maxE);
    ASSERT(isfinite(maxf) && isfinite(maxM),
        "Eccentric, mean anomaly not NaN");

    double M = t * maxM;
    double E = t * maxE;
    double f = t * maxf;

    double fasymptote = conic_max_true_anomaly(e);

    double f2 = anomaly_eccentric_to_true(e, E);
    ASSERT(isfinite(f2), "True anomaly not NaN");
    ASSERT_RANGEF(f2, -maxf, maxf, "True anomaly within range");
    ASSERT_RANGEF(f2, -fasymptote, fasymptote, "True anomaly asymptote");
    ASSERT_EQF(E, anomaly_true_to_eccentric(e, f2), "Eccentric -> True");

    double E2 = anomaly_true_to_eccentric(e, f);
    ASSERT(isfinite(E2), "Eccentric anomaly not NaN");
    ASSERT_RANGEF(E2, -maxE, maxE, "Eccentric anomaly within range");
    ASSERT_EQF(f, anomaly_eccentric_to_true(e, E2), "True -> Eccentric");

    double E3 = anomaly_mean_to_eccentric(e, M);
    ASSERT(isfinite(E3), "Eccentric anomaly not NaN");
    ASSERT_RANGEF(E3, -maxE, maxE, "Eccentric anomaly within range");
    ASSERT_EQF(M, anomaly_eccentric_to_mean(e, E3), "True -> Eccentric");

    double M2 = anomaly_eccentric_to_mean(e, E);
    ASSERT(isfinite(M2), "Mean anomaly not NaN");
    ASSERT_RANGEF(M2, -maxM, maxM, "Mean anomaly within range");
    ASSERT_EQF(E, anomaly_mean_to_eccentric(e, M2), "Eccentric -> Mean");

    double M3 = anomaly_true_to_mean(e, f);
    ASSERT(isfinite(M3), "Mean anomaly not NaN");
    ASSERT_RANGEF(M3, -maxM, maxM, "Mean anomaly within range");
    ASSERT_EQF(f, anomaly_mean_to_true(e, M3), "True -> Mean");

    double f3 = anomaly_mean_to_true(e, M);
    ASSERT(isfinite(f3), "True anomaly not NaN");
    ASSERT_RANGEF(f3, -maxf, maxf, "True anomaly within range");
    ASSERT_EQF(M, anomaly_true_to_mean(e, f3), "Mean -> True");

    double dEdM = anomaly_dEdM(e, E);
    double MM = anomaly_eccentric_to_mean(e, E);
    double dM = 1.0e-9 * maxM;
    double Eplus = anomaly_mean_to_eccentric(e, MM+dM);
    double Eminus = anomaly_mean_to_eccentric(e, MM-dM);
    ASSERT(isfinite(dEdM), "dE/dM not NaN");
    ASSERT(dEdM > 0.0, "dE/dM is positive");
    ASSERT_EQF(dEdM * 2.0 * dM, (Eplus-Eminus), "dE/dM");

    double dE = (2.0*M_PI) * (1.0 / 36000.0); // 1/100 degree
    double dfdE = anomaly_dfdE(e, E);
    double fplus = anomaly_eccentric_to_true(e, E+dE);
    double fminus = anomaly_eccentric_to_true(e, E-dE);
    if(fplus < fminus)
        fplus += 2.0 * M_PI;
    ASSERT(isfinite(dfdE), "df/dE not NaN");
    ASSERT(dfdE > 0.0, "df/dE is positive");
    ASSERT_EQF(dfdE * 2.0 * dE, (fplus - fminus), "df/dE");

    double ff = anomaly_eccentric_to_true(e, E);
    double sinf = anomaly_true_sin(e, E);
    double cosf = anomaly_true_cos(e, E);
    ASSERT(isfinite(sinf) && isfinite(cosf),
        "True anomaly sine, cosine not NaN");
    ASSERT_EQF(sinf, sin(ff), "True anomaly sine");
    ASSERT_EQF(cosf, cos(ff), "True anomaly cosine");

    if(fabs(ff) < M_PI*0.9) { // tan(pi/2) goes to infinity
        double tanhalff = anomaly_true_tan_half(e, E);
        ASSERT(isfinite(tanhalff), "True anomaly tan half not NaN");
        ASSERT_EQF(tanhalff, tan(ff/2.0), "True anomaly tan half");
    }

    double EE = anomaly_true_to_eccentric(e, f);
    double sinE = anomaly_eccentric_sin(e, f);
    double cosE = anomaly_eccentric_cos(e, f);
    double tanhalfE = anomaly_eccentric_tan_half(e, f);

    if(!conic_parabolic(e)) {
        ASSERT(isfinite(sinE) && isfinite(cosE),
            "Eccentric anomaly sine, cosine not NaN");
    }

    if(conic_parabolic(e)) {
        // no sine and cosine for parabolic anomaly
    } else if(conic_hyperbolic(e)) {
        // hyperbola
        ASSERT_EQF(cosE, cosh(EE), "Hyperbolic anomaly cosine");
        ASSERT_EQF(sinE, sinh(EE), "Hyperbolic anomaly sine");

        ASSERT(isfinite(tanhalfE), "Hyperbolic anomaly tan half not NaN");
        ASSERT_EQF(tanhalfE, tanh(EE/2.0), "Hyperbolic anomaly tan half");
    } else {
        // ellipse
        ASSERT_EQF(cosE, cos(EE), "Eccentric anomaly cosine");
        ASSERT_EQF(sinE, sin(EE), "Eccentric anomaly sine");

        if(fabs(EE) < M_PI*0.9) { // tan(pi/2) goes to infinity
            ASSERT(isfinite(tanhalfE), "Eccentric anomaly tan half not NaN");
            ASSERT_EQF(tanhalfE, tan(EE/2.0), "Eccentric anomaly tan half");
        }
    }
}
