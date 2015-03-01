#include <twobody/conic.h>
#include <twobody/anomaly.h>
#include <twobody/true_anomaly.h>
#include <twobody/eccentric_anomaly.h>
#include <twobody/stumpff.h>
#include <twobody/universal.h>
#include <twobody/fg.h>
#include <twobody/simd4d.h>

#include "../numtest.h"

void fg_test(
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

    double E1 = (-1.0 + params[3] * 2.0) * M_PI;
    double f1 = anomaly_eccentric_to_true(e, E1);
    double M1 = anomaly_eccentric_to_mean(e, E1);
    double t1 = t0 + M1 / n;
    double r1 = eccentric_radius(p, e, E1);

    double E2 = (-1.0 + params[4] * 2.0) * M_PI;
    double f2 = anomaly_eccentric_to_true(e, E2);
    double M2 = anomaly_eccentric_to_mean(e, E2);
    double t2 = t0 + M2 / n;
    double r2 = eccentric_radius(p, e, E2);

    double x1 = eccentric_x(p, e, E1), y1 = eccentric_y(p, e, E1);
    double xdot1 = eccentric_xdot(mu, p, e, E1), ydot1 = eccentric_ydot(mu, p, e, E1);
    double x2 = eccentric_x(p, e, E2), y2 = eccentric_y(p, e, E2);
    double xdot2 = eccentric_xdot(mu, p, e, E2), ydot2 = eccentric_ydot(mu, p, e, E2);

    double v1 = eccentric_velocity(mu, p, e, E1);
    double sigma1 = (x1*xdot1 + y1*ydot1) / sqrt(mu);

    vec4d pos1 = { x1, y1, 0.0, 0.0 }, vel1 = { xdot1, ydot1, 0.0, 0.0 };
    vec4d pos2 = { x2, y2, 0.0, 0.0 }, vel2 = { xdot2, ydot2, 0.0, 0.0 };

    double f, g, fdot, gdot;
    fg(x1, y1, xdot1, ydot1,
        x2, y2, xdot2, ydot2,
        &f, &g, &fdot, &gdot);

    ASSERT(isfinite(f) && isfinite(g) && isfinite(fdot) && isfinite(gdot),
        "fg not NaN");
    ASSERT_EQF(f*gdot, 1.0 + fdot*g,
        "fg identity (sanity)");

    ASSERT(eqv4d(splat4d(f) * pos1 + splat4d(g) * vel1, pos2),
        "fg position identity");
    ASSERT(eqv4d(splat4d(fdot) * pos1 + splat4d(gdot) * vel1, vel2),
        "fg velocity identity");

    // true anomaly
    double ff = true_f(mu, p, r1, r2, f2-f1);
    double gf = true_g(mu, p, r1, r2, f2-f1);
    double fdotf = true_fdot(mu, p, r1, r2, f2-f1);
    double gdotf = true_gdot(mu, p, r1, r2, f2-f1);

    ASSERT(isfinite(ff) && isfinite(gf) && isfinite(fdotf) && isfinite(gdotf),
        "fg not NaN (true anomaly)");
    ASSERT_EQF(ff*gdotf, 1.0 + fdotf*gf,
        "fg identity (true anomaly)");

    vec4d posf = splat4d(ff) * pos1 + splat4d(gf) * vel1;
    vec4d velf = splat4d(fdotf) * pos1 + splat4d(gdotf) * vel1;
    ASSERT(eqv4d(posf, pos2),
        "fg position identity (true anomaly)");
    if(fabs(f2-f1) < M_PI) // XXX: accuracy is bad at 180 degrees
        ASSERT(eqv4d(velf, vel2),
            "fg velocity identity (true anomaly)");

    // eccentric anomaly
    double fE = eccentric_f(mu, p, e, r1, E2-E1);
    double gE = eccentric_g(mu, p, e, r1, sigma1, E2-E1);
    double gE_t = eccentric_g_t(mu, p, e, E2-E1, t2-t1);
    double fdotE = eccentric_fdot(mu, p, e, r1, r2, E2-E1);
    double gdotE = eccentric_gdot(mu, p, e, r2, E2-E1);

    ASSERT(isfinite(fE) && isfinite(gE) && isfinite(fdotE) && isfinite(gdotE),
        "fg not NaN (eccentric anomaly)");
    ASSERT(isfinite(gE_t),
        "fg g time function not NaN (eccentric anomaly)");
    ASSERT_EQF(fE*gdotE, 1.0 + fdotE*gE,
        "fg identity (eccentric anomaly)");

    vec4d posE = splat4d(fE) * pos1 + splat4d(gE) * vel1;
    vec4d velE = splat4d(fdotE) * pos1 + splat4d(gdotE) * vel1;

    ASSERT(eqv4d(posE, pos2),
        "fg position identity (eccentric anomaly)");
    ASSERT(eqv4d(velE, vel2),
        "fg velocity identity (eccentric anomaly)");

    ASSERT_EQF(fE*gdotE, 1.0 + fdotE*gE_t,
        "fg identity with time (eccentric anomaly)");
    ASSERT_EQF(t2-t1, -eccentric_g_t(mu, p, e, E2-E1, -gE),
        "fg g function time identity (eccentric anomaly)");

    // universal variables
    double alpha = universal_alpha(mu, r1, v1*v1);
    double s = universal_from_eccentric(p, e, E2-E1);
    double z = alpha * s*s;
    double cs[4];
    stumpff_fast(z, cs);

    double fs = universal_f(mu, r1, s, cs);
    double gs = universal_g(mu, r1, sigma1, s, cs);
    double gs_t = universal_g_t(mu, t2-t1, s, cs);
    double fdots = universal_fdot(mu, r1, r2, s, cs);
    double gdots = universal_gdot(mu, r2, s, cs);

    ASSERT(isfinite(fs) && isfinite(gs) && isfinite(fdots) && isfinite(gdots),
        "fg not NaN (universal)");
    ASSERT(isfinite(gs_t),
        "fg g time function not NaN (universal)");
    ASSERT_EQF(fs*gdots, 1.0 + fdots*gs,
        "fg identity (universal)");

    vec4d poss = splat4d(fs) * pos1 + splat4d(gs) * vel1;
    vec4d vels = splat4d(fdots) * pos1 + splat4d(gdots) * vel1;

    ASSERT(eqv4d(poss, pos2),
        "fg position identity (universal)");
    ASSERT(eqv4d(vels, vel2),
        "fg velocity identity (universal)");

    ASSERT_EQF(fs*gdots, 1.0 + fdots*gs_t,
        "fg identity with time (universal)");
    ASSERT_EQF(t2-t1, -universal_g_t(mu, -gs, s, cs),
        "fg g function time identity (universal)");
}
