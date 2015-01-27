#include <twobody/conic.h>

#include <float.h>
#include <math.h>

#include "../numtest.h"

void conic_test(
    double *params,
    int num_params,
    void *extra_args,
    struct numtest_ctx *test_ctx) {
    (void)extra_args;
    ASSERT(num_params == 3, "");

    double mu = 1.0 + params[0] * 1.0e10;
    double p = 1.0 + params[1] * 1.0e10;
    double e = params[2] * 2.0;

    ASSERT(conic_elliptic(e) + conic_parabolic(e) + conic_hyperbolic(e) == 1,
        "Conic is elliptic, parabolic or hyperbolic");
    ASSERT(!conic_circular(e) || conic_elliptic(e),
        "Circles are ellipses");
    ASSERT(conic_closed(e) == conic_elliptic(e),
        "Only elliptic orbits are closed");

    double h = conic_specific_angular_momentum(mu, p, e);
    ASSERT(isfinite(h) && h > 0.0,
        "Specific angular momentum is finite and positive");

    double visviva = conic_specific_orbital_energy(mu, p, e);
    ASSERT(isfinite(visviva),
        "Specific orbital energy not NaN");

    ASSERT_EQF(sqrt(fmax(0.0, 1.0 + (2.0*visviva*h*h / (mu*mu)))), e,
        "Eccentricity");

    double a = conic_semi_major_axis(p, e);
    double b = conic_semi_minor_axis(p, e);
    double c = conic_focal_distance(p, e);

    if(conic_parabolic(e)) {
        ASSERT(!isfinite(a) && !isfinite(b) && !isfinite(c),
            "Parabolic semi-major axis, semi-minor axis and focal distance "
            "are infinite");
    } else {
        ASSERT(isfinite(a) && isfinite(b) && isfinite(c),
            "Semi-major axis, semi-minor axis and focal distance not NaN");
    }

    if(conic_parabolic(e)) {
        ASSERT(ZEROF(visviva),
            "Parabola orbital energy is zero");
    } else if(conic_hyperbolic(e)) {
        ASSERT(a < 0.0 && b > 0.0,
            "Hyperbola major axis is negative, minor axis is positive");
        ASSERT_EQF(a*a + b*b, c*c,
            "Hyperbola focal distance");

        ASSERT(visviva > 0.0,
            "Hyperbola orbital energy is positive");
    } else if(conic_elliptic(e)) {
        ASSERT(a > 0.0 && b > 0.0,
            "Ellipse major and minor axis are positive");
        ASSERT_LTF(b, a,
            "Ellipse semi-minor axis is smaller than semi-major axis");
        ASSERT_EQF(a*a - b*b, c*c,
            "Ellipse focal distance");

        ASSERT(visviva < 0.0,
            "Ellipse orbital energy is negative");
    }

    double n = conic_mean_motion(mu, p, e);
    ASSERT(isfinite(n) && n > 0.0,
        "Mean motion is finite and positive");

    double rp = conic_periapsis(p, e);
    double vp = conic_periapsis_velocity(mu, p, e);
    double h_p = rp*vp, visviva_p = vp*vp/2.0 - mu/rp;
    ASSERT(isfinite(rp) && rp > 0.0 &&
        isfinite(vp) && vp > 0.0,
        "Periapsis radius and velocity are finite and positive");
    ASSERT_EQF(h_p, h,
        "Periapsis angular momentum");
    if(conic_parabolic(e)) { // accuracy is bad for parabolic trajectory
        ASSERT(ZEROF(visviva_p*visviva_p),
            "Periapsis orbital energy (parabolic)");
    } else {
        ASSERT_EQF(visviva_p, visviva,
            "Periapsis orbital energy");
    }

    double maxf = conic_max_true_anomaly(e);
    if(conic_hyperbolic(e)) {
        ASSERT_EQF(tan(maxf), b/a,
            "True anomaly asymptote (hyperbola)");
    } else {
        ASSERT_EQF(maxf, M_PI,
            "True anomaly asymptote");
    }

    double ra = conic_apoapsis(p, e);
    double va = conic_apoapsis_velocity(mu, p, e);

    if(conic_closed(e)) {
        ASSERT(isfinite(ra) && isfinite(va),
            "Closed orbit apoapsis radius and velocity are finite");
        ASSERT_LTF(rp, ra,
            "Periapsis radius is less than apoapsis radius");
        ASSERT_LTF(va, vp,
            "Periapsis velocity is greater than apoapsis radius");

        double h_a = ra*va, visviva_a = va*va/2.0 - mu/ra;
        ASSERT_EQF(h_a, h,
            "Apoapsis angular momentum");
        ASSERT_EQF(visviva_a, visviva,
            "Apoapsis orbital energy");

        double P = conic_period(mu, p, e);
        ASSERT_EQF(P*n, 2.0*M_PI,
            "Period and mean motion");
    } else {
        ASSERT(!isfinite(ra) && !isfinite(va),
            "Open orbit apoapsis radius and velocity are infinite");
    }
}
