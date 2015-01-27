#include <twobody/conic.h>
#include <twobody/anomaly.h>
#include <twobody/true_anomaly.h>
#include <twobody/math_utils.h>

#include <math.h>
#include <float.h>

#include "../numtest.h"

void true_anomaly_test(
    double *params,
    int num_params,
    void *extra_args,
    struct numtest_ctx *test_ctx) {
    (void)extra_args;

    ASSERT(num_params == 4, "");

    double mu = 1.0 + params[0] * 1.0e10;
    double p = 1.0 + params[1] * 1.0e10;
    double e = params[2] * 4.0;

    double maxf = anomaly_eccentric_to_true(e, M_PI);
    double f = (-1.0 + params[3] * 2.0) * maxf;

    double a = conic_semi_major_axis(p, e);
    double b = conic_semi_minor_axis(p, e);
    double c = conic_focal_distance(p, e);
    double q = conic_periapsis(p, e);

    double dfdt = true_dfdt(mu, p, e, f);
    ASSERT(isfinite(dfdt), "df/dt not NaN");
    ASSERT(dfdt > 0.0, "df/dt is positive");

    double n = conic_mean_motion(mu, p, e);
    double dt = (2.0 * M_PI / n) * (1.0 / 36000.0); // 1/100 degree
    double df = dfdt * dt;

    double r = true_radius(p, e, f);
    ASSERT(isfinite(r), "Radius not NaN");
    ASSERT(r > 0, "Radius is positive");

    ASSERT_LTF(q, r,
        "Radius is larger than periapsis");
    if(conic_closed(e))
        ASSERT_LTF(r, conic_apoapsis(p, e),
            "Radius is smaller than apoapsis");

    if(!conic_circular(e)) {
        double ff = true_anomaly_from_radius(p, e, r);
        ASSERT(isfinite(ff), "True anomaly not NaN");
        ASSERT_RANGEF(ff, 0, M_PI, "True anomaly range");

        if(zero(f)) // Accuracy is bad near periapsis
            ASSERT(zero(ff*ff), "True anomaly radius identity (zero)");
        else
            ASSERT_EQF(fabs(f), ff, "True anomaly radius identity");
    }

    double v = true_velocity(mu, p, e, f);
    double rdot = true_velocity_radial(mu, p, e, f);
    double rfdot = true_velocity_horizontal(mu, p, e, f);
    ASSERT(isfinite(v) && isfinite(rdot) && isfinite(rfdot),
        "Velocity not NaN");
    ASSERT_EQF(rdot*rdot + rfdot*rfdot, v*v,
        "Velocity magitude");
    ASSERT_EQF(rfdot, r*dfdt, "Horizontal velocity = r * fdot");
    ASSERT(rfdot > 0, "Horizontal velocity is positive");

    double h = conic_specific_angular_momentum(mu, p, e);
    ASSERT_EQF(h, r * rfdot,
        "Specific relative angular momentum");

    double energy = v*v/2.0 - mu/r;
    if(conic_parabolic(e)) { // Accuracy is bad for parabolic trajectory
        ASSERT(zero(energy*energy),
            "Specific orbital energy (vis-viva) is zero (parabolic)");
    } else {
        double visviva = conic_specific_orbital_energy(mu, p, e);
        ASSERT_EQF(visviva, energy,
            "Specific orbital energy (vis-viva)");
    }

    double rplus = true_radius(p, e, f+df);
    double rminus = true_radius(p, e, f-df);
    ASSERT_EQF(rplus - rminus, 2.0*rdot*dt,
        "rdot = dr/dt");

    double tan_phi = true_tan_phi(e, f), phi = true_flight_path_angle(e, f);
    ASSERT(isfinite(tan_phi) && isfinite(phi),
        "Flight path angle not NaN");
    ASSERT_EQF(rdot / rfdot, tan_phi,
        "Flight path angle");
    ASSERT_EQF(tan(phi), tan_phi,
        "Flight path angle tangent");

    ASSERT_EQF(h, r*v*cos(phi),
        "Specific angular momentum and flight path angle");

    double x = true_x(p, e, f);
    double y = true_y(p, e, f);
    double xdot = true_xdot(mu, p, e, f);
    double ydot = true_ydot(mu, p, e, f);
    ASSERT(isfinite(x) && isfinite(y) &&
        isfinite(xdot) && isfinite(ydot),
        "Position and velocity not NaN");

    ASSERT_EQF(x*x + y*y, r*r,
        "Position magnitude");
    ASSERT_EQF(xdot*xdot + ydot*ydot, v*v,
        "Velocity magnitude (xy)");

    ASSERT_EQF(atan2(y, x), f,
        "True anomaly angle");

    ASSERT_EQF(h, x * ydot - y * xdot,
        "Specific relative angular momentum (xy)");

    if(!conic_circular(e))
        ASSERT_EQF(p/e, x + r/e,
            "Focus-Directrix property");

    if(!conic_parabolic(e)) {
        double x2 = e < 1.0 ? x+2.0*c : x-2.0*c;
        double r2 = sqrt(x2*x2 + y*y);
        double sum = e < 1.0 ? r2 + r : r2 - r;

        ASSERT_EQF(e < 1.0 ? 2.0*a : -2.0*a, sum,
            "Focal-radii property");
    }

    if(conic_parabolic(e))
        ASSERT_EQF(2.0*p*(x - q), -(y*y),
            "Trajectory is a parabola");
    else if(conic_hyperbolic(e))
        ASSERT_EQF((x-c)*(x-c)/(a*a), 1.0 + y*y/(b*b),
            "Trajectory is a hyperbola");
    else
        ASSERT_EQF((x+c)*(x+c)/(a*a), 1.0 - y*y/(b*b),
            "Trajectory is an ellipse");

    double xplus = true_x(p, e, f+df);
    double yplus = true_y(p, e, f+df);
    double xminus = true_x(p, e, f-df);
    double yminus = true_y(p, e, f-df);
    double dx = xplus - xminus, dy = yplus - yminus;
    double dx2 = 2.0 * xdot * dt, dy2 = 2.0 * ydot * dt;
    double xx = dx - dx2, yy = dy - dy2;

    ASSERT(ZEROF((xx*xx + yy*yy) / (dx2*dx2 + dy2*dy2)),
        "v = dr/dt");

    double acc = mu / (r*r);
    double dvx2 = -(acc * x / r) * 2.0*dt, dvy2 = -(acc * y / r) * 2.0*dt;
    double xdotplus = true_xdot(mu, p, e, f+df);
    double ydotplus = true_ydot(mu, p, e, f+df);
    double xdotminus = true_xdot(mu, p, e, f-df);
    double ydotminus = true_ydot(mu, p, e, f-df);
    double dvx = xdotplus - xdotminus, dvy = ydotplus - ydotminus;
    double vxx = dvx - dvx2, vyy = dvy - dvy2;

    ASSERT(ZEROF((vxx*vxx + vyy*vyy) / (dvx2*dvx2 + dvy2*dvy2)),
        "a = dv/dt");
}
