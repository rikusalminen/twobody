#include <twobody/conic.h>
#include <twobody/eccentric_anomaly.h>
#include <twobody/math_utils.h>

#include <math.h>
#include <float.h>

#include "../numtest.h"

void eccentric_anomaly_test(
    double *params,
    int num_params,
    void *extra_args,
    struct numtest_ctx *test_ctx) {
    (void)extra_args;

    ASSERT(num_params == 4, "");

    double mu = 1.0 + params[0] * 1.0e10;
    double p = 1.0 + params[1] * 1.0e10;
    double e = params[2] * 4.0;

    double E = (-1.0 + params[3] * 2.0) * M_PI;

    double a = conic_semi_major_axis(p, e);
    double b = conic_semi_minor_axis(p, e);
    double c = conic_focal_distance(p, e);
    double q = conic_periapsis(p, e);

    double n = conic_mean_motion(mu, p, e);
    double dt = (2.0 * M_PI / n) * (1.0 / 36000.0); // 1/100 degree

    double r = eccentric_radius(p, e, E);
    ASSERT(isfinite(r), "Radius not NaN");
    ASSERT(r > 0, "Radius is positive");

    ASSERT_LTF(q, r, "Radius is larger than periapsis");
    if(conic_closed(e))
        ASSERT_LTF(r, p / (1.0 - e), "Radius is smaller than apoapsis");

    if(!conic_circular(e)) {
        double EE = eccentric_anomaly_from_radius(p, e, r);
        ASSERT(isfinite(EE), "Eccentric anomaly not NaN");
        ASSERT_RANGEF(EE, 0, M_PI, "Eccentric anomaly range");

        if(ZEROF(E)) // Accuracy is bad near periapsis
            ASSERT(ZEROF(EE*EE), "Eccentric anomaly radius identity (zero)");
        else
            ASSERT_EQF(fabs(E), EE, "Eccentric anomaly radius identity");
    }

    double dEdt = eccentric_dEdt(mu, p, e, E);
    double dE = dEdt * dt;
    ASSERT(isfinite(dEdt), "dE/dt not NaN");
    ASSERT(dEdt > 0.0, "dE/dt is positive");

    double t = eccentric_time(mu, p, e, E);
    ASSERT(isfinite(t), "Eccentric anomaly time of flight not NaN");

    double tplus = eccentric_time(mu, p, e, E+dE);
    double tminus = eccentric_time(mu, p, e, E-dE);
    ASSERT_EQF(tplus - tminus, 2.0*dt, "dE/dt");

    double v = eccentric_velocity(mu, p, e, E);
    double rdot = eccentric_velocity_radial(mu, p, e, E);
    double rfdot = eccentric_velocity_horizontal(mu, p, e, E);
    ASSERT(isfinite(v) && isfinite(rdot) && isfinite(rfdot),
        "Velocity not NaN");
    ASSERT_EQF(rdot*rdot + rfdot*rfdot, v*v,
        "Velocity magitude");
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

    double rplus = eccentric_radius(p, e, E+dE);
    double rminus = eccentric_radius(p, e, E-dE);
    ASSERT_EQF(rplus - rminus, 2.0*rdot*dt,
        "rdot = dr/dt");

    double tan_phi = eccentric_tan_phi(e, E);
    double phi = eccentric_flight_path_angle(e, E);

    ASSERT(isfinite(tan_phi) && isfinite(phi),
        "Flight path angle not NaN");
    ASSERT_EQF(rdot / rfdot, tan_phi,
        "Flight path angle");
    ASSERT_EQF(tan(phi), tan_phi,
        "Flight path angle tangent");

    ASSERT_EQF(h, r*v*cos(phi),
        "Specific angular momentum and flight path angle");

    double x = eccentric_x(p, e, E);
    double y = eccentric_y(p, e, E);
    double xdot = eccentric_xdot(mu, p, e, E);
    double ydot = eccentric_ydot(mu, p, e, E);
    ASSERT(isfinite(x) && isfinite(y) &&
        isfinite(xdot) && isfinite(ydot),
        "Position and velocity not NaN");

    ASSERT_EQF(x*x + y*y, r*r,
        "Position magnitude");
    ASSERT_EQF(xdot*xdot + ydot*ydot, v*v,
        "Velocity magnitude (xy)");

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

    double dotrv = x*xdot + y*ydot;
    if(conic_parabolic(e)) {
        ASSERT_EQF(dotrv / sqrt(mu*p), E,
            "Parabolic anomaly r dot v");
        ASSERT_EQF(r/q - 1.0, E*E,
            "Parabolic anomaly squared");
    } else if(conic_hyperbolic(e)) {
        ASSERT_EQF(dotrv / sqrt(-mu*a), e*sinh(E),
            "Hyperbolic anomaly sine");
        ASSERT_EQF(1.0 - r/a, e*cosh(E),
            "Hyperbolic anomaly cosine");
    } else {
        ASSERT_EQF(dotrv / sqrt(mu*a), e*sin(E),
            "Eccentric anomaly sine");
        ASSERT_EQF(1.0 - r/a, e*cos(E),
            "Eccentric anomaly cosine");
    }

    double xplus = eccentric_x(p, e, E+dE);
    double yplus = eccentric_y(p, e, E+dE);
    double xminus = eccentric_x(p, e, E-dE);
    double yminus = eccentric_y(p, e, E-dE);
    double dx = xplus - xminus, dy = yplus - yminus;
    double dx2 = 2.0 * xdot * dt, dy2 = 2.0 * ydot * dt;
    double xx = dx - dx2, yy = dy - dy2;

    ASSERT(ZEROF((xx*xx + yy*yy) / (dx2*dx2 + dy2*dy2)),
        "v = dr/dt");

    double acc = mu / (r*r);
    double dvx2 = -(acc * x / r) * 2.0*dt, dvy2 = -(acc * y / r) * 2.0*dt;
    double xdotplus = eccentric_xdot(mu, p, e, E+dE);
    double ydotplus = eccentric_ydot(mu, p, e, E+dE);
    double xdotminus = eccentric_xdot(mu, p, e, E-dE);
    double ydotminus = eccentric_ydot(mu, p, e, E-dE);
    double dvx = xdotplus - xdotminus, dvy = ydotplus - ydotminus;
    double vxx = dvx - dvx2, vyy = dvy - dvy2;

    ASSERT(ZEROF((vxx*vxx + vyy*vyy) / (dvx2*dvx2 + dvy2*dvy2)),
        "a = dv/dt");
}
