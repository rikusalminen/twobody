#include <twobody/conic.h>
#include <twobody/eccentric_anomaly.h>
#include <twobody/math_utils.h>

#include <math.h>

double eccentric_radius(double p, double e, double E) {
    double a = conic_semi_major_axis(p, e);
    double q = conic_periapsis(p, e);

    if(conic_parabolic(e))
        return q * (E*E + 1.0);
    else if(conic_hyperbolic(e))
        return a * (1.0 - e*cosh(E));
    else
        return a * (1.0 - e*cos(E));
}

double eccentric_anomaly_from_radius(double p, double e, double r) {
    double a = conic_semi_major_axis(p, e);
    double q = conic_periapsis(p, e);

    if(conic_parabolic(e))
        return sqrt(fmax(0.0, r/q - 1.0));
    else if(conic_hyperbolic(e))
        return acosh(fmax(1.0, (1.0 - r/a) / e));
    else
        return acos(clamp(-1.0, 1.0, (1.0 - r/a) / e));
}

double eccentric_dEdt(double mu, double p, double e, double E) {
    double a = conic_semi_major_axis(p, e);

    if(conic_parabolic(e))
        return sqrt(mu / p) / eccentric_radius(p, e, E);
    else if(conic_hyperbolic(e))
        return sqrt(mu / -a) / eccentric_radius(p, e, E);
    else
        return sqrt(mu / a) / eccentric_radius(p, e, E);
}

double eccentric_time(double mu, double p, double e, double E) {
    double a = conic_semi_major_axis(p, e);

    if(conic_parabolic(e))
        return sqrt(p*p*p / mu) * (E*E*E / 6.0 + E / 2.0);
    else if(conic_hyperbolic(e))
        return sqrt(-a*a*a / mu) * (e*sinh(E) - E);
    else
        return sqrt(a*a*a / mu) * (E - e*sin(E));
}

double eccentric_velocity(double mu, double p, double e, double E) {
    double a = conic_semi_major_axis(p, e);

    if(conic_parabolic(e))
        return sqrt((mu/p) * 4.0 / (E*E + 1.0));
    else if(conic_hyperbolic(e))
        return sqrt((mu/-a) * (e*cosh(E) + 1.0) / (e*cosh(E) - 1.0));
    else
        return sqrt((mu/a) * (1.0 + e*cos(E)) / (1.0 - e*cos(E)));
}

double eccentric_velocity_radial(double mu, double p, double e, double E) {
    double a = conic_semi_major_axis(p, e);

    if(conic_parabolic(e))
        return sqrt(mu/p) * 2.0*E / (E*E + 1.0);
    else if(conic_hyperbolic(e))
        return sqrt(mu/-a) * e*sinh(E) / (e*cosh(E) - 1.0);
    else
        return sqrt(mu/a) * e*sin(E) / (1.0 - e*cos(E));
}

double eccentric_velocity_horizontal(double mu, double p, double e, double E) {
    double a = conic_semi_major_axis(p, e);

    if(conic_parabolic(e))
        return sqrt(mu/p) * 2.0 / (E*E + 1.0);
    else if(conic_hyperbolic(e))
        return sqrt((mu/-a) * (e*e - 1.0)) / (e*cosh(E) - 1.0);
    else
        return sqrt((mu/a) * (1.0 - e*e)) / (1.0 - e*cos(E));
}

double eccentric_tan_phi(double e, double E) {
    if(conic_parabolic(e))
        return E;
    else if(conic_hyperbolic(e))
        return e*sinh(E) / sqrt(e*e - 1.0);
    else
        return e*sin(E) / sqrt(1.0 - e*e);
}

double eccentric_flight_path_angle(double e, double E) {
    return atan(eccentric_tan_phi(e, E));
}

double eccentric_x(double p, double e, double E) {
    double a = conic_semi_major_axis(p, e);
    double q = conic_periapsis(p, e);

    if(conic_parabolic(e))
        return q * (1.0 - E*E);
    else if(conic_hyperbolic(e))
        return a * (cosh(E) - e);
    else
        return a * (cos(E) - e);
}

double eccentric_y(double p, double e, double E) {
    double b = conic_semi_minor_axis(p, e);

    if(conic_parabolic(e))
        return p * E;
    else if(conic_hyperbolic(e))
        return b * sinh(E);
    else
        return b * sin(E);
}

double eccentric_xdot(double mu, double p, double e, double E) {
    double a = conic_semi_major_axis(p, e);

    if(conic_parabolic(e))
        return sqrt(mu / p) * -2.0*E / (E*E + 1.0);
    else if(conic_hyperbolic(e))
        return sqrt(mu / -(a*a*a)) * a*sinh(E) / (e*cosh(E) - 1.0);
    else
        return sqrt(mu / (a*a*a)) * -a*sin(E) / (1.0 - e*cos(E));
}

double eccentric_ydot(double mu, double p, double e, double E) {
    double a = conic_semi_major_axis(p, e);
    double b = conic_semi_minor_axis(p, e);

    if(conic_parabolic(e))
        return sqrt(mu / p) * 2.0 / (E*E + 1.0);
    else if(conic_hyperbolic(e))
        return sqrt(mu / -(a*a*a)) * b*cosh(E) / (e*cosh(E) - 1.0);
    else
        return sqrt(mu / (a*a*a)) * b*cos(E) / (1.0 - e*cos(E));
}

double eccentric_f(
    double mu, double p, double e,
    double r0,
    double dE) {
    (void)mu;
    double a = conic_semi_major_axis(p, e);
    double q = conic_periapsis(p, e);

    if(conic_parabolic(e))
        return 1.0 - q/r0 * dE*dE;
    else if(conic_hyperbolic(e))
        return 1.0 - a/r0 * (1.0 - cosh(dE));
    else
        return 1.0 - a/r0 * (1.0 - cos(dE));
}

double eccentric_g(
    double mu, double p, double e,
    double r0, double sigma0,
    double dE) {
    double a = conic_semi_major_axis(p, e);
    double q = conic_periapsis(p, e);

    if(conic_parabolic(e))
        return r0 * sqrt(p/mu) * dE +
            sigma0/sqrt(mu) * q * dE*dE;
    else if(conic_hyperbolic(e))
        return r0 * sqrt(-a/mu) * sinh(dE) +
            sigma0/sqrt(mu) * a * (1.0 - cosh(dE));
    else
        return r0 * sqrt(a/mu) * sin(dE) +
            sigma0/sqrt(mu) * a * (1.0 - cos(dE));
}

double eccentric_g_t(
    double mu, double p, double e,
    double dE, double dt) {
    double a = conic_semi_major_axis(p, e);

    if(conic_parabolic(e))
        return dt - sqrt(p*p*p/mu) * (1.0/6.0) * dE*dE*dE;
    else if(conic_hyperbolic(e))
        return dt - sqrt(-a*a*a/mu) * (sinh(dE) - dE);
    else
        return dt - sqrt(a*a*a/mu) * (dE - sin(dE));
}

double eccentric_fdot(
    double mu, double p, double e,
    double r0, double r,
    double dE) {
    double a = conic_semi_major_axis(p, e);

    if(conic_parabolic(e))
        return -sqrt(mu*p)/(r*r0) * dE;
    else if(conic_hyperbolic(e))
        return -sqrt(mu*-a)/(r*r0) * sinh(dE);
    else
        return -sqrt(mu*a)/(r*r0) * sin(dE);
}

double eccentric_gdot(
    double mu, double p, double e,
    double r,
    double dE) {
    (void)mu;
    double a = conic_semi_major_axis(p, e);
    double q = conic_periapsis(p, e);

    if(conic_parabolic(e))
        return 1.0 - q/r * dE*dE;
    else if(conic_hyperbolic(e))
        return 1.0 - a/r * (1.0 - cosh(dE));
    else
        return 1.0 - a/r * (1.0 - cos(dE));
}
