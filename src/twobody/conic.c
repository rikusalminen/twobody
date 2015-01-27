#include <twobody/conic.h>
#include <twobody/math_utils.h>

#include <math.h>
#include <float.h>

int conic_circular(double e) {
    return zero(e);
}

int conic_parabolic(double e) {
    return zero(e - 1.0);
}

int conic_elliptic(double e) {
    return !conic_parabolic(e) && e < 1.0;
}

int conic_hyperbolic(double e) {
    return !conic_parabolic(e) && e > 1.0;
}

int conic_closed(double e) {
    return !conic_parabolic(e) && e < 1.0;
}

double conic_semi_major_axis(double p, double e) {
    if(conic_parabolic(e))
        return INFINITY;
    return p / (1.0 - e*e);
}

double conic_semi_minor_axis(double p, double e) {
    if(conic_parabolic(e))
        return INFINITY;
    else if(conic_hyperbolic(e))
        return p / sqrt(e*e - 1.0);
    else
        return p / sqrt(1.0 - e*e);
}

double conic_focal_distance(double p, double e) {
    double a = conic_semi_major_axis(p, e);

    if(conic_parabolic(e))
        return INFINITY;
    else if(conic_hyperbolic(e))
        return -a * e;
    else
        return a * e;
}

double conic_periapsis(double p, double e) {
    return p / (1.0 + e);
}

double conic_apoapsis(double p, double e) {
    if(!conic_closed(e))
        return INFINITY;
    return p / (1.0 - e);
}

double conic_periapsis_velocity(double mu, double p, double e) {
    return sqrt(mu / p) * (1.0 + e);
}

double conic_apoapsis_velocity(double mu, double p, double e) {
    if(!conic_closed(e))
        return NAN;
    return sqrt(mu / p) * (1.0 - e);
}

double conic_max_true_anomaly(double e) {
    if(conic_hyperbolic(e))
        return M_PI - acos(fmin(1.0, 1.0/e));
    return M_PI;
}

double conic_mean_motion(double mu, double p, double e) {
    double a = conic_semi_major_axis(p, e);

    if(conic_parabolic(e))
        return sqrt(mu / (p*p*p));
    else if(conic_hyperbolic(e))
        return sqrt(mu / -(a*a*a));
    else
        return sqrt(mu / (a*a*a));
}

double conic_period(double mu, double p, double e) {
    if(!conic_closed(e))
        return INFINITY;

    return 2.0 * M_PI / conic_mean_motion(mu, p, e);
}


double conic_specific_orbital_energy(double mu, double p, double e) {
    if(conic_parabolic(e))
        return 0.0;

    double a = conic_semi_major_axis(p, e);
    return -mu / (2.0 * a);
}

double conic_specific_angular_momentum(double mu, double p, double e) {
    (void)e;
    return sqrt(mu * p);
}
