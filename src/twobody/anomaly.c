#include <twobody/conic.h>
#include <twobody/anomaly.h>
#include <twobody/math_utils.h>

#include <math.h>
#include <float.h>

double anomaly_eccentric_iterate(double e, double M, double E0, int max_steps) {
    if(max_steps <= 0)
        max_steps = e < 1.0 ? 10 : 20;
    double threshold = DBL_EPSILON;

    double Mperiod = 0.0;
    if(conic_elliptic(e)) {
        // mean anomaly -pi..pi (faster covergence), Mperiod is multiple of 2pi
        double MM = angle_clamp(M);
        Mperiod = M - MM;
        M = MM;
        E0 = E0 - Mperiod;
    }

    double E = E0;
    for(int step = 0; step < max_steps; ++step) {
        double f0, f1, f2;

        if(e < 1.0) { // elliptic
            f0 = E - e*sin(E) - M;
            f1 = 1.0 - e*cos(E);
            f2 = e*sin(E);
        } else { // hyperbolic
            f0 = e*sinh(E) - E - M;
            f1 = e*cosh(E) - 1.0;
            f2 = e*sinh(E);
        }

        double N = 5.0; // laguerre-conway magic constant
        double dE = -N * f0 /
            (f1 + sign(f1) * sqrt(fabs(square(N-1.0) * f1*f1 - N*(N-1.0) * f0*f2)));
        E = E + dE;

        if(dE*dE < threshold)
            break;
    }

    return E + Mperiod;
}

double anomaly_mean_to_eccentric(double e, double M) {
    if(conic_parabolic(e)) {
        // parabolic anomaly
        double x = pow(sqrt(9.0*M*M + 1.0) + 3.0*M, 1.0/3.0);
        return x - 1.0/x;
    }

    double E0 = M;
    if(e > 1.0) // hyperbolic orbit
        E0 = sign(M) * log(2.0 * fabs(M) / e + 1.85);
    else if(e > 0.9) // high eccentricity
        E0 = M + 0.85 * e * sign(angle_clamp(M));

    // eccentric or hyperbolic anomaly
    return anomaly_eccentric_iterate(e, M, E0, 0);
}

double anomaly_eccentric_to_mean(double e, double E) {
    if(conic_parabolic(e))
        return E*E*E/6.0 + E/2.0;
    else if(conic_hyperbolic(e))
        return e * sinh(E) - E;
    else
        return E - e * sin(E);
}

double anomaly_eccentric_to_true(double e, double E) {
    if(conic_parabolic(e))
        return 2.0 * atan(E);
    else if(conic_hyperbolic(e))
        return 2.0 * atan(sqrt((e+1.0) / (e-1.0)) * tanh(E/2.0));
    else
        return atan2(sqrt(1.0-e*e) * sin(E), cos(E) - e);
}

double anomaly_true_to_eccentric(double e, double f) {
    if(conic_parabolic(e))
        return tan(f / 2.0);
    else if(conic_hyperbolic(e))
        return 2.0 * atanh(sqrt((e-1.0) / (e+1.0)) * tan(f/2.0));
    else
        return atan2(sqrt(1.0-e*e) * sin(f), cos(f) + e);
}

double anomaly_true_to_mean(double e, double f) {
    return anomaly_eccentric_to_mean(e, anomaly_true_to_eccentric(e, f));
}

double anomaly_mean_to_true(double e, double M) {
    return anomaly_eccentric_to_true(e, anomaly_mean_to_eccentric(e, M));
}

double anomaly_dEdM(double e, double E) {
    if(conic_parabolic(e))
        return 2.0 / (E*E + 1.0);
    else if(conic_hyperbolic(e))
        return 1.0 / (e*cosh(E) - 1.0);
    else
        return 1.0 / (1.0 - e * cos(E));
}

double anomaly_dfdE(double e, double E) {
    if(conic_parabolic(e))
        return 2.0 / (E*E + 1.0);
    else if(conic_hyperbolic(e))
        return sqrt(e*e - 1.0) / (e*cosh(E) - 1.0);
    else
        return sqrt(1.0 - e*e) / (1.0 - e * cos(E));
}

double anomaly_true_sin(double e, double E) {
    if(conic_parabolic(e))
        return 2.0*E / (E*E + 1.0);
    else if(conic_hyperbolic(e))
        return sqrt(e*e - 1.0) * sinh(E) / (e*cosh(E) - 1.0);
    else
        return sqrt(1.0 - e*e) * sin(E) / (1.0 - e*cos(E));
}

double anomaly_true_cos(double e, double E) {
    if(conic_parabolic(e))
        return (1.0 - E*E) / (1.0 + E*E);
    else if(conic_hyperbolic(e))
        return (e - cosh(E)) / (e*cosh(E) - 1.0);
    else
        return (cos(E) - e) / (1.0 - e*cos(E));
}

double anomaly_true_tan_half(double e, double E) {
    if(conic_parabolic(e))
        return E;
    else if(conic_hyperbolic(e))
        return sqrt((e + 1.0)/(e - 1.0)) * tanh(E/2.0);
    else
        return sqrt((1.0 + e)/(1.0 - e)) * tan(E/2.0);
}

double anomaly_eccentric_sin(double e, double f) {
    if(conic_parabolic(e))
        return 1.0 / 0.0; // TODO: parabolic
    else if(conic_hyperbolic(e))
        return sqrt(e*e - 1.0) * sin(f) / (1.0 + e*cos(f));
    else
        return sqrt(1.0 - e*e) * sin(f) / (1.0 + e*cos(f));
}

double anomaly_eccentric_cos(double e, double f) {
    if(conic_parabolic(e))
        return 1.0 / 0.0; // TODO: parabolic
    else
        // hyperbolic or elliptic (equal)
        return (e + cos(f)) / (1.0 + e*cos(f));
}

double anomaly_eccentric_tan_half(double e, double f) {
    if(conic_parabolic(e))
        return 1.0 / 0.0; // TODO: parabolic
    else if(conic_hyperbolic(e))
        return sqrt((e - 1.0)/(e + 1.0)) * tan(f/2.0);
    else
        return sqrt((1.0 - e)/(1.0 + e)) * tan(f/2.0);
}
