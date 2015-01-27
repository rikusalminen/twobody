#include <twobody/conic.h>
#include <twobody/anomaly.h>
#include <twobody/math_utils.h>

#include <math.h>
#include <float.h>

static inline double anomaly_iter1(double e, double M, double x) {
    return M + e * sin(x);
}

static inline double anomaly_iter2(double e, double M, double x) {
    return x + (M + e * sin(x) - x) / (1.0 - e * cos(x));
}

static inline double anomaly_iter3(double e, double M, double x) {
    double s = e * sin(x);
    double c = e * cos(x);
    double f0 = x - s - M;
    double f1 = 1.0 - c;
    double f2 = s;

    return x + (-5.0) * f0 /
        (f1 + sign(f1) * sqrt(fabs(16.0 * f1 * f1 - 20.0 * f0 * f2)));
}

static inline double anomaly_iter4(double e, double M, double x) {
    double s = e * sinh(x);
    double c = e * cosh(x);
    double f0 = s - x - M;
    double f1 = c - 1.0;
    double f2 = s;

    return x + (-5.0) * f0 /
        (f1 + sign(f1) * sqrt(fabs(16.0 * f1 * f1 - 20.0 * f0 * f2)));
}

double anomaly_eccentric_iterate(double e, double M, double E0, int max_steps) {
    typedef double (*iter_func)(double, double, double);
    iter_func iter = 0;
    if(e < 0.3) iter = anomaly_iter1;
    else if(e < 0.9) iter = anomaly_iter2;
    else if(e < 1.0) iter = anomaly_iter3;
    else iter = anomaly_iter4;   // e > 1.0

    int num_steps = 10;
    if(max_steps > 0) num_steps = max_steps;
    else if(e > 0.0 && e < 0.3) num_steps = 10;
    else if(e < 0.9) num_steps = 20;
    else if(e < 1.0) num_steps = 20;
    else num_steps = 31;    // e > 1.0

    double threshold = DBL_EPSILON;

    double x = E0, x0 = x;
    do {
        x0 = x;
        x = iter(e, M, x);
    } while(num_steps >= 0 && (x0-x)*(x0-x) > threshold);

    return x;
}

double anomaly_mean_to_eccentric(double e, double M) {
    if(conic_parabolic(e)) {
        // parabolic anomaly
        double x = pow(sqrt(9.0*M*M + 1.0) + 3.0*M, 1.0/3.0);
        return x - 1.0/x;
    }

    // eccentric or hyperbolic anomaly
    return anomaly_eccentric_iterate(e, M, M, 0);
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
