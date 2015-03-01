#ifndef TWOBODY_UNIVERSAL_H
#define TWOBODY_UNIVERSAL_H

double universal_alpha(double mu, double r, double v2);
double universal_period(double mu, double alpha);

int universal_parabolic(double alpha);
int universal_hyperbolic(double alpha);
int universal_elliptic(double alpha);

double universal_from_eccentric(double p, double e, double dE);

double universal_half_sin_true(
    double p,
    double r0,
    double r,
    double s_half,
    const double *cs_half);
double universal_half_cos_true(
    double p,
    double r0, double sigma0,
    double r,
    double s_half,
    const double *cs_half);
double universal_half_tan_true(
    double p,
    double r0, double sigma0,
    double s_half,
    const double *cs_half);
double universal_to_true(
    double p,
    double r0, double sigma0,
    double r,
    double s_half,
    const double *cs_half);

double universal_time(
    double mu,
    double r0,
    double sigma0,
    double s,
    const double *cs);
double universal_radius(
    double r0,
    double sigma0,
    double s,
    const double *cs);
double universal_sigma(
    double alpha,
    double r0, double sigma0,
    double s,
    const double *cs);

double universal_guess_s(
    double mu,
    double alpha,
    double r0, double sigma0,
    double time);
double universal_iterate_s(
    double mu,
    double alpha,
    double r0, double sigma0,
    double s0, double time,
    int max_steps);

#endif
