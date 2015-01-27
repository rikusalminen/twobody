#ifndef TWOBODY_CONIC_H
#define TWOBODY_CONIC_H

int conic_circular(double e);
int conic_elliptic(double e);
int conic_parabolic(double e);
int conic_hyperbolic(double e);
int conic_closed(double e);

double conic_semi_major_axis(double p, double e);
double conic_semi_minor_axis(double p, double e);
double conic_focal_distance(double p, double e);
double conic_periapsis(double p, double e);
double conic_apoapsis(double p, double e);
double conic_periapsis_velocity(double mu, double p, double e);
double conic_apoapsis_velocity(double mu, double p, double e);

double conic_max_true_anomaly(double e);

double conic_mean_motion(double mu, double p, double e);
double conic_period(double mu, double p, double e);

double conic_specific_orbital_energy(double mu, double p, double e);
double conic_specific_angular_momentum(double mu, double p, double e);

#endif
