#ifndef TWOBODY_ECCENTRIC_ANOMALY_H
#define TWOBODY_ECCENTRIC_ANOMALY_H

double eccentric_radius(double p, double e, double E);
double eccentric_anomaly_from_radius(double p, double e, double r);

double eccentric_dEdt(double mu, double p, double e, double E);

double eccentric_time(double mu, double p, double e, double E);

double eccentric_velocity(double mu, double p, double e, double E);
double eccentric_velocity_radial(double mu, double p, double e, double E);
double eccentric_velocity_horizontal(double mu, double p, double e, double E);

double eccentric_tan_phi(double e, double E);
double eccentric_flight_path_angle(double e, double E);

double eccentric_x(double p, double e, double E);
double eccentric_y(double p, double e, double E);
double eccentric_xdot(double mu, double p, double e, double E);
double eccentric_ydot(double mu, double p, double e, double E);

double eccentric_f(
    double mu, double p, double e,
    double r0,
    double dE);
double eccentric_g(
    double mu, double p, double e,
    double r0, double sigma0,
    double dE);
double eccentric_g_t(
    double mu, double p, double e,
    double dE, double dt);
double eccentric_fdot(
    double mu, double p, double e,
    double r0, double r,
    double dE);
double eccentric_gdot(
    double mu, double p, double e,
    double r,
    double dE);

#endif
