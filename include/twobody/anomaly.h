#ifndef TWOBODY_ANOMALY_H
#define TWOBODY_ANOMALY_H

double anomaly_eccentric_iterate(double e, double M, double E0, int max_steps);
double anomaly_mean_to_eccentric(double e, double M);
double anomaly_eccentric_to_mean(double e, double E);
double anomaly_eccentric_to_true(double e, double E);
double anomaly_true_to_eccentric(double e, double f);
double anomaly_true_to_mean(double e, double f);
double anomaly_mean_to_true(double e, double M);
double anomaly_dEdM(double e, double E);
double anomaly_dfdE(double e, double E);

double anomaly_true_sin(double e, double E);
double anomaly_true_cos(double e, double E);
double anomaly_true_tan_half(double e, double E);
double anomaly_eccentric_sin(double e, double f);
double anomaly_eccentric_cos(double e, double f);
double anomaly_eccentric_tan_half(double e, double f);

#endif
