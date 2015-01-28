#ifndef TWOBODY_STUMPFF_H
#define TWOBODY_STUMPFF_H

double stumpff_c0(double alpha, double s);
double stumpff_c1(double alpha, double s);
double stumpff_c2(double alpha, double s);
double stumpff_c3(double alpha, double s);
double stumpff_series(int k, double z);

void stumpff_fast(double z, double *cs);

#endif
