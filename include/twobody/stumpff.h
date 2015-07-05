#ifndef TWOBODY_STUMPFF_H
#define TWOBODY_STUMPFF_H

double stumpff_c0(double z);
double stumpff_c1(double z);
double stumpff_c2(double z);
double stumpff_c3(double z);

double stumpff_dc0dz(double z);
double stumpff_dc1dz(double z);
double stumpff_dc2dz(double z);
double stumpff_dc3dz(double z);

double stumpff_series(int k, double z);
double stumpff_series_dcdz(int k, double z);

void stumpff_fast(double z, double *cs);

#endif
