#ifndef TWOBODY_GAUSS_H
#define TWOBODY_GAUSS_H

double gauss_iterate_z(
    double mu,
    double r1, double r2,
    double df,
    double dt,
    double z0,
    int max_iterations);

#endif
