#include <math.h>
#include <twobody/stumpff.h>
#include <twobody/gauss.h>

double gauss_iterate_z(
    double mu,
    double r1, double r2,
    double df,
    double dt,
    double z0,
    int max_iterations) {
    double A = sqrt(r1*r2) * sin(df) / sqrt(1.0 - cos(df));

    double z = z0;
    for(int iter = 0; iter < max_iterations; ++iter) {
        double cs[4];
        stumpff_fast(z, cs);

        double y = r1 + r2 - A * (1.0 - z * cs[3]) / sqrt(cs[2]);
        double s = sqrt(y / cs[2]);

        double t = (s*s*s * cs[3] + A * sqrt(y)) / sqrt(mu);
    }

    return z;
}
