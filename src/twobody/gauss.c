#include <math.h>
#include <float.h>

#include <twobody/stumpff.h>
#include <twobody/gauss.h>

#include <assert.h> // XXX: kill me

double gauss_iterate_z(
    double mu,
    double r1, double r2,
    double df,
    double time,
    double z0,
    int max_iterations) {
    double A = sqrt(r1*r2) * sin(df) / sqrt(fmax(0.0, 1.0 - cos(df)));
    double threshold = DBL_EPSILON;

    if(!isfinite(A)) return -1;

    assert(isfinite(A));

    double z = z0;
    for(int iter = 0; iter < max_iterations; ++iter) {
        double cs[4];
        stumpff_fast(z, cs);

        double y = r1 + r2 - A * cs[1] / sqrt(cs[2]);
        double s = sqrt(fabs(y / cs[2]));
        double t = (s*s*s * cs[3] + A * sqrt(fabs(y))) / sqrt(mu);

        assert(isfinite(y));
        assert(isfinite(s));
        assert(isfinite(t));

        // XXX: div by zero if z goes 0 (parabolic / degenerate);
        // TODO: calculate stumpff dc/dz derivatives with series
        double dSdz = (cs[2] - 3.0*cs[3]) / (2.0 * z);
        double dCdz = (cs[1] - 2.0*cs[2]) / (2.0 * z);

        assert(isfinite(dSdz));
        assert(isfinite(dCdz));

        double dtdz = (s*s*s * (dSdz - 3.0*cs[3]*dCdz / (2.0 * cs[2])) +
            A/8.0 * (3.0 * cs[3] * sqrt(fabs(y)) / cs[2] + A/s)) /
            sqrt(mu);

        assert(isfinite(dtdz));

        double dz = (time - t) / dtdz;
        z = z + dz;

        assert(isfinite(z));

        if(dz*dz < threshold)
            break;
    }

    return z;
}
