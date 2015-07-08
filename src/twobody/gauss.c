#include <math.h>
#include <float.h>

#include <twobody/math_utils.h>
#include <twobody/stumpff.h>
#include <twobody/gauss.h>

#include <assert.h> // XXX: kill me

double gauss_iterate_z(
    double mu,
    double r1, double r2,
    double df,
    double time,
    double z0,
    double *f, double *g,
    double *fdot, double *gdot,
    int max_iterations) {
    double A = sqrt(r1*r2) * sin(df) / sqrt(fmax(0.0, 1.0 - cos(df)));
    double threshold = DBL_EPSILON;

    if(!isfinite(A)) return -1;

    assert(isfinite(A));

    // lower bound value for z, solved from y == 0
    // c1/sqrt(c2) = (r1+r2)/A
    // c1^2/c2 = sinh^2(sqrt(-z)) / (cosh(sqrt(-z)) - 1) = ((r1+r2)/A)^2
    // 2 cosh^2 (sqrt(-z)/2) = ((r1+r2)/A)^2
    double zmin = -square(2.0 * acosh(sign(A)*(r1+r2)/(sqrt(2.0)*A)));
    assert(A >= 0 || (isfinite(zmin) && zmin <= 0.0));

    double z = z0, y = 1.0/0.0;
    for(int iter = 0; iter < max_iterations; ++iter) {
        double cs[4];
        stumpff_fast(z, cs);

        y = r1 + r2 - A * cs[1] / sqrt(cs[2]);
        assert(y >= 0);

        double s = sqrt(fabs(y / cs[2]));
        double t = (s*s*s * cs[3] + A * sqrt(y)) / sqrt(mu);

        assert(isfinite(y));
        assert(isfinite(s));
        assert(isfinite(t));

        // XXX: div by zero if z goes 0 (parabolic / degenerate);
        // TODO: calculate stumpff dc/dz derivatives with series
        //double dSdz = (cs[2] - 3.0*cs[3]) / (2.0 * z);
        //double dCdz = (cs[1] - 2.0*cs[2]) / (2.0 * z);
        double dCdz = stumpff_series_dcdz(2, z);
        double dSdz = stumpff_series_dcdz(3, z);

        assert(isfinite(dSdz));
        assert(isfinite(dCdz));

        double dtdz =
            (s*s*s * (dSdz - 3.0*cs[3]*dCdz / (2.0 * cs[2])) +
            A/8.0 * (3.0 * cs[3] * sqrt(y) / cs[2] + A/s)) /
            sqrt(mu);

        assert(isfinite(dtdz));

        double dz = (time - t) / dtdz;
        if(dz*dz < threshold)
            break;

        if(A < 0.0 || z+dz > zmin)
            z = z + dz;
        else  // don't exceed minimum value (time of flight = 0)
            z = (z + zmin) / 2;
    }

    *f = 1.0 - y/r1;
    *g = A * sqrt(fabs(y) / mu);
    *gdot = 1.0 - y/r2;
    *fdot = (*f * *gdot - 1.0) / *g;

    return z;
}
