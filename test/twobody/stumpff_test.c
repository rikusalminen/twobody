#include <twobody/stumpff.h>

#include <math.h>

#include "../numtest.h"

void stumpff_test(
    double *params,
    int num_params,
    void *extra_args,
    struct numtest_ctx *test_ctx) {
    (void)extra_args;
    ASSERT(num_params == 2, "");

    const double maxs = 4.0*M_PI, maxalpha = 1.0;
    double s = (-1.0 + 2.0 * params[0]) * maxs;
    double alpha = (-1.0 + 2.0 * params[1]) * maxalpha;

    double z = alpha * s*s;

    // test Stumpff series
    double cs[4];
    for(int i = 0; i < 4; ++i)
        cs[i] = stumpff_series(i, z);

    for(int i = 0; i < 4; ++i)
        ASSERT(isfinite(cs[i]),
            "Stumpff series c%d not NaN", i);

    for(int i = 0; i < 2; ++i)
        ASSERT_EQF(1.0 - cs[i], z*cs[i+2],
            "Stumpff series c%d recurrence relation", i);

    // derivatives of Stumpff series
    for(int i = 0; i < 4; ++i) {
        double ds = 1.0e-4;

        double dcds = (i == 0) ?
            -alpha * s * cs[1] : // d/ds c_0(z) = -alpha * s * c_1(z)
            pow(s, i-1) * cs[i-1]; // d/ds s^i c_i(z) = s^(i-1) c_i-1(z)

        double splus = s + ds, sminus = s - ds;
        double zplus = alpha * splus * splus, zminus = alpha * sminus * sminus;
        double cplus = pow(splus, i) * stumpff_series(i, zplus),
           cminus = pow(sminus, i) * stumpff_series(i, zminus);
        double dc = cplus - cminus;

        ASSERT_EQF(2.0 * dcds * ds, dc,
            "d/ds c%d(z)", i);
    }

    // test Stumpff functions (trigonometric)
    double cs_fun[4]  = {
        stumpff_c0(z),
        stumpff_c1(z),
        stumpff_c2(z),
        stumpff_c3(z)
    };

    for(int i = 0; i < 4; ++i)
        ASSERT(isfinite(cs_fun[i]),
            "Stumpff function c%d not NaN", i);

    for(int i = 0; i < 2; ++i)
        ASSERT_EQF(1.0 - cs_fun[i], z*cs_fun[i+2],
            "Stumpff function c%d recurrence relation", i);

    for(int i = 0; i < 4; ++i)
        ASSERT_EQF(cs[i], cs_fun[i],
            "Stumpff function and series c%d are equal", i);

    // test Stumpff derivatives (trigonometric)
    double cs_dz[4] = {
        stumpff_dc0dz(z),
        stumpff_dc1dz(z),
        stumpff_dc2dz(z),
        stumpff_dc3dz(z)
    };

    for(int i = 0; i < 4; ++i)
        ASSERT(isfinite(cs_dz[i]),
            "Stumpff derivative function c%d not NaN", i);

    for(int i = 0; i < 4; ++i) {
        double dz = 1.0e-5;
        double zplus = z + dz, zminus = z - dz;

        double cplus = stumpff_series(i, zplus);
        double cminus = stumpff_series(i, zminus);
        double dc = cplus - cminus;
        double dcdz = cs_dz[i];
        double dcf = 2.0 * dcdz * dz;

        if(fabs(z) >= 1.0e-4) // accuracy is bad for small z
            ASSERT_EQF(dcf, dc,
                "d/dz c%d(z) function", i);
        else
            ASSERT(ZEROF((dcf-dc)*(dcf-dc)),
                "d/dz c%d(z) function (near zero)", i);
    }

    // test Stumpff derivatives (series)
    double cs_dz_s[4];
    for(int i = 0; i < 4; ++i)
        cs_dz_s[i] = stumpff_series_dcdz(i, z);

    for(int i = 0; i < 4; ++i)
        ASSERT(isfinite(cs_dz_s[i]),
            "Stumpff derivative series c%d not NaN", i);

    for(int i = 0; i < 4; ++i) {
        if(fabs(z) >= 1.0e-4) // accuracy is bad for small z
            ASSERT_EQF(cs_dz_s[i], cs_dz[i],
                "Stumpff derivative function and series c%d are equal", i);
        else
            ASSERT(fabs(cs_dz_s[i]-cs_dz[i])/cs_dz_s[i] < 1.0e-5,
                "Stumpff derivative function and series c%d are equal"
                " (near zero)", i);
    }

    for(int i = 0; i < 4; ++i) {
        double dz = 1.0e-4;
        double zplus = z + dz, zminus = z - dz;

        double cplus = stumpff_series(i, zplus);
        double cminus = stumpff_series(i, zminus);
        double dc = cplus - cminus;
        double dcdz = cs_dz_s[i];

        ASSERT_EQF(2.0 * dcdz * dz, dc,
            "d/dz c%d(z) series", i);
    }

    // test fast Stumpff series
    double cs_fast[4];
    stumpff_fast(z, cs_fast);

    for(int i = 0; i < 4; ++i)
        ASSERT(isfinite(cs_fast[i]),
            "Stumpff fast series c%d not NaN", i);

    for(int i = 0; i < 2; ++i)
        ASSERT_EQF(1.0 - cs_fast[i], z*cs_fast[i+2],
            "Stumpff fast series c%d recurrence relation", i);

    for(int i = 0; i < 4; ++i)
        ASSERT_EQF(cs[i], cs_fast[i],
            "Stumpff fast series and series c%d are equal", i);

    double cs_four[4];
    stumpff_fast(4.0 * z, cs_four);

    // quadruple angle formulae
    ASSERT_EQF(cs_four[0], 2.0*cs[0]*cs[0] - 1.0,
        "c_0(4z) = 2*c_0(z)^2 - 1");
    ASSERT_EQF(cs_four[1], cs[0]*cs[1],
        "c_1(4z) = c_0(z)*c_1(z)");
    ASSERT_EQF(cs_four[2], cs[1]*cs[1] / 2.0,
        "c_2(4z) = c_1(z)^2 / 2");
    ASSERT_EQF(cs_four[3], (cs[2] + cs[0]*cs[3])/4.0,
        "c_3(4z) = (c_2(z) + c_0(z)*c_3(z))/4");
}
