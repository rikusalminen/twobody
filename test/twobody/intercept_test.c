#include <twobody/intercept.h>
#include <twobody/conic.h>
#include <twobody/anomaly.h>
#include <twobody/eccentric_anomaly.h>
#include <twobody/orbit.h>

#include "../numtest.h"

void intercept_test(
    double *params,
    int num_params,
    void *extra_args,
    struct numtest_ctx *test_ctx) {
    (void)extra_args;
    ASSERT(num_params == 10, "");

    double mu = 1.0 + params[0] * 1.0e5;
    double p1 = 1.0 + params[1] * 1.0e5;
    double e1 = params[2] * 2.0;
    double i = params[3] * M_PI;
    double an = (ZEROF(i) || ZEROF(i - M_PI))  ? 0.0 :
        (-1.0 + 2.0*params[4]) * M_PI;
    double arg = (-1.0 + 2.0*params[5]) * M_PI;
    double E1 = (-1.0 + 2.0*params[6]) *
        conic_closed(e1) ? M_PI : M_PI*0.5;
    double f1 = anomaly_eccentric_to_true(e1, E1);

    struct orbit orbit1;
    orbit_from_elements(&orbit1, mu, p1, e1, i, an, arg, 0.0);

    double M1 = anomaly_eccentric_to_mean(e1, E1);
    double n1 = conic_mean_motion(mu, p1, e1);
    double t = M1/n1 + orbit_periapsis_time(&orbit1);

    double r1 = eccentric_radius(p1, e1, E1);
    vec4d pos1 = orbit_position_eccentric(&orbit1, E1);
    vec4d normal = orbit1.normal_axis;
    vec4d radial = unit4d(pos1), horizontal = cross(normal, radial);

    double e2 = params[7] * 2.0;
    double E2 = (-1.0 + 2.0*params[8]) *
        conic_closed(e2) ? M_PI : M_PI*0.5;
    double f2 = anomaly_eccentric_to_true(e2, E2);
    double reli = (-1.0 + 2.0*params[8]) * M_PI;
    double p2 = r1 * (1.0 + e2 * cos(f2));

    double vr = eccentric_velocity_radial(mu, p2, e2, E2);
    double vh = eccentric_velocity_horizontal(mu, p2, e2, E2);
    vec4d vel2 = splat4d(vr) * radial +
        splat4d(cos(reli) * vh) * horizontal +
        splat4d(sin(reli) * vh) * normal;

    struct orbit orbit2;
    orbit_from_state(&orbit2, mu, pos1, vel2, t);

    if(zero(e2)) { // "fix" true anomaly for circular orbits
        f2 = -sign(dot(vel2, orbit2.major_axis)) *
            acos(clamp(-1.0, 1.0, dot(orbit2.major_axis, radial)));
        E2 = anomaly_true_to_eccentric(e2, f2);
    }

    double n2 = conic_mean_motion(mu, p2, e2);
    double M2 = (t-orbit_periapsis_time(&orbit2)) * n2;
    double EE2 = anomaly_mean_to_eccentric(e2, M2);

    ASSERT(EQF(E2, EE2) || (EQF(fabs(E2), M_PI) && EQF(fabs(EE2), M_PI)),
        "Eccentric anomaly sanity");

    vec4d pos2 = orbit_position_eccentric(&orbit2, E2);
    ASSERT(eqv4d(pos1, pos2), "Intercept position is sane");

    double threshold = (p1 + p2) / 1000.0;
    double intersect_fs[8];
    const struct orbit *orbits[2] = { &orbit1, &orbit2 };
    for(int i = 0; i < 2; ++i) {
        double *fs = intersect_fs + 4*i;
        int is = intercept_intersect(orbits[i], orbits[!i], threshold, fs);

        ASSERT(is == 1 || is == 2,
            "1 or 2 intersects");

        if(is <= 0 || is > 2)
            continue;

        for(int j = 0; j < is; ++j) {
            double f0 = fs[j*2], f1 = fs[j*2+1];

            ASSERT_RANGEF(f0, j == 0 ? -2.0*M_PI : -M_PI, M_PI,
                "True anomaly range begin");
            ASSERT_RANGEF(f1, -M_PI, M_PI,
                "True anomaly range end");

            ASSERT_LTF(f0, f1, "Range begin is less than range end");
        }

        for(int j = is; j < 2; ++j) {
            double f0 = fs[j*2], f1 = fs[j*2+1];
            ASSERT_LTF(f1, f0, "Range is empty");
        }

        ASSERT(fs[0] >= -M_PI || conic_closed(i == 0 ? e1 : e2),
            "Intersect may overlap apoapsis on closed orbits only");

        double f = i == 0 ? f1 : f2;
        ASSERT(((LTF(fs[0], f) && LTF(f, fs[1])) ||
            ((fs[0] <= -M_PI || ZEROF(fs[0]+M_PI)) &&
                LTF(fs[0], f-2.0*M_PI) && LTF(f-2.0*M_PI, fs[1]))) ||
            (is == 2 && LTF(fs[2], f) && LTF(f, fs[3])),
            "True anomaly is within range");

        if(is == 2) {
            ASSERT_LTF(fs[1], fs[2], "First range is before second");
        }
    }

    double Mmax1 = conic_closed(e1) ? 2.0*M_PI :
        anomaly_eccentric_to_mean(e1, M_PI);
    double Mmax2 = conic_closed(e2) ? 2.0*M_PI :
        anomaly_eccentric_to_mean(e2, M_PI);
    double t0 = fmax(
        orbit_periapsis_time(&orbit1) - Mmax1/n1,
        orbit_periapsis_time(&orbit2) - Mmax2/n2);
    double t1 = fmin(
        orbit_periapsis_time(&orbit1) + Mmax1/n1,
        orbit_periapsis_time(&orbit2) + Mmax2/n2);

    int max_times = 8;
    double times[max_times * 2];

    int ts = intercept_times(
        &orbit1, &orbit2,
        t0, t1,
        intersect_fs,
        times, max_times);

    ASSERT(ts >= 1, "1 or more intervals of intercept");
    ASSERT(ts <= max_times, "At most max_times intervals of intercept");

    if(ts < 0 || ts > max_times) return;

    int time_interval = -1;
    for(int i = 0; i < ts; ++i) {
        ASSERT_LTF(times[2*i+0], times[2*i+1],
            "Time range not empty");

        ASSERT_RANGEF(times[2*i+0], t0, t1,
            "Time range in search interval (begin)");
        ASSERT_RANGEF(times[2*i+1], t0, t1,
            "Time range in search interval (end)");

        if(LTF(times[2*i+0], t) && LTF(t, times[2*i+1]))
            time_interval = i;
    }

    ASSERT(time_interval >= 0 && time_interval < ts,
        "Intercept time interval found");

    if(time_interval == -1) return;

    t0 = times[time_interval*2+0]; t1 = times[time_interval*2+1];

    int search_steps = 30;
    double target_distance = 0.0;
    struct intercept intercept;
    double t_end = intercept_search(
        &orbit1, &orbit2,
        t0, t1,
        threshold, target_distance,
        search_steps,
        &intercept);
    (void)t_end; // TODO: loop over range and search for correct intercept?

    ASSERT_LTF(intercept.distance, threshold,
        "Intercept distance is less than threshold");
}
