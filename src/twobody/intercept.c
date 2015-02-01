#include <twobody/conic.h>
#include <twobody/anomaly.h>
#include <twobody/true_anomaly.h>
#include <twobody/orbit.h>
#include <twobody/intercept.h>

static int intersect_ranges(
    const double *fs1, const double *fs2,
    int closed,
    double *fs) {
    // calculate two angle ranges fs[0]..fs[1] and fs[2]..fs[3]
    // where ranges fs1 (two pairs) and fs2 (two pairs) overlap (-2pi..2pi)
    // fs[0] = -2pi..pi, fs[1..3] = -pi..pi
    // (only first range may overlap apoapsis)
    for(int i = 0; i < 2; ++i) {
        for(int j = 0; j < 2; ++j) {
            double f0 = fs1[(fs1[2] < fs1[3] ? 2*i : 0) + j];
            double f1 = fs2[(fs2[2] < fs2[3] ? 2*i : 0) + j];
            fs[i*2+j] = j ? fmin(f0, f1) : fmax(f0, f1);
        }
    }

    if(closed &&
        (fs[0] <= -M_PI || zero(fs[0]+M_PI)) &&
        (fs[3] >= M_PI || zero(fs[3]-M_PI))) {
        // ranges overlap at apoapsis -> union ranges -2pi
        fs[0] = fmin(fs[0], fs[2]-2.0*M_PI);
        fs[1] = fmax(fs[1], fs[3]-2.0*M_PI);
        fs[2] = 1.0; fs[3] = -1.0;
    }

    if((fs[1] >= fs[2] || zero(fs[1]-fs[2])) &&
        fs[2] < fs[3]) {
        // ranges overlap at periapsis -> union ranges
        fs[1] = fs[3];
        fs[2] = 1.0; fs[3] = -1.0;
    }

    if(fs[2] < fs[3] && fs[3] > M_PI) {
        // one or two ranges, second overlaps apoapsis -> swap ranges -2pi
        double f0 = fs[0], f1 = fs[1];
        fs[0] = fs[2]-2.0*M_PI;
        fs[1] = fs[3]-2.0*M_PI;
        fs[2] = f0; fs[3] = f1;
    }

    if(fs[2] < fs[3] && !(fs[0] < fs[1])) {
        // first range is empty, second is not -> swap ranges
        fs[0] = fs[2]; fs[1] = fs[3];
        fs[2] = 1.0; fs[3] = -1.0;
    }

    if(fs[1] - fs[0] >= 2.0*M_PI) {
        // first range covers full orbit -> -pi..pi
        // (second range must be empty)
        fs[0] = -M_PI; fs[1] = M_PI;
    }

    return (fs[0] < fs[1]) + (fs[2] < fs[3]);
}

int intercept_intersect(
    const struct orbit *orbit1,
    const struct orbit *orbit2,
    double threshold,
    double *fs) {
    // find 0, 1 or 2 ranges of true anomaly (fs[0]..fs[1]) and (fs[2]..fs[3])
    // where orbit1 is between the periapsis and apoapsis of orbit2
    // and closer than threshold to the orbital plane of orbit2

    // radial orbits not handled
    if(orbit_radial(orbit1) || orbit_radial(orbit2))
        return 0;

    // find 0, 1, or 2 ranges of true anomaly (w.r.t. orbit1)
    // where orbit1 and orbit2 may be than threshold

    double p1 = orbit_semi_latus_rectum(orbit1);
    double p2 = orbit_semi_latus_rectum(orbit2);
    double e1 = orbit_eccentricity(orbit1);
    double e2 = orbit_eccentricity(orbit2);

    // apoapsis-periapsis test
    double ap1 = conic_apoapsis(p1, e1), pe1 = conic_periapsis(p1, e1);
    double ap2 = conic_apoapsis(p2, e2), pe2 = conic_periapsis(p2, e2);
    if( (conic_closed(e1) && ap1 <= pe2 - threshold) ||
        (conic_closed(e2) && ap2 <= pe1 - threshold))
        return 0;

    // true anomaly between target apoapsis/periapsis +/- threshold
    double maxf = M_PI; //conic_max_true_anomaly(e1); // XXX: maxf
    double fpe = conic_circular(e1) ? 0.0 :
        true_anomaly_from_radius(p1, e1, pe2 - threshold);
    double fap = fmin(maxf, (conic_circular(e1) || !conic_closed(e2)) ? M_PI :
        true_anomaly_from_radius(p1, e1, ap2 + threshold));

    double f1 = fmin(fap, fpe), f2 = fmax(fap, fpe);

    double fs1[4] = { 2.0*M_PI, -2.0*M_PI, 2.0*M_PI, -2.0*M_PI };
    if(conic_closed(e1) && zero(f1) && !(f2 < M_PI)) {
        // intersects anywhere on orbit (f = -pi .. pi)
        fs1[0] = -2.0*M_PI; fs1[1] = 2.0*M_PI;
    } else if(zero(f1)) {
        // intersect near periapsis (f = -f2 .. f2)
        fs1[0] = -f2; fs1[1] = f2;
    } else if(conic_closed(e1) && !(f2 < M_PI)) {
        // intersect near apoapsis (f < -f1, f > f1)
        fs1[0] = -2.0*M_PI; fs1[1] = -f1;
        fs1[2] = f1; fs1[3] = 2.0*M_PI;
    } else {
        // two intersects (-f2 < f < -f1, f1 < f < f2)
        fs1[0] = -f2; fs1[1] = -f1;
        fs1[2] = f1; fs1[3] = f2;
    }

    // ascending node
    vec4d nodes = cross(orbit1->normal_axis, orbit2->normal_axis);
    double N2 = dot(nodes, nodes), N = sqrt(N2);
    int coplanar = N2 < DBL_EPSILON;

    double fs2[4] = { -M_PI, M_PI, 1.0, -1.0 };
    if(!coplanar) {
        // relative inclination
        double reli = sign(dot(orbit1->normal_axis, orbit2->normal_axis)) *
            asin(clamp(-1.0, 1.0, N));

        // true anomaly of ascending/descending node
        double f_an = sign(dot(orbit1->minor_axis, nodes)) *
            acos(clamp(-1.0, 1.0, dot(orbit1->major_axis, nodes)/N));
        double f_dn = f_an - sign(f_an) * M_PI;

        for(int i = 0; i < 2; ++i) {
            double f_node = i == 0 ? fmin(f_an, f_dn) : fmax(f_an, f_dn);
#if 1
            // distance at node
            double r = p1 / (1.0 + e1*cos(f_node));
#else
            // periapsis = conservative estimate
            double r = conic_periapsis(p1, e1);
#endif
            // spherical trigonometry sine law
            double delta_f = asin(
                clamp(-1.0, 1.0,
                    sin(threshold / (2.0*r)) / sin(fabs(reli) / 2.0)));

            fs2[i*2+0] = f_node - delta_f;
            fs2[i*2+1] = f_node + delta_f;
        }
    }

    return intersect_ranges(fs1, fs2, conic_closed(e1), fs);
}

#include <assert.h> // XXX: kill me

int intercept_times(
    const struct orbit *orbit1,
    const struct orbit *orbit2,
    double t0, double t1,
    const double *fs,
    double *intercept_times,
    int max_times) {
    // find at most max_times time ranges between (t0..t1) where
    // orbit1 is between (fs[0]..fs[1]) or (fs[2]..fs[3]) and
    // orbit2 is between (fs[4]..fs[5]) or (fs[6]..fs[7])

    double mu = orbit_gravity_parameter(orbit1);
    const struct orbit *orbits[2] = { orbit1, orbit2 };

    // find time ranges corresponding to true anomaly ranges
    double times[2][4] = { { 1.0, -1.0, 1.0, -1.0 }, { 1.0, -1.0, 1.0, -1.0 } };
    double periods[2] = { 0.0, 0.0 }; // orbital period
    int n_orbit[2] = { 0, 0 }; // number of complete orbits to t0
    for(int o = 0; o < 2; ++o) {
        double p = orbit_semi_latus_rectum(orbits[o]);
        double e = orbit_eccentricity(orbits[o]);
        double t_pe = orbit_periapsis_time(orbits[o]);
        double n = conic_mean_motion(mu, p, e);

        for(int i = 0; i < 2; ++i) {
            double f0 = fs[4*o+2*i+0], f1 = fs[4*o+2*i+1];
            if(f0 >= f1)
                continue;

            for(int j = 0; j < 2; ++j) {
                double f = fs[4*o+2*i+j];
                double E = anomaly_true_to_eccentric(e, f);
                double M = anomaly_eccentric_to_mean(e, E);

                if(f < -M_PI) {
                    assert(conic_closed(e));
                    assert(E >= -M_PI && E <= M_PI); // XXX: kill
                    assert(M >= -M_PI && M <= M_PI); // XXX: kill
                    M -= 2.0*M_PI;
                }

                times[o][2*i+j] = t_pe + M/n;
            }
        }

        if(conic_closed(e)) {
            double P = conic_period(mu, p, e);
            periods[o] = P;
            n_orbit[o] = (int)trunc((t0 - t_pe)/P + (t_pe > t0 ? -0.5 : 0.5));
        }
    }

    int isect[2] = { 0, 0 };
    int num_times = 0;
    double t = t0;
    while(t < t1 && num_times < max_times) {
        double trange[2][2];

        // time interval on this orbital period
        for(int o = 0; o < 2; ++o) {
            double period = n_orbit[o] * periods[o];
            trange[o][0] = times[o][2*isect[o]+0] + period;
            trange[o][1] = times[o][2*isect[o]+1] + period;
        }

        // overlapping time interval
        double t_begin = fmax(t, fmax(trange[0][0], trange[1][0]));
        double t_end = fmin(t1, fmin(trange[0][1], trange[1][1]));
        t = t_end;

        // non-empty interval found
        if(t_begin < t_end) {
            intercept_times[2*num_times+0] = t_begin;
            intercept_times[2*num_times+1] = t_end;
            num_times += 1;
        }

        // advance to next intersect range
        int advance = trange[0][1] < trange[1][1] ? 0 : 1;
        isect[advance] += 1;

        if(isect[advance] == 2 ||
            times[advance][2*isect[advance]+0] >= // XXX: use true anomaly range instead?
            times[advance][2*isect[advance]+1]) {
            // advance to next orbit
            if(!conic_closed(orbit_eccentricity(orbits[advance]))) // XXX: orbit_closed!
                break; // open orbit, search exhausted

            isect[advance] = 0;
            n_orbit[advance] += 1;
        }
    }

    return num_times;
}


#if 0

int intercept_search(
    const struct orbit *orbit1,
    const struct orbit *orbit2,
    double t0, double t1,
    double threshold,
    double *times) {
    // find 0, 1, or 2 time ranges
    // where distance between orbit1 and orbit2 is less than threshold
}

struct intercept {
    double time;
    vec4d relative_position;
    vec4d relative_velocity;
};

double intercept_minimize(
    const struct orbit *orbit1,
    const struct orbit *orbit2,
    double t0, double t1,
    double target_distance,
    struct intercept *intercept) {
    // find 1 time
    // when distance between orbit1 and orbit2 is closest to target_distance
    // there must be an intercept between t0 and t1
}

int intercept_orbit(
    const struct orbit *orbit1,
    const struct orbit *orbit2,
    double t0, double t1,
    double threshold,
    double target_distance,
    struct intercept *intercepts,
    int max_intercepts) {
    // find at most max_intercepts times
    // when distance between orbit1 and orbit2 is closest to target_distance

    // geometry prefilter
    // XXX: this must be done twice!
    double fs[4];
    int intersects = intercept_intersect(orbit1, orbit2, threshold, fs);
    if(intersects == 0)
        return 0;

    int num_intercepts = 0;
    while(t0 < t1 && num_intercepts < max_intercepts) {
        // time prefilter
        int max_times = max_intercepts;
        double times[2*max_times];
        int num_times = intercept_times(
            orbit1, orbit2,
            t0, t1,
            fs,
            times, max_times);

        if(num_times == 0)
            break;

        for(int i = 0; i < num_times; ++i) {
            double ti[4];
            int is = intercept_search(
                orbit1, orbit2,
                times[2*i+0], times[2*i+1],
                threshold, ti);

            for(int j = 0; j < is; ++i) {
                intercept_minimize(
                    orbit1, orbit2,
                    ti[2*j+0], ti[2*j+1],
                    target_distance,
                    intercepts + num_intercepts);
                t0 = ti[2*j+1];
            }

            num_intercepts += is;
        }
    }

    return num_intercepts;
}

#endif
