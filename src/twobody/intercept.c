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
