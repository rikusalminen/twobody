#ifndef TWOBODY_ORBIT_H
#define TWOBODY_ORBIT_H

#ifndef TWOBODY_NO_SIMD
#include <twobody/simd4d.h>
#endif

struct orbit {
    double gravity_parameter;
    double orbital_energy;
    double angular_momentum;
    double periapsis_time;

#ifndef TWOBODY_NO_SIMD
    vec4d major_axis;
    vec4d minor_axis;
    vec4d normal_axis;
#else
    double major_axis[4];
    double minor_axis[4];
    double normal_axis[4];
#endif
};

void orbit_from_state_ptr(
    struct orbit *orbit,
    double mu,
    const double *pos, const double *vel,
    double epoch);

void orbit_from_elements(
    struct orbit *orbit,
    double mu,
    double p, double e,
    double i, double an, double arg,
    double periapsis_time);

double orbit_gravity_parameter(const struct orbit *orbit);
double orbit_orbital_energy(const struct orbit *orbit);
double orbit_angular_momentum(const struct orbit *orbit);
double orbit_periapsis_time(const struct orbit *orbit);

int orbit_zero(const struct orbit *orbit);
int orbit_radial(const struct orbit *orbit);
int orbit_parabolic(const struct orbit *orbit);
int orbit_hyperbolic(const struct orbit *orbit);
int orbit_elliptic(const struct orbit *orbit);

double orbit_semi_latus_rectum(const struct orbit *orbit);
double orbit_eccentricity(const struct orbit *orbit);

void orbit_state_true(
    const struct orbit *orbit,
    double *pos, double *vel,
    double f);
void orbit_state_eccentric(
    const struct orbit *orbit,
    double *pos, double *vel,
    double E);
void orbit_state_time(
    const struct orbit *orbit,
    double *pos, double *vel,
    double t);

#ifndef TWOBODY_NO_SIMD
#include <twobody/math_utils.h>

#include <twobody/conic.h>
#include <twobody/anomaly.h>
#include <twobody/true_anomaly.h>
#include <twobody/eccentric_anomaly.h>

static inline void orbit_from_state(
    struct orbit *orbit,
    double mu,
    vec4d pos, vec4d vel,
    double epoch)
    __attribute__((always_inline));
static inline void orbit_from_state(
    struct orbit *orbit,
    double mu,
    vec4d pos, vec4d vel,
    double epoch) {
    double r2 = dot(pos, pos);
    double r = sqrt(r2);
    double v2 = dot(vel, vel);

    // specific orbital energy
    double visviva = v2/2.0 - mu/r;

    // specific relative angular momentum
    vec4d h = cross(pos, vel);

    if(dot(h, h) < DBL_EPSILON) {       // radial trajectory
        vec4d major = r2 < DBL_EPSILON ?
            (vec4d){ 1.0, 0.0, 0.0, 0.0 } : unit4d(pos);
        vec4d up = (major[0]*major[0] + major[1]*major[1] < DBL_EPSILON) ?
            (vec4d){ 0.0, 1.0, 0.0, 0.0 } :
            (vec4d){ 0.0, 0.0, 1.0, 0.0 };
        vec4d minor = cross(up, major);
        vec4d normal = cross(major, minor);

        orbit->gravity_parameter = mu;
        orbit->orbital_energy = visviva;
        orbit->angular_momentum = 0.0;
        orbit->periapsis_time = NAN; // TODO: radial trajectory equation
        orbit->major_axis = major;
        orbit->minor_axis = minor;
        orbit->normal_axis = normal;
    } else {                            // conic trajectory
        // semi-latus rectum
        double p = dot(h, h) / mu;

        // eccentricity vector
        vec4d ecc = splat4d(1.0/mu) *
            (splat4d(v2 - mu/r)*pos - dot4d(pos, vel) * vel);
        double e2 = dot(ecc, ecc);
        double e = sqrt(e2);
        int circular = e2 < DBL_EPSILON;

        // orbit normal vector
        vec4d normal = unit4d(h);

        // ascending node
        vec4d nodes = (vec4d){ -normal[1], normal[0], 0.0, 0.0 };
        int equatorial = dot(nodes, nodes) < DBL_EPSILON;

        vec4d major = circular && equatorial ?
            (vec4d){ 1.0, 0.0, 0.0, 0.0 } : // circular & equatorial -> x-axis
            (circular ?  unit4d(nodes) :    // circular -> line of nodes
                unit4d(ecc));               // eccentricity vector
        vec4d minor = cross(normal, major);

        // true anomaly
        double f0 = -sign(dot(vel, major)) *
            acos(clamp(-1.0, 1.0, dot(major, pos) / r));

        // mean anomaly and mean motion
        double M0 = anomaly_true_to_mean(e, f0);
        double n = conic_mean_motion(mu, p, e);

        // periapsis time
        double t0 = epoch - M0 / n;

        orbit->gravity_parameter = mu;
        orbit->orbital_energy = visviva;
        orbit->angular_momentum = mag(h);
        orbit->periapsis_time = t0;
        orbit->major_axis = major;
        orbit->minor_axis = minor;
        orbit->normal_axis = normal;
    }
}

static inline vec4d orbit_position_true(const struct orbit *orbit, double f)
    __attribute__((always_inline));
static inline vec4d orbit_position_true(const struct orbit *orbit, double f) {
    double p = orbit_semi_latus_rectum(orbit);
    double e = orbit_eccentricity(orbit);

    return splat4d(true_x(p, e, f)) * orbit->major_axis +
        splat4d(true_y(p, e, f)) * orbit->minor_axis;
}

static inline vec4d orbit_velocity_true(const struct orbit *orbit, double f)
    __attribute__((always_inline));
static inline vec4d orbit_velocity_true(const struct orbit *orbit, double f) {
    double mu = orbit_gravity_parameter(orbit);
    double p = orbit_semi_latus_rectum(orbit);
    double e = orbit_eccentricity(orbit);

    return splat4d(true_xdot(mu, p, e, f)) * orbit->major_axis +
        splat4d(true_ydot(mu, p, e, f)) * orbit->minor_axis;
}

static inline vec4d orbit_position_eccentric(const struct orbit *orbit, double E)
    __attribute__((always_inline));
static inline vec4d orbit_position_eccentric(const struct orbit *orbit, double E) {
    double p = orbit_semi_latus_rectum(orbit);
    double e = orbit_eccentricity(orbit);

    return splat4d(eccentric_x(p, e, E)) * orbit->major_axis +
        splat4d(eccentric_y(p, e, E)) * orbit->minor_axis;
}

static inline vec4d orbit_velocity_eccentric(const struct orbit *orbit, double E)
    __attribute__((always_inline));
static inline vec4d orbit_velocity_eccentric(const struct orbit *orbit, double E) {
    double mu = orbit_gravity_parameter(orbit);
    double p = orbit_semi_latus_rectum(orbit);
    double e = orbit_eccentricity(orbit);

    return splat4d(eccentric_xdot(mu, p, e, E)) * orbit->major_axis +
        splat4d(eccentric_ydot(mu, p, e, E)) * orbit->minor_axis;
}

#endif

#endif
