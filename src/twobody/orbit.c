#include <twobody/orbit.h>
#include <twobody/conic.h>
#include <twobody/orientation.h>

#include <twobody/math_utils.h>

void orbit_from_elements(
    struct orbit *orbit,
    double mu,
    double p, double e,
    double i, double an, double arg,
    double periapsis_time) {
    orbit->gravity_parameter = mu;
    orbit->orbital_energy = conic_specific_orbital_energy(mu, p, e);
    orbit->angular_momentum = conic_specific_angular_momentum(mu, p, e);
    orbit->periapsis_time = periapsis_time;
    orbit->major_axis = orientation_major_axis(i, an, arg);
    orbit->minor_axis = orientation_minor_axis(i, an, arg);
    orbit->normal_axis = orientation_normal_axis(i, an, arg);
}

double orbit_gravity_parameter(const struct orbit *orbit) {
    return orbit->gravity_parameter;
}

double orbit_orbital_energy(const struct orbit *orbit) {
    return orbit->orbital_energy;
}

double orbit_angular_momentum(const struct orbit *orbit) {
    return orbit->angular_momentum;
}

double orbit_periapsis_time(const struct orbit *orbit) {
    return orbit->periapsis_time;
}

int orbit_zero(const struct orbit *orbit) {
    return !isfinite(orbit->orbital_energy) &&
        zero(orbit->angular_momentum);
}

int orbit_radial(const struct orbit *orbit) {
    return zero(orbit->angular_momentum);
}

int orbit_parabolic(const struct orbit *orbit) {
    return zero(orbit->orbital_energy);
}

int orbit_hyperbolic(const struct orbit *orbit) {
    return !orbit_parabolic(orbit) &&
        orbit->orbital_energy > 0;
}

int orbit_elliptic(const struct orbit *orbit) {
    return !orbit_parabolic(orbit) &&
        orbit->orbital_energy < 0;
}

double orbit_semi_latus_rectum(const struct orbit *orbit) {
    double mu = orbit->gravity_parameter;
    double h = orbit->angular_momentum;
    return h*h / mu;
}

double orbit_eccentricity(const struct orbit *orbit) {
    double mu = orbit->gravity_parameter;
    double h = orbit->angular_momentum;
    double ee = orbit->orbital_energy;
    return sqrt(fmax(0.0, 1.0 + 2.0*ee*h*h / (mu*mu)));
}

void orbit_state_true(
    const struct orbit *orbit,
    double *pos, double *vel,
    double f) {
    *(vec4d*)pos = orbit_position_true(orbit, f);
    *(vec4d*)vel = orbit_velocity_true(orbit, f);
}

void orbit_state_eccentric(
    const struct orbit *orbit,
    double *pos, double *vel,
    double E) {
    *(vec4d*)pos = orbit_position_eccentric(orbit, E);
    *(vec4d*)vel = orbit_velocity_eccentric(orbit, E);
}

void orbit_state_time(
    const struct orbit *orbit,
    double *pos, double *vel,
    double t) {
    double mu = orbit_gravity_parameter(orbit);
    double p = orbit_semi_latus_rectum(orbit);
    double e = orbit_eccentricity(orbit);

    double dt = t - orbit->periapsis_time;
    double M = dt * conic_mean_motion(mu, p, e);
    double E = anomaly_mean_to_eccentric(e, M);

    orbit_state_eccentric(orbit, pos, vel, E);
}
