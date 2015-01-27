#include <twobody/orbit.h>
#include <twobody/conic.h>
#include <twobody/anomaly.h>
#include <twobody/true_anomaly.h>
#include <twobody/orientation.h>

#include "../numtest.h"

void orbit_from_state_test(
    double *params,
    int num_params,
    void *extra_args,
    struct numtest_ctx *test_ctx) {
    (void)extra_args;
    ASSERT(num_params == 7, "");

    double mu = 1.0 + params[0] * 1.0e10;
    double p = 1.0 + params[1] * 1.0e10;
    double e = params[2] * 4.0;

    double maxf = anomaly_eccentric_to_true(e, M_PI);
    double f = (-1.0 + 2.0 * params[3]) * maxf;

    double rx = (-1.0 + 2.0 * params[4]) * M_PI;
    double ry = (-1.0 + 2.0 * params[5]) * M_PI;
    double rz = (-1.0 + 2.0 * params[6]) * M_PI;

    double t = 0.0;

    vec4d major_axis = {
        cos(ry)*cos(rz),
        cos(ry)*sin(rz),
        -sin(ry),
        0.0 };
    vec4d minor_axis = {
        sin(rx)*sin(ry)*cos(rz) - cos(rx)*sin(rz),
        sin(rx)*sin(ry)*sin(rz) + cos(rx)*cos(rz),
        sin(rx)*cos(ry),
        0.0 };

    vec4d radial =
        splat4d(cos(f)) * major_axis +
        splat4d(sin(f)) * minor_axis;
    vec4d horizontal =
        -splat4d(sin(f)) * major_axis +
        splat4d(cos(f)) * minor_axis;

    if(conic_circular(e)) { // "fix" circular orbits
        vec4d normal = cross(major_axis, minor_axis);
        vec4d nodes = { -normal[1], normal[0], 0.0, 0.0 };

        if(zero(dot(nodes, nodes))) { // circular and equatorial
            major_axis = (vec4d){ 1.0, 0.0, 0.0, 0.0 };
            minor_axis = (vec4d){ 0.0, sign(normal[2]), 0.0, 0.0 };
        } else {
            major_axis = unit4d(nodes);
            minor_axis = cross(normal, major_axis);
        }

        f = -sign(dot(major_axis, horizontal)) *
            acos(clamp(-1.0, 1.0, dot(major_axis, radial)));
    }

    double r = true_radius(p, e, f);
    double vr = true_velocity_radial(mu, p, e, f);
    double vh = true_velocity_horizontal(mu, p, e, f);

    vec4d pos = splat4d(r) * radial;
    vec4d vel = splat4d(vr) * radial + splat4d(vh) * horizontal;

    struct orbit orbit;
    orbit_from_state(&orbit, mu, pos, vel, t);

    ASSERT(!orbit_zero(&orbit) && !orbit_radial(&orbit),
        "Orbit is not degenerate");

    ASSERT(ZEROF(dot(orbit.major_axis, orbit.minor_axis)) &&
        ZEROF(dot(orbit.major_axis, orbit.normal_axis)) &&
        ZEROF(dot(orbit.minor_axis, orbit.normal_axis)),
        "Axes are orthogonal");

    ASSERT(eqv4d(
            cross(orbit.major_axis, orbit.minor_axis),
            orbit.normal_axis),
        "Axis cross product");

    ASSERT(eqv4d(orbit.major_axis, major_axis) &&
        eqv4d(orbit.minor_axis, minor_axis) &&
        eqv4d(orbit.normal_axis, cross(major_axis, minor_axis)),
        "Orbit orientation");

    if(conic_parabolic(e)) { // accuracy is bad for parabola
        ASSERT(ZEROF(square(orbit_orbital_energy(&orbit))),
            "Specific orbital energy is zero (parabola)");
    } else {
        ASSERT_EQF(
            conic_specific_orbital_energy(mu, p, e),
            orbit_orbital_energy(&orbit),
            "Specific orbital energy");
    }

    ASSERT_EQF(
        conic_specific_angular_momentum(mu, p, e),
        orbit_angular_momentum(&orbit),
        "Specific relative angular momentum energy");

    ASSERT_EQF(p, orbit_semi_latus_rectum(&orbit),
        "Semi latus rectum");
    ASSERT_EQF(e, orbit_eccentricity(&orbit),
        "Semi latus rectum");

    if(!conic_circular(e)) {
        double v2 = dot(vel, vel);
        vec4d ecc = splat4d(1.0/mu) *
            (splat4d(v2 - mu/r)*pos - dot4d(pos, vel) * vel);

        ASSERT_EQF(dot(ecc, orbit.major_axis), e,
            "Eccentricity vector");
    }

    double dt = t - orbit_periapsis_time(&orbit);
    double M = dt * conic_mean_motion(mu, p, e);
    double E = anomaly_mean_to_eccentric(e, M);
    double f2 = anomaly_eccentric_to_true(e, E);

    if(EQF(fabs(f), M_PI)) {
        ASSERT_EQF(fabs(f), fabs(f2),
            "True anomaly (apoapsis)");
    } else {
        ASSERT_EQF(f, f2,
            "True anomaly");
    }

    ASSERT(eqv4d(pos, orbit_position_true(&orbit, f)),
        "Orbit position (true anomaly)");
    ASSERT(eqv4d(vel, orbit_velocity_true(&orbit, f)),
        "Orbit velocity (true anomaly)");

    ASSERT(eqv4d(pos, orbit_position_eccentric(&orbit, E)),
        "Orbit position (eccentric anomaly)");
    ASSERT(eqv4d(vel, orbit_velocity_eccentric(&orbit, E)),
        "Orbit velocity (eccentric anomaly)");
}

void orbit_from_elements_test(
    double *params,
    int num_params,
    void *extra_args,
    struct numtest_ctx *test_ctx) {
    (void)extra_args;
    ASSERT(num_params == 6, "");

    double mu = 1.0 + params[0] * 1.0e5;
    double p = 1.0 + params[1] * 1.0e5;
    double e = params[2] * 4.0;
    double i = params[3] * M_PI;
    double an = (ZEROF(i) || ZEROF(i - M_PI))  ? 0.0 :
        (-1.0 + 2.0*params[4]) * M_PI;
    double arg = (-1.0 + 2.0*params[5]) * M_PI;
    double t0 = 0.0;

    struct orbit orbit;
    orbit_from_elements(&orbit, mu, p, e, i, an, arg, t0);

    ASSERT(!orbit_radial(&orbit) && !orbit_zero(&orbit),
        "Orbit is not degenerate");

    ASSERT_EQF(orbit_gravity_parameter(&orbit), mu,
        "Orbit gravity parameter");
    ASSERT_EQF(orbit_semi_latus_rectum(&orbit), p,
        "Orbit semi-latus rectum");
    ASSERT_EQF(orbit_eccentricity(&orbit), e,
        "Orbit eccentricity");
    ASSERT_EQF(orbit_periapsis_time(&orbit), t0,
        "Orbit eccentricity");

    ASSERT(orbit_parabolic(&orbit) == conic_parabolic(e),
        "Orbit is parabolic");
    ASSERT(orbit_hyperbolic(&orbit) == conic_hyperbolic(e),
        "Orbit is hyperbolic");
    ASSERT(orbit_elliptic(&orbit) == conic_elliptic(e),
        "Orbit is elliptic");

    ASSERT_EQF(i,
        orientation_inclination(
            orbit.major_axis,
            orbit.minor_axis,
            orbit.normal_axis),
        "Orbit inclination");
    ASSERT_EQF(an,
        orientation_longitude_of_ascending_node(
            orbit.major_axis,
            orbit.minor_axis,
            orbit.normal_axis),
        "Orbit longitude of ascending node");
    ASSERT_EQF(arg,
        orientation_argument_of_periapsis(
            orbit.major_axis,
            orbit.minor_axis,
            orbit.normal_axis),
        "Orbit argument of perigee");

    ASSERT_EQF(orbit_orbital_energy(&orbit),
        conic_specific_orbital_energy(mu, p, e),
        "Orbit specific orbital energy");
    ASSERT_EQF(orbit_angular_momentum(&orbit),
        conic_specific_angular_momentum(mu, p, e),
        "Orbit specific angular momentum");
}

void orbit_radial_test(
    double *params,
    int num_params,
    void *extra_args,
    struct numtest_ctx *test_ctx) {
    (void)extra_args;
    ASSERT(num_params == 5, "");

    double mu = 1.0 + params[0] * 1.0e10;
    double r = params[1] * 1.0e5;
    double coeff = params[2];
    double t0 = 0.0;

    double ry = (-1.0 + 2.0 * params[3]) * M_PI;
    double rz = (-1.0 + 2.0 * params[4]) * M_PI;

    vec4d major_axis = {
        cos(ry)*cos(rz),
        cos(ry)*sin(rz),
        -sin(ry),
        0.0 };

    double v = zero(r) ?
        coeff * sqrt(2.0 * mu) :
        (0.5 + coeff) * sqrt(2.0 * mu/r);

    vec4d pos = splat4d(r) * major_axis;
    vec4d vel = splat4d(v) * major_axis;

    struct orbit orbit;
    orbit_from_state(&orbit, mu, pos, vel, t0);

    ASSERT(orbit_radial(&orbit),
        "Orbit is radial");
    ASSERT(orbit_zero(&orbit) == zero(r),
        "Orbit is zero iff position is zero");

    ASSERT(ZEROF(dot(orbit.major_axis, orbit.minor_axis)) &&
        ZEROF(dot(orbit.major_axis, orbit.normal_axis)) &&
        ZEROF(dot(orbit.minor_axis, orbit.normal_axis)),
        "Axes are orthogonal");

    ASSERT(ZEROF(orbit_angular_momentum(&orbit)),
        "Radial trajectory orbit has zero angular momentum");

    if(orbit_zero(&orbit))
        return;

    ASSERT_EQF(dot(orbit.major_axis, pos), r,
        "Radial trajectory major axis");

    ASSERT(eqv4d(major_axis, orbit.major_axis),
        "Radial trajectory major axis");

    double visviva = v*v/2.0 - mu/r;
    ASSERT_EQF(visviva, orbit_orbital_energy(&orbit),
        "Radial trajectory energy");

    // TODO: implement and test radial trajectories
}
