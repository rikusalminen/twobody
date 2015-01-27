#include <twobody/orientation.h>

void orientation_major_axis_ptr(double *axis, double i, double an, double arg) {
    *(vec4d*)axis = orientation_major_axis(i, an, arg);
}

void orientation_minor_axis_ptr(double *axis, double i, double an, double arg) {
    *(vec4d*)axis = orientation_minor_axis(i, an, arg);
}

void orientation_normal_axis_ptr(double *axis, double i, double an, double arg) {
    *(vec4d*)axis = orientation_normal_axis(i, an, arg);
}


double orientation_inclination_ptr(
    const double *major,
    const double *minor,
    const double *normal) {
    return orientation_inclination(
        *(vec4d*)major, *(vec4d*)minor, *(vec4d*)normal);
}

double orientation_longitude_of_ascending_node_ptr(
    const double *major,
    const double *minor,
    const double *normal) {
    return orientation_longitude_of_ascending_node(
        *(vec4d*)major, *(vec4d*)minor, *(vec4d*)normal);
}

double orientation_argument_of_periapsis_ptr(
    const double *major,
    const double *minor,
    const double *normal) {
    return orientation_argument_of_periapsis(
        *(vec4d*)major, *(vec4d*)minor, *(vec4d*)normal);
}
