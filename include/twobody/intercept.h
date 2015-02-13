#ifndef TWOBODY_INTERCEPT_H
#define TWOBODY_INTERCEPT_H

struct orbit;

int intercept_intersect(
    const struct orbit *orbit1,
    const struct orbit *orbit2,
    double threshold,
    double *fs);

int intercept_times(
    const struct orbit *orbit1,
    const struct orbit *orbit2,
    double t0, double t1,
    const double *fs,
    double *intercept_times,
    int max_times);

#include <twobody/simd4d.h> // XXX: need only vec4d typedef

struct intercept {
#ifndef TWOBODY_NO_SIMD
    vec4d position[2];
    vec4d velocity[2];
    vec4d relative_position;
    vec4d relative_velocity;
#else
    double position[8];
    double velocity[8];
    double relative_position[4];
    double relative_velocity[4];
#endif

    double mu;
    double time;
    double distance;
    double speed;

    double E1, E2;
    double xxx1, xxx2; // XXX: padding
};

double intercept_search(
    const struct orbit *orbit1,
    const struct orbit *orbit2,
    double t0, double t1,
    double threshold,
    double target_distance,
    int max_steps,
    struct intercept *intercept);

#endif
