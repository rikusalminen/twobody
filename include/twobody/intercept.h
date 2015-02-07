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

int intercept_search(
    const struct orbit *orbit1,
    const struct orbit *orbit2,
    double t0, double t1,
    double threshold,
    int max_steps,
    double *times,
    int max_times);

#endif
