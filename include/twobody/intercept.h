#ifndef TWOBODY_INTERCEPT_H
#define TWOBODY_INTERCEPT_H

struct orbit;

int intercept_intersect(
    const struct orbit *orbit1,
    const struct orbit *orbit2,
    double threshold,
    double *fs);

#endif
