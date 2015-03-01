#ifndef TWOBODY_FG_H
#define TWOBODY_FG_H

void fg(
    double x0, double y0,
    double xdot0, double ydot0,
    double x1, double y1,
    double xdot1, double ydot1,
    double *f, double *g,
    double *fdot, double *gdot);

#endif
