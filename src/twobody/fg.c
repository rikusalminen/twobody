#include <twobody/fg.h>

void fg(
    double x0, double y0,
    double xdot0, double ydot0,
    double x1, double y1,
    double xdot1, double ydot1,
    double *f, double *g,
    double *fdot, double *gdot) {
    // specific angular momentum
    double h = x0*ydot0 - y0*xdot0;

    *f = (x1*ydot0 - xdot0*y1) / h;
    *g = (x0*y1 - x1*y0) / h;

    *fdot = (xdot1*ydot0 - xdot0*ydot1) / h;
    *gdot = (x0*ydot1 - xdot1*y0) / h;
}
