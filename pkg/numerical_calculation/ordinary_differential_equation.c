#include <stdio.h>
#include <stdlib.h>
#include "ordinary_differential_equation.h"


void validate_args(double t0, double t) {
    if (t0 > t) {
        fprintf(stderr, "invalid args. t0 < t must hold. t0: %16.15f, t: %16.15f", t0, t);
        exit(1);
    }
}

/*
 * euler_method returns x(p->t) value resolving dx/dt = p->f(x, t)
 * the initial condition is x(p->t0) = p->x0
 * */
double euler_method(ODE_params *p) {
    validate_args(p->t0, p->t);
    double t = p->t0, x = p->x0;
    while (t < p->t) {
        x += p->delta_t * p->f(x, t);
        t += p->delta_t;
    }
    return x;
}

/*
 * heun_method returns x(p->t) value resolving dx/dt = p->f(x, t)
 * the initial condition is x(p->t0) = p->x0
 * */
double heun_method(ODE_params *p) {
    validate_args(p->t0, p->t);
    double t = p->t0, x = p->x0;
    double k1, k2;
    while (t < p->t) {
        k1 = p->delta_t * p->f(x, p->t);
        k2 = p->delta_t * p->f(x + k1, p->t + p->delta_t);
        x += (k1 + k2) / 2.0;
        t += p->delta_t;
    }
    return x;
}

/*
 * runge_kutta_method returns x(p->t) value resolving dx/dt = p->f(x, p->t)
 * the initial condition is x(p->t0) = p->x0
 * */
double runge_kutta_method(ODE_params *p) {
    validate_args(p->t0, p->t);
    double t = p->t0, x = p->x0;
    double k1, k2, k3, k4;
    while (t < p->t) {
        k1 = p->delta_t * p->f(x, p->t);
        k2 = p->delta_t * p->f(x + k1 / 2.0, p->t + p->delta_t / 2.0);
        k3 = p->delta_t * p->f(x + k2 / 2.0, p->t + p->delta_t / 2.0);
        k4 = p->delta_t * p->f(x + k3, p->t + p->delta_t);
        x += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        t += p->delta_t;
    }
    return x;
}