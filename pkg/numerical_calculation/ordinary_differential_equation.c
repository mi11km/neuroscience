#include <stdio.h>
#include <stdlib.h>

/*
 * f is used as dx/dt = f(x, t) in the following methods for resoling ordinary differential equation
 */
double f(double x, double t) {
    // ex. dx/dt = x
    return x;  // TODO 可変にする（今だと dx/dt = x 固定になっている）
}

void validate_args(double t0, double t) {
    if (t0 > t) {
        fprintf(stderr, "invalid args. t0 < t must hold. t0: %16.15f, t: %16.15f", t0, t);
        exit(1);
    }
}

/*
 * euler_method returns x(t) value resolving dx/dt = f(x, t), whose initial condition is x(t0) = x0
 * */
double euler_method(double t0, double x0, double delta_t, double t) {
    validate_args(t0, t);
    double tt = t0, x = x0;
    while (tt < t) {
        x += delta_t * f(x, tt);
        tt += delta_t;
    }
    return x;
}

/*
 * heun_method returns x(t) value resolving dx/dt = f(x, t), whose initial condition is x(t0) = x0
 * */
double heun_method(double t0, double x0, double delta_t, double t) {
    validate_args(t0, t);
    double tt = t0, x = x0;
    double k1, k2;
    while (tt < t) {
        k1 = delta_t * f(x, t);
        k2 = delta_t * f(x + k1, t + delta_t);
        x += (k1 + k2) / 2.0;
        tt += delta_t;
    }
    return x;
}

/*
 * runge_kutta_method returns x(t) value resolving dx/dt = f(x, t), whose initial condition is x(t0) = x0
 * */
double runge_kutta_method(double t0, double x0, double delta_t, double t) {
    validate_args(t0, t);
    double tt = t0, x = x0;
    double k1, k2, k3, k4;
    while (tt < t) {
        k1 = delta_t * f(x, t);
        k2 = delta_t * f(x + k1 / 2.0, t + delta_t / 2.0);
        k3 = delta_t * f(x + k2 / 2.0, t + delta_t / 2.0);
        k4 = delta_t * f(x + k3, t + delta_t);
        x += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        tt += delta_t;
    }
    return x;
}