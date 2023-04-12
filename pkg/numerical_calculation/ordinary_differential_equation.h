#ifndef NEUROSCIENCE_ORDINARY_DIFFERENTIAL_EQUATION_H
#define NEUROSCIENCE_ORDINARY_DIFFERENTIAL_EQUATION_H

typedef double (*f_ptr)(double, double);

/*
 * ODE_params is the ordinary differential equation parameters.
 *   f: dx/dt = f(x, t)
 *   t0: 初期値 x(t0) = x0
 *   x0: 初期値 x(t0) = x0
 *   delta_t: 刻み幅
 *   t: 求めたい値 x(t)
 */
typedef struct {
    f_ptr f;
    double t0, x0, delta_t;
    double t;
} ODE_params;

double euler_method(ODE_params *p);

double heun_method(ODE_params *p);

double runge_kutta_method(ODE_params *p);

#endif //NEUROSCIENCE_ORDINARY_DIFFERENTIAL_EQUATION_H
