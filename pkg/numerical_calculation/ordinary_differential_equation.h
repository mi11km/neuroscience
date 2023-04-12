#ifndef NEUROSCIENCE_ORDINARY_DIFFERENTIAL_EQUATION_H
#define NEUROSCIENCE_ORDINARY_DIFFERENTIAL_EQUATION_H

double euler_method(double t0, double x0, double delta_t, double t);

double heun_method(double t0, double x0, double delta_t, double t);

double runge_kutta_method(double t0, double x0, double delta_t, double t);

#endif //NEUROSCIENCE_ORDINARY_DIFFERENTIAL_EQUATION_H
