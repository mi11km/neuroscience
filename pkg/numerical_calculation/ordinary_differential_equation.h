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

/* example usage

#include <stdio.h>
#include <math.h>
#include "pkg/numerical_calculation/ordinary_differential_equation.h"

double f(double x, double t) {
    // ex. dx/dt = x
    return x;
}

int main(void) {
    // 以下の1階常微分方程式の数値計算
    // dx/dt = f(x, t)
    // x(t0) = x0
    ODE_params p = {f, 0.0, 1.0, 1.0, 1.0,};
    double v_euler, v_heun, v_runge_kutta;  // それぞれの方法の解を入れる変数
    for (int i = 0; i < 12; ++i) {
        p.delta_t = 1.0 / pow(2, i); // 刻み幅
        v_euler = euler_method(&p);
        v_heun = heun_method(&p);
        v_runge_kutta = runge_kutta_method(&p);
        printf("%16.15f %16.15f %16.15f %16.15f\n", p.delta_t, v_euler, v_heun, v_runge_kutta);
    }
    return 0;
}

*/
