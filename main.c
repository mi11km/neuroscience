#include <stdio.h>
#include <math.h>

#include "pkg/numerical_calculation/ordinary_differential_equation.h"


int main() {
    const double t0 = 0.0;  // 初期値 t = 0.0
    const double x0 = 1.0;  // 初期値 x(t0) = 1.0
    const double t = 1.0;   // 求める数値 t = 1.0 のときの x(t)
    double v_euler, v_heun, v_runge_kutta;  // それぞれの方法の解を入れる変数
    for (int i = 0; i < 12; ++i) {
        const double delta_t = 1.0 / pow(2, i); // 刻み幅
        v_euler = euler_method(t0, x0, delta_t, t);
        v_heun = heun_method(t0, x0, delta_t, t);
        v_runge_kutta = runge_kutta_method(t0, x0, delta_t, t);
        printf("%16.15f %16.15f %16.15f %16.15f\n", delta_t, v_euler, v_heun, v_runge_kutta);
    }
    return 0;
}
