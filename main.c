#include <stdio.h>

#include "pkg/neuron_model/hodgkin_huxley.h"


#define DT ( 0.01 )   // 10 micro s
#define T  ( 1000 )   // 1000 ms; unused
#define NT ( T / DT ) // T / DT


int main(void) {
    double v = E_REST;
    double m = m0(v);
    double h = h0(v);
    double n = n0(v);

    double i_ext = 9.0; // micro A / cm^2

    for (int32_t nt = 0; nt < NT; nt++) {
        double t = DT * nt;
        printf("%f %f %f %f %f\n", t, v, m, h, n);

        double dmdt1 = dmdt(v, m);
        double dhdt1 = dhdt(v, h);
        double dndt1 = dndt(v, n);
        double dvdt1 = dvdt(v, m, h, n, i_ext);

        double dmdt2 = dmdt(v + .5 * DT * dvdt1, m + .5 * DT * dmdt1);
        double dhdt2 = dhdt(v + .5 * DT * dvdt1, h + .5 * DT * dhdt1);
        double dndt2 = dndt(v + .5 * DT * dvdt1, n + .5 * DT * dndt1);
        double dvdt2 = dvdt(v + .5 * DT * dvdt1, m + .5 * DT * dmdt1, h + .5 * DT * dhdt1, n + .5 * DT * dndt1, i_ext);

        double dmdt3 = dmdt(v + .5 * DT * dvdt2, m + .5 * DT * dmdt2);
        double dhdt3 = dhdt(v + .5 * DT * dvdt2, h + .5 * DT * dhdt2);
        double dndt3 = dndt(v + .5 * DT * dvdt2, n + .5 * DT * dndt2);
        double dvdt3 = dvdt(v + .5 * DT * dvdt2, m + .5 * DT * dmdt2, h + .5 * DT * dhdt2, n + .5 * DT * dndt2, i_ext);

        double dmdt4 = dmdt(v + DT * dvdt3, m + DT * dmdt3);
        double dhdt4 = dhdt(v + DT * dvdt3, h + DT * dhdt3);
        double dndt4 = dndt(v + DT * dvdt3, n + DT * dndt3);
        double dvdt4 = dvdt(v + DT * dvdt3, m + DT * dmdt3, h + DT * dhdt3, n + DT * dndt3, i_ext);

        m += DT * (dmdt1 + 2 * dmdt2 + 2 * dmdt3 + dmdt4) / 6.;
        h += DT * (dhdt1 + 2 * dhdt2 + 2 * dhdt3 + dhdt4) / 6.;
        n += DT * (dndt1 + 2 * dndt2 + 2 * dndt3 + dndt4) / 6.;
        v += DT * (dvdt1 + 2 * dvdt2 + 2 * dvdt3 + dvdt4) / 6.;
    }
}