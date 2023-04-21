#include <math.h>

#include "hodgkin_huxley.h"


#define C      (   1.0 )              // micro F / cm^2
#define G_LEAK (   0.3 )              // mS / cm^2
#define E_LEAK (  10.6 + ( E_REST ) ) // mV
#define G_NA   ( 120.0 )              // mS / cm^2
#define E_NA   ( 115.0 + ( E_REST ) ) // mV
#define G_K    (  36.0 )              // mS / cm^2
#define E_K    ( -12.0 + ( E_REST ) ) // mV


static inline double alpha_m(const double v) {
    return (2.5 - 0.1 * (v - E_REST)) / (exp(2.5 - 0.1 * (v - E_REST)) - 1.0);
}

static inline double beta_m(const double v) {
    return 4.0 * exp(-(v - E_REST) / 18.0);
}

static inline double alpha_h(const double v) {
    return 0.07 * exp(-(v - E_REST) / 20.0);
}

static inline double beta_h(const double v) {
    return 1.0 / (exp(3.0 - 0.1 * (v - E_REST)) + 1.0);
}

static inline double alpha_n(const double v) {
    return (0.1 - 0.01 * (v - E_REST)) / (exp(1 - 0.1 * (v - E_REST)) - 1.0);
}

static inline double beta_n(const double v) {
    return 0.125 * exp(-(v - E_REST) / 80.0);
}


double m0(const double v) { return alpha_m(v) / (alpha_m(v) + beta_m(v)); }

double h0(const double v) { return alpha_h(v) / (alpha_h(v) + beta_h(v)); }

double n0(const double v) { return alpha_n(v) / (alpha_n(v) + beta_n(v)); }

static inline double tau_m(const double v) { return 1. / (alpha_m(v) + beta_m(v)); }

static inline double tau_h(const double v) { return 1. / (alpha_h(v) + beta_h(v)); }

static inline double tau_n(const double v) { return 1. / (alpha_n(v) + beta_n(v)); }


double dmdt(const double v, const double m) { return (1.0 / tau_m(v)) * (-m + m0(v)); }

double dhdt(const double v, const double h) { return (1.0 / tau_h(v)) * (-h + h0(v)); }

double dndt(const double v, const double n) { return (1.0 / tau_n(v)) * (-n + n0(v)); }

double dvdt(const double v, const double m, const double h, const double n, const double i_ext) {
    return (-G_LEAK * (v - E_LEAK) - G_NA * m * m * m * h * (v - E_NA) - G_K * n * n * n * n * (v - E_K) + i_ext) / C;
}