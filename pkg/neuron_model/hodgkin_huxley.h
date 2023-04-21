#ifndef NEUROSCIENCE_HODGKIN_HUXLEY_H
#define NEUROSCIENCE_HODGKIN_HUXLEY_H

#define E_REST ( -65.0 )              // mV

double m0(double v);

double h0(double v);

double n0(double v);


double dmdt(double v, double m);

double dhdt(double v, double h);

double dndt(double v, double n);

double dvdt(double v, double m, double h, double n, double i_ext);

#endif //NEUROSCIENCE_HODGKIN_HUXLEY_H
