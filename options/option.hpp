#ifndef __OPTION_HPP
#define __OPTION_HPP

#include "../payoffs/payoff.hpp"


class Option{
public:
    double _K;
    double _r;
    double _T;
    double _v;
    double _q;

    PayOff *_pPayOff;

    Option();
    Option(double K,
           double r, double T,
           double v, double q,
           PayOff *pPayOff);

    virtual ~Option();
};




#endif