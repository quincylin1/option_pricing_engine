#ifndef __EUROPEAN_OPTION_H
#define __EUROPEAN_OPTION_H

#include "../payoffs/payoff.h"
#include "option.hpp"

class EuropeanOption : public Option
{
public:
    EuropeanOption(double K,
                   double r, double T,
                   double v, double q,
                   PayOff *pPayOff);
    virtual ~EuropeanOption();
};
#endif