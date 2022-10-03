#ifndef __AMERICAN_OPTION_HPP
#define __AMERICAN_OPTION_HPP

#include "../payoffs/payoff.h"
#include "option.hpp"

class AmericanOption : public Option
{
public:

    AmericanOption(double K,
                   double r, double T,
                   double v, double q,
                   PayOff *pPayOff);
    virtual ~AmericanOption();
};
#endif