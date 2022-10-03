#ifndef __AMERICAN_OPTION_CPP
#define __AMERICAN_OPTION_CPP

#include "american_option.hpp"

AmericanOption::AmericanOption(double K,
                               double r, double T,
                               double v, double q,
                               PayOff *pPayOff)
: Option(K, r, T, v, q, pPayOff) {}

AmericanOption::~AmericanOption() {}

#endif