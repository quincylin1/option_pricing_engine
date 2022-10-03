#ifndef __EUROPEAN_OPTION_CPP
#define __EUROPEAN_OPTION_CPP

#include "european_option.hpp"

EuropeanOption::EuropeanOption(double K,
                               double r, double T,
                               double v, double q,
                               PayOff *pPayOff)
 : Option(K, r, T, v, q, pPayOff) {}

EuropeanOption::~EuropeanOption() {}

#endif