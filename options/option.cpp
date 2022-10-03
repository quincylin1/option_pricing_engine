#ifndef __OPTION_CPP
#define __OPTION_CPP

#include "option.hpp"

Option::Option(){}
Option::Option(double K,
               double r, double T,
               double v, double q,
               PayOff *pPayOff) 
: _K(K), _r(r), _T(T), _v(v), _q(q), _pPayOff(pPayOff) {}

Option::~Option(){}

#endif 