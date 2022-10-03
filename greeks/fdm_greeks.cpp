#ifndef __FDM_GREEKS_CPP
#define __FDM_GREEKS_CPP

#include "fdm_greeks.h"

FDMGreeks::FDMGreeks() {}
FDMGreeks::FDMGreeks(
                     Option *pOption,
                     AnalyticEuropean *_pAnalytic,
                     StatisticalDistribution *sd)
: _pOption(_pOption), _pAnalytic(_pAnalytic), _sd(sd) {}
FDMGreeks::~FDMGreeks() {}

double FDMGreeks::calc_call_delta(const double &S_0, const double &delta_S) const
{
    return (_pAnalytic->calc_call_price(S_0 + delta_S, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T) 
      - _pAnalytic->calc_call_price(S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T)) / delta_S;
}

double FDMGreeks::calc_call_gamma(const double &S_0, const double &delta_S) const
{
    return (_pAnalytic->calc_call_price(S_0 + delta_S, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T) 
      - 2 * _pAnalytic->calc_call_price(S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T) + 
      _pAnalytic->calc_call_price(S_0 - delta_S, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T)) /
      (delta_S * delta_S);
}

double FDMGreeks::calc_call_vega(const double &S_0, const double &delta_v) const
{
    return (_pAnalytic->calc_call_price(S_0, _pOption->_K, _pOption->_r, _pOption->_v + delta_v, _pOption->_T) - _pAnalytic->calc_call_price(S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T)) / delta_v;
}

double FDMGreeks::calc_call_theta(const double &S_0, const double &delta_T) const
{
    return (_pAnalytic->calc_call_price(S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T + delta_T) - _pAnalytic->calc_call_price(S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T)) / delta_T;
}

double FDMGreeks::calc_call_rho(const double &S_0, const double &delta_r) const
{
    return (_pAnalytic->calc_call_price(S_0, _pOption->_K, _pOption->_r + delta_r, _pOption->_v, _pOption->_T) - _pAnalytic->calc_call_price(S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T)) / delta_r;
}

double FDMGreeks::calc_put_delta(const double &S_0, const double &delta_S) const
{
    return (_pAnalytic->calc_put_price(S_0 + delta_S, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T) - _pAnalytic->calc_put_price(S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T)) / delta_S;
}

double FDMGreeks::calc_put_gamma(const double &S_0, const double &delta_S) const
{
    return (_pAnalytic->calc_put_price(S_0 + delta_S, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T) - 2 * _pAnalytic->calc_put_price(S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T) +
            _pAnalytic->calc_put_price(S_0 - delta_S, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T)) /
           (delta_S * delta_S);
}

double FDMGreeks::calc_put_vega(const double &S_0, const double &delta_v) const
{
    return (_pAnalytic->calc_put_price(S_0, _pOption->_K, _pOption->_r, _pOption->_v + delta_v, _pOption->_T) - _pAnalytic->calc_put_price(S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T)) / delta_v;
}

double FDMGreeks::calc_put_theta(const double &S_0, const double &delta_T) const
{
    return (_pAnalytic->calc_put_price(S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T + delta_T) - _pAnalytic->calc_put_price(S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T)) / delta_T;
}

double FDMGreeks::calc_put_rho(const double &S_0, const double &delta_r) const
{
    return (_pAnalytic->calc_put_price(S_0, _pOption->_K, _pOption->_r + delta_r, _pOption->_v, _pOption->_T) - _pAnalytic->calc_put_price(S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T)) / delta_r;
}

#endif