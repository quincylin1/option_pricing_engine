#ifndef __ANALYTIC_GREEKS_CPP
#define __ANALYTIC_GREEKS_CPP

#include "analytic_greeks.hpp"


AnalyticGreeks::AnalyticGreeks(){}
AnalyticGreeks::AnalyticGreeks(Option *_pOption,
                               AnalyticEuropean *_pAnalytic,
                               StatisticalDistribution *_sd) 
: _pOption(_pOption), _pAnalytic(_pAnalytic), _sd(_sd) {}
AnalyticGreeks::~AnalyticGreeks(){}

double AnalyticGreeks::calc_call_delta(const double& S_0) const {
    return _sd->cdf(_pAnalytic->d_j(1, S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T));
}

double AnalyticGreeks::calc_call_gamma(const double &S_0) const
{
    return _sd->pdf(_pAnalytic->d_j(1, S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T))
     / (S_0 * _pOption->_v * sqrt(_pOption->_T));
}

double AnalyticGreeks::calc_call_vega(const double &S_0) const
{
    return S_0 * _sd->pdf(_pAnalytic->d_j(1, S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T)) * sqrt(_pOption->_T);
}

double AnalyticGreeks::calc_call_theta(const double &S_0) const
{
    return -(S_0 * _sd->pdf(_pAnalytic->d_j(1, S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T)) * _pOption->_v) / (2.0 * sqrt(_pOption->_T)) 
      - _pOption->_r * _pOption->_K * exp(-_pOption->_r * _pOption->_T) * _sd->cdf(_pAnalytic->d_j(2, S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T));
}

double AnalyticGreeks::calc_call_rho(const double &S_0) const
{
    return _pOption->_K * _pOption->_T * exp(-_pOption->_r * _pOption->_T) * _sd->cdf(_pAnalytic->d_j(2, S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T));
}

double AnalyticGreeks::calc_put_delta(const double &S_0) const
{
    return _sd->cdf(_pAnalytic->d_j(1, S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T)) - 1.0;
}

double AnalyticGreeks::calc_put_gamma(const double &S_0) const
{
    return _sd->pdf(_pAnalytic->d_j(1, S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T)) / (S_0 * _pOption->_v * sqrt(_pOption->_T));
}

double AnalyticGreeks::calc_put_vega(const double &S_0) const
{
    return S_0 * _sd->pdf(_pAnalytic->d_j(1, S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T)) * sqrt(_pOption->_T);
}

double AnalyticGreeks::calc_put_theta(const double &S_0) const
{
    return -(S_0 * _sd->pdf(_pAnalytic->d_j(1, S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T)) * _pOption->_v) / (2.0 * sqrt(_pOption->_T)) 
      + _pOption->_r * _pOption->_K * exp(-_pOption->_r * _pOption->_T) * _sd->cdf(-_pAnalytic->d_j(2, S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T));
}

double AnalyticGreeks::calc_put_rho(const double &S_0) const
{
    return -_pOption->_K * _pOption->_T * exp(-_pOption->_r * _pOption->_T) * _sd->cdf(-_pAnalytic->d_j(2, S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T));
}

#endif 