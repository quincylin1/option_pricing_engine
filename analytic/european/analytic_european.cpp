#ifndef __ANALYTIC_EUROPEAN_CPP
#define __ANALYTIC_EUROPEAN_CPP

#include "analytic_european.hpp"


AnalyticEuropean::AnalyticEuropean(){}
AnalyticEuropean::AnalyticEuropean(StatisticalDistribution *sd) : _sd(sd) {}
AnalyticEuropean::~AnalyticEuropean(){}

double AnalyticEuropean::d_j(const int &j, const double &S, const double &K, const double &r,
                             const double &v, const double &T, const double& q) const
{
    return (log(S / K) + (r - q + pow(-1, j - 1) * 0.5 * v * v) * T / (v * pow(T, 0.5)));
}

double AnalyticEuropean::calc_call_price(const double &S, const double &K, const double &r,
                                              const double &v, const double &T, const double& q) const
{
    return S * exp(-q * T) * _sd->cdf(d_j(1, S, K, r, v, T, q)) - K * exp(-r * T) * _sd->cdf(d_j(2, S, K, r, v, T, q));
}

double AnalyticEuropean::calc_put_price(const double &S, const double &K, const double &r,
                                             const double &v, const double &T, const double& q) const
{
    return -S * exp(-q * T) * _sd->cdf(-d_j(1, S, K, r, v, T, q)) + K * exp(-r * T) * _sd->cdf(-d_j(2, S, K, r, v, T, q));
}

#endif 