#ifndef __ANALYTIC_DIGITAL_CPP
#define __ANALYTIC_DIGITAL_CPP

#include "analytic_digital.hpp"

AnalyticDigital::AnalyticDigital(StatisticalDistribution *sd) : sd_(sd) {}

double AnalyticDigital::d_j(const int &j, const double &S, 
                            const double &K, const double &r,
                            const double &v, const double &T) const
{
    return (log(S / K) + (r + pow(-1, j - 1) * 0.5 * v * v) * T / (v * pow(T, 0.5)));
}

double AnalyticDigital::calc_call_price(const double &S, const double &K, const double &r,
                                        const double &v, const double &T) const
{
    return exp(-r * T) * sd_->cdf(d_j(2, S, K, r, v, T));
}

double AnalyticDigital::calc_put_price(const double &S, const double &K, const double &r,
                                       const double &v, const double &T) const
{
    return exp(-r * T) * sd_->cdf(-d_j(2, S, K, r, v, T));
}

#endif 