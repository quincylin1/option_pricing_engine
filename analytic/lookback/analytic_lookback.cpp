#ifndef __ANALYTIC_LOOKBACK_CPP
#define __ANALYTIC_LOOKBACK_CPP

#include "analytic_lookback.hpp"


AnalyticLookBack::AnalyticLookBack(StatisticalDistribution* sd, unsigned long num_intervals)
 : sd_(sd), num_intervals_(num_intervals) {}

double AnalyticLookBack::calc_a_1(const double &S, 
                             const double &H, 
                             const double &r, 
                             const double &v, 
                             const double &T) const
{
    return (log(S / H) + (r + 0.5 * v * v) * T) / (v * sqrt(T));
}

double AnalyticLookBack::calc_a_2(const double &S,
                             const double &H,
                             const double &r,
                             const double &v,
                             const double &T) const
{
    return calc_a_1(S, H, r, v, T) - v * sqrt(T);
}

double AnalyticLookBack::calc_a_3(const double &S,
                             const double &H,
                             const double &r,
                             const double &v,
                             const double &T) const
{
    return calc_a_1(S, H, r, v, T) - 2.0 * r * sqrt(T) /v;
}

double AnalyticLookBack::calc_call_price(const double &S, const double& m, const double &r,
                                  const double &v, const double &T) const
{
    double a_1 = calc_a_1(S, m, r, v, T);
    double a_2 = calc_a_2(S, m, r, v, T);
    double a_3 = calc_a_3(S, m, r, v, T);

    return S * sd_->cdf(a_1) 
      - m * exp(-r * T) * sd_->cdf(a_2) 
      - ((S * v * v) / (2 * r)) 
      * (sd_->cdf(-a_1) - exp(-r * T) * pow((m / S), ((2 * r) / (v * v))) * sd_->cdf(-a_3));
}

double AnalyticLookBack::calc_put_price(const double &S, const double &M, const double &r,
                                         const double &v, const double &T) const
{
    double a_1 = calc_a_1(S, M, r, v, T);
    double a_2 = calc_a_2(S, M, r, v, T);
    double a_3 = calc_a_3(S, M, r, v, T);

    return - S * sd_->cdf(-a_1) 
      + M * exp(-r * T) * sd_->cdf(-a_2) 
      + ((S * v * v) / (2 * r)) 
      * (sd_->cdf(a_1) - exp(-r * T) * pow((M / S), ((2 * r) / (v * v))) * sd_->cdf(a_3));
}

#endif 