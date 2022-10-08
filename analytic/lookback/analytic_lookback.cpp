#ifndef __ANALYTIC_LOOKBACK_CPP
#define __ANALYTIC_LOOKBACK_CPP

#include "analytic_lookback.hpp"


AnalyticLookBack::AnalyticLookBack(StatisticalDistribution* sd) : sd_(sd) {}

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
    return calc_a_1(S, H, r, v, T) - (2.0 * r * sqrt(T) / v);
}

double AnalyticLookBack::calc_call_price(const double &S, const double& m, const double &r,
                                  const double &v, const double &T) const
{
    double a_1 = calc_a_1(S, m, r, v, T);
    double a_2 = calc_a_2(S, m, r, v, T);
    double a_3 = calc_a_3(S, m, r, v, T);

    double term_1 = S * sd_->cdf(a_1);
    double term_2 = m * exp(-r * T) * sd_->cdf(a_2);
    double term_3 = ((S * v * v) / (2.0 * r)) * (sd_->cdf(-a_1) - exp(-r * T) * pow((m / S), ((2 * r) / (v * v))) * sd_->cdf(-a_3));

    return term_1 - term_2 - term_3;
}

double AnalyticLookBack::calc_put_price(const double &S, const double &M, const double &r,
                                         const double &v, const double &T) const
{
    double a_1 = calc_a_1(S, M, r, v, T);
    double a_2 = calc_a_2(S, M, r, v, T);
    double a_3 = calc_a_3(S, M, r, v, T);

    double term_1 = S * sd_->cdf(-a_1);
    double term_2 = M * exp(-r * T) * sd_->cdf(-a_2);
    double term_3 = ((S * v * v) / (2 * r)) * (sd_->cdf(a_1) - exp(-r * T) * pow((M / S), ((2 * r) / (v * v))) * sd_->cdf(a_3));

    return - term_1 + term_2 + term_3;
}

#endif 