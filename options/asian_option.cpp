#ifndef __ASIAN_OPTION_CPP
#define __ASIAN_OPTION_CPP

#include <numeric>
#include <cmath>
#include "asian_option.hpp"
#include <iostream>

AsianOption::AsianOption(double K,
                         double r, double T,
                         double v, double q,
                         PayOff *pPayOff)
    : Option(K, r, T, v, q, pPayOff) {}

// =====================
// AsianOptionArithmetic
// =====================
AsianOptionArithmetic::AsianOptionArithmetic(double K,
                                             double r, double T,
                                             double v, double q,
                                             PayOff *pPayOff)
    : AsianOption(K, r, T, v, q, pPayOff) {}

double AsianOptionArithmetic::pay_off_price(
    const std::vector<double> &spot_prices) const
{

    unsigned num_times = spot_prices.size();
    double sum = std::accumulate(spot_prices.begin(), spot_prices.end(), 0.0);
    double arith_mean = sum / static_cast<double>(num_times);
    return (*_pPayOff)(arith_mean);
}

AsianOptionGeometric::AsianOptionGeometric(double K,
                                           double r, double T,
                                           double v, double q,
                                           PayOff *pPayOff)
    : AsianOption(K, r, T, v, q, pPayOff) {}

double AsianOptionGeometric::pay_off_price(
    const std::vector<double> &spot_prices) const
{

    unsigned num_times = spot_prices.size();
    double log_sum = 0.0;

    for (int i = 0; i < spot_prices.size(); i++)
    {
        log_sum += log(spot_prices[i]);
    }
    double geom_mean = log_sum / static_cast<double>(num_times);
    return (*_pPayOff)(geom_mean);
}

#endif
