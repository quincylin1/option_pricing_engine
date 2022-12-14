#ifndef __ANALYTIC_LOOKBACK_HPP
#define __ANALYTIC_LOOKBACK_HPP

#include "../../statistics/statistics.hpp"
#include <vector>

class AnalyticLookBack
{
private:
    StatisticalDistribution *sd_;

    double calc_a_1(const double &S, // Spot price
                    const double &H, // Min/max of asset price over period
                    const double &r, // Risk free rate
                    const double &v, // Volatility of underlying asset
                    const double &T) const;

    double calc_a_2(const double &S, // Spot price
                    const double &H, // Min/max of asset price over period
                    const double &r, // Risk free rate
                    const double &v, // Volatility of underlying asset
                    const double &T) const;

    double calc_a_3(const double &S, // Spot price
                    const double &H, // Min/max of asset price over period
                    const double &r, // Risk free rate
                    const double &v, // Volatility of underlying asset
                    const double &T) const;

public:
    AnalyticLookBack(){};
    AnalyticLookBack(StatisticalDistribution *sd);
    ~AnalyticLookBack(){};

    double calc_call_price(const double &S_0, const double &m, const double &r,
                           const double &v, const double &T) const;

    double calc_put_price(const double &S_0, const double &M, const double &r,
                           const double &v, const double &T) const;
};

#endif 