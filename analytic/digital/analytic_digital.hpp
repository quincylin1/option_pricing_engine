#ifndef __ANALYTIC_DIGITAL_HPP
#define __ANALYTIC_DIGITAL_HPP

#include "../../statistics/statistics.hpp"

class AnalyticDigital
{
private:
    StatisticalDistribution *sd_;

    double d_j(const int &j, const double &S,
               const double &K, const double &r,
               const double &v, const double &T, const double &q) const;

public:
    AnalyticDigital(){};
    AnalyticDigital(StatisticalDistribution *sd);
    ~AnalyticDigital(){};

    double calc_call_price(const double &S, const double &K, const double &r,
                           const double &v, const double &T, const double &q) const;
    double calc_put_price(const double &S, const double &K, const double &r,
                          const double &v, const double &T, const double& q) const;
};

#endif