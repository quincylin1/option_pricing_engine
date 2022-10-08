#ifndef __ANALYTIC_EUROPEAN_HPP
#define __ANALYTIC_EUROPEAN_HPP

#include "../../statistics/statistics.hpp"


class AnalyticEuropean {
private:
    StatisticalDistribution *_sd;

public:
    AnalyticEuropean();
    AnalyticEuropean(StatisticalDistribution *sd);
    ~AnalyticEuropean();

    virtual double d_j(const int &j, const double &S, const double &K, const double &r,
                       const double &v, const double &T, const double& q) const;

    virtual double calc_call_price(const double &S, const double &K, const double &r,
                                   const double &v, const double &T, const double &q) const;
    virtual double calc_put_price(const double &S, const double &K, const double &r,
                                  const double &v, const double &T, const double &q) const;
};

#endif