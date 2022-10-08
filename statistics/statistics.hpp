#ifndef __STATISTIC_HPP
#define __STATISTIC_HPP

#include <cmath>
#include <vector>


class StatisticalDistribution{
public:
    StatisticalDistribution();
    virtual ~StatisticalDistribution();

    // Distribution functions
    virtual double pdf(const double &x) = 0;
    virtual double cdf(const double &x) = 0;

    // Inverse cdf (a.k.a quantile function)
    virtual double inv_cdf(const double &quantile) = 0;

    // Descriptive stats
    virtual double mean() const = 0;
    virtual double var() const = 0;
    virtual double stdev() const = 0;

    // Obtain a sequence of random numbers from this distribution
    virtual void random_draws(const std::vector<double> &uniform_draws,
                              std::vector<double> &dist_draws) = 0;
};

#endif