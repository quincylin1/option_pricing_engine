#ifndef __STANDARD_NORMAL_HPP
#define __STANDARD_NORMAL_HPP

#include "statistics.hpp"

class StandardNormalDistribution : public StatisticalDistribution
{
public:
    StandardNormalDistribution();
    virtual ~StandardNormalDistribution();

    virtual double pdf(const double &x) const;
    virtual double cdf(const double &x) const;

    virtual double inv_cdf(const double &quantile) const;

    virtual double mean() const;
    virtual double var() const;
    virtual double stdev() const;

    virtual void random_draws(const std::vector<double> &uniform_draws,
                              std::vector<double> &dist_draws);
};

#endif 