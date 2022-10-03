#ifndef __STANDARD_GAMMA_HPP
#define __STANDARD_GAMMA_HPP

#include "statistics.hpp"

class StandardGammaDistribution : public StatisticalDistribution
{
public:
    StandardGammaDistribution();
    virtual ~StandardGammaDistribution();

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