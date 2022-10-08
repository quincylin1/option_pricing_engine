#ifndef __STANDARD_GAMMA_HPP
#define __STANDARD_GAMMA_HPP

#include "statistics.hpp"
#include <complex>

class StandardGammaDistribution : public StatisticalDistribution
{
private:
    double alpha_;

public:
    StandardGammaDistribution(double alpha);
    virtual ~StandardGammaDistribution();

    virtual double gamma_function(double z);

    virtual double pdf(const double &x);
    virtual double cdf(const double &x);

    virtual double inv_cdf(const double &quantile);

    virtual double mean() const;
    virtual double var() const;
    virtual double stdev() const;

    virtual void random_draws(const std::vector<double> &uniform_draws,
                              std::vector<double> &dist_draws);
};

#endif