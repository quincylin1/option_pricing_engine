#ifndef __LOG_NORMAL_HPP
#define __LOG_NORMAL_HPP

#include "statistics.hpp"
#include "standard_normal.hpp"

class LogNormalDistribution : public StatisticalDistribution
{
private:
    StandardNormalDistribution *_snd;

public:
    LogNormalDistribution(StandardNormalDistribution *snd);
    virtual ~LogNormalDistribution();

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