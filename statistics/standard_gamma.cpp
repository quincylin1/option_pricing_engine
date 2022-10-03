#ifndef __STANDARD_GAMMA_CPP
#define __STANDARD_GAMMA_CPP

#include "standard_gamma.hpp"
#include <iostream>

StandardGammaDistribution::StandardGammaDistribution(){}
StandardGammaDistribution::~StandardGammaDistribution() {}

double StandardGammaDistribution::pdf(const double &x) const
{
}

double StandardGammaDistribution::cdf(const double &x) const
{
}

double StandardGammaDistribution::inv_cdf(const double &x) const
{
}

double StandardGammaDistribution::mean() const {}

double StandardGammaDistribution::var() const {}

double StandardGammaDistribution::stdev() const {}

void StandardGammaDistribution::random_draws(const std::vector<double> &uniform_draws,
                                         std::vector<double> &dist_draws)
{
}

#endif