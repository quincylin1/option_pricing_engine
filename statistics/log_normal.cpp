#ifndef __LOG_NORMAL_CPP
#define __LOG_NORMAL_CPP

#include "log_normal.hpp"
#include <iostream>

LogNormalDistribution::LogNormalDistribution(StandardNormalDistribution *snd) : _snd(snd) {}
LogNormalDistribution::~LogNormalDistribution(){}

double LogNormalDistribution::pdf(const double &x) const {
    if (x <= 0){
        std::cout << "x must be greater than zero" << std::endl;
    }
    return (1 / (sqrt(2.0 * M_PI) * x)) * exp(-0.5 * (log(x) * log(x)));
}

double LogNormalDistribution::cdf(const double &x) const {
    return _snd->cdf(log(x));
}

double LogNormalDistribution::inv_cdf(const double &x) const {
    return exp(_snd->inv_cdf(x));
}

double LogNormalDistribution::mean() const { return exp(0.5); }

double LogNormalDistribution::var() const { return exp(1.0) * (exp(1.0) - 1.0); }

double LogNormalDistribution::stdev() const { return sqrt(exp(1.0) * (exp(1.0) - 1.0)); }

void LogNormalDistribution::random_draws(const std::vector<double> &uniform_draws,
                                         std::vector<double> &dist_draws) 
{
    // The simplest method to calculate this is with the Box−Muller method,
    // which has been used procedurally in many other chapters

    // Check that the uniform draws and dist draws are the same size and
    // have an even number of elements (necessary for B−M)
    if (uniform_draws.size() != dist_draws.size())
    {
        std::cout << "Draws vectors are of unequal size" << std::endl;
        return;
    }

    if (uniform_draws.size() % 2 != 0)
    {
        std::cout << "Uniform draw vector size not an even number" << std::endl;
        return;
    }

    for (int i = 0; i < uniform_draws.size() / 2; i++)
    {
        dist_draws[2 * i] = sqrt(-2.0 * log(uniform_draws[2 * i])) *
                            sin(2 * M_PI * uniform_draws[2 * i + 1]);
        dist_draws[2 * i + 1] = sqrt(-2.0 * log(uniform_draws[2 * i])) *
                                cos(2 * M_PI * uniform_draws[2 * i + 1]);
    }

    // Convert normal distribution to log-normal distribution
    for (int i = 0; i < dist_draws.size(); i++){
        dist_draws[i] = exp(dist_draws[i]);
    }
}

#endif 