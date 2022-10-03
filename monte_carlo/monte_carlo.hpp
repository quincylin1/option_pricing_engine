#ifndef __MONTE_CARLO_HPP
#define __MONTE_CARLO_HPP

#include "../options/option.hpp"
#include "../statistics/statistics.hpp"

class MonteCarloBase {
protected:
    unsigned long _num_sims;
    Option *_pOption;
    StatisticalDistribution *_sd;

public:
    MonteCarloBase();
    MonteCarloBase(const unsigned long num_sims,
                   Option *pOption,
                   StatisticalDistribution *sd);
    virtual ~MonteCarloBase();

    virtual double monte_carlo_sim(const double& S_0) = 0;
};

#endif