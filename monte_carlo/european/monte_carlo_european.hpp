#ifndef __MONTE_CARLO_EUROPEAN_HPP
#define __MONTE_CARLO_EUROPEAN_HPP

#include "../monte_carlo.hpp"
#include <vector>

class MonteCarloEuropean : public MonteCarloBase
{
public:
    MonteCarloEuropean(const unsigned long num_sims,
                       Option *pOption,
                       StatisticalDistribution *sd);

    virtual ~MonteCarloEuropean();

    virtual double calc_S_drift(const double &S) const;
    virtual double calc_S_stochastic(const double &spot_normal) const;

    virtual double monte_carlo_sim(const double& S_0);
};

#endif