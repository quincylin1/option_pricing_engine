#ifndef __MONTE_CARLO_DOUBLE_DIGITAL_HPP
#define __MONTE_CARLO_DOUBLE_DIGITAL_HPP

#include "../monte_carlo.hpp"
#include <vector>

class MonteCarloDoubleDigital : public MonteCarloBase
{
public:
    MonteCarloDoubleDigital(const unsigned long num_sims,
                      Option *pOption,
                      StatisticalDistribution *sd);

    virtual ~MonteCarloDoubleDigital();

    virtual double calc_S_drift(const double &S) const;
    virtual double calc_S_stochastic(const double &spot_normal) const;

    virtual double monte_carlo_sim(const double &S_0);
};

#endif