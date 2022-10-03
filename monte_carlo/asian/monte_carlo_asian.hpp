#ifndef __MONTE_CARLO_ASIAN_HPP
#define __MONTE_CARLO_ASIAN_HPP

#include "../monte_carlo.hpp"
#include "../../options/asian_option.hpp"
#include <vector>

class MonteCarloAsian : public MonteCarloBase
{
private:
    unsigned long _num_sims;
    unsigned long _num_intervals;
    AsianOption *_pAsianOption;
    StatisticalDistribution *_sd;

public:
    MonteCarloAsian(const unsigned long _num_sims,
                    const unsigned long _num_intervals,
                    AsianOption *_pAsianOption,
                    StatisticalDistribution *_sd);

    virtual ~MonteCarloAsian();

    virtual void calc_path_spot_prices(std::vector<double> &spot_prices);

    virtual double monte_carlo_sim(const double& S_0);
};

#endif