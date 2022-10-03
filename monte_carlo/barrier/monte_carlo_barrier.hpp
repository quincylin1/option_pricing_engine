#ifndef __MONTE_CARLO_BARRIER_HPP
#define __MONTE_CARLO_BARRIER_HPP

#include "../monte_carlo.hpp"
#include "../../options/barrier_option.hpp"
#include <vector>

class MonteCarloBarrier : public MonteCarloBase
{
private:
    unsigned long _num_sims;
    unsigned long _num_intervals;
    BarrierOption *_pBarrierOption;
    StatisticalDistribution *_sd;

public:
    MonteCarloBarrier(const unsigned long num_sims,
                    const unsigned long num_intervals,
                    BarrierOption *pBarrierOption,
                    StatisticalDistribution *sd);

    virtual ~MonteCarloBarrier();

    virtual void calc_path_spot_prices(std::vector<double> &spot_prices);

    virtual double monte_carlo_sim(const double& S_0);
};

#endif