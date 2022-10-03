#ifndef __MONTE_CARLO_BARRIER_CPP
#define __MONTE_CARLO_BARRIER_CPP

#include "monte_carlo_barrier.hpp"
#include <cmath>
#include <vector>

MonteCarloBarrier::MonteCarloBarrier(const unsigned long num_sims,
                                 const unsigned long num_intervals,
                                 BarrierOption *pBarrierOption,
                                 StatisticalDistribution *sd)
    : _num_sims(num_sims), _num_intervals(num_intervals), _pBarrierOption(pBarrierOption), _sd(sd) {}
MonteCarloBarrier::~MonteCarloBarrier() {}


void MonteCarloBarrier::calc_path_spot_prices(std::vector<double> &spot_prices)
{
    double dt = _pBarrierOption->_T / static_cast<double>(spot_prices.size());
    double drift = exp(dt * (_pBarrierOption->_r - 0.5 * (_pBarrierOption->_v) * (_pBarrierOption->_v)));
    double vol = _pBarrierOption->_v * sqrt(dt);

    std::vector<double> spot_uniforms(_num_intervals, 0.0);
    std::vector<double> spot_normals(_num_intervals, 0.0);

    for (int i = 0; i < spot_uniforms.size(); i++)
    {
        spot_uniforms[i] = rand() / static_cast<double>(RAND_MAX);
    }

    _sd->random_draws(spot_uniforms, spot_normals);

    for (int i = 1; i < spot_prices.size(); i++)
    {
        spot_prices[i] = spot_prices[i - 1] * drift * exp(vol * spot_normals[i]);
    }
}

double MonteCarloBarrier::monte_carlo_sim(const double& S_0)
{

    double payoff_sum = 0.0;
    std::vector<double> spot_prices(_num_intervals, S_0);

    for (int i = 0; i < _num_sims; i++)
    {
        calc_path_spot_prices(spot_prices);
        payoff_sum += _pBarrierOption->calc_payoff(spot_prices);
    }

    double discount_payoff_avg = (payoff_sum / static_cast<double>(_num_sims)) * exp(-_pBarrierOption->_r * _pBarrierOption->_T);
    return discount_payoff_avg;
}

#endif
