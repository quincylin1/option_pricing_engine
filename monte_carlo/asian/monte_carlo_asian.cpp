#ifndef __MONTE_CARLO_ASIAN_CPP
#define __MONTE_CARLO_ASIAN_CPP

#include "monte_carlo_asian.hpp"
#include <cmath>
#include <vector>

MonteCarloAsian::MonteCarloAsian(const unsigned long num_sims,
                                 const unsigned long num_intervals,
                                 AsianOption *pAsianOption,
                                 StatisticalDistribution *sd)
: _num_sims(num_sims), _num_intervals(num_intervals), _pAsianOption(pAsianOption), _sd(sd) {}
MonteCarloAsian::~MonteCarloAsian(){}

void MonteCarloAsian::calc_path_spot_prices(std::vector<double> &spot_prices) {
    double dt = _pAsianOption->_T / static_cast<double>(spot_prices.size());
    double drift = exp(dt * (_pAsianOption->_r - 0.5 * (_pAsianOption->_v) * (_pAsianOption->_v)));
    double vol = _pAsianOption->_v * sqrt(dt);

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

double MonteCarloAsian::monte_carlo_sim(const double& S_0) {

    double payoff_sum = 0.0;
    std::vector<double> spot_prices(_num_intervals, S_0);

    for (int i = 0; i < _num_sims; i++){
        calc_path_spot_prices(spot_prices);
        payoff_sum += _pAsianOption->pay_off_price(spot_prices);
    }

    double discount_payoff_avg = (payoff_sum / static_cast<double>(_num_sims)) * exp(-_pAsianOption->_r * _pAsianOption->_T);
    return discount_payoff_avg;
}

#endif 
