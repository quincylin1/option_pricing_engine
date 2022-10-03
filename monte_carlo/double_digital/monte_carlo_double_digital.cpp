#ifndef __MONTE_CARLO_DOUBLE_DIGITAL_CPP
#define __MONTE_CARLO_DOUBLE_DIGITAL_CPP

#include "monte_carlo_double_digital.hpp"
#include <iostream>

MonteCarloDoubleDigital::MonteCarloDoubleDigital(const unsigned long num_sims,
                                     Option *pOption,
                                     StatisticalDistribution *sd)
    : MonteCarloBase(num_sims, pOption, sd) {}
MonteCarloDoubleDigital::~MonteCarloDoubleDigital() {}

double MonteCarloDoubleDigital::calc_S_drift(const double &S) const
{
    return S * exp((_pOption->_r - 0.5 * (_pOption->_v) * (_pOption->_v)) * (_pOption->_T));
}

double MonteCarloDoubleDigital::calc_S_stochastic(const double &spot_normal) const
{
    return exp(_pOption->_v * pow(_pOption->_T, 0.5) * spot_normal);
}

double MonteCarloDoubleDigital::monte_carlo_sim(const double &S_0)
{
    double S_cur = 0.0;
    double payoff_sum = 0.0;

    std::vector<double> spot_uniforms(_num_sims, 0.0);
    std::vector<double> spot_normals(_num_sims, 0.0);

    for (int i = 0; i < spot_uniforms.size(); i++)
    {
        spot_uniforms[i] = rand() / static_cast<double>(RAND_MAX);
    }

    _sd->random_draws(spot_uniforms, spot_normals);
    double S_drift = calc_S_drift(S_0);

    for (int i = 0; i < _num_sims; i++)
    {
        double S_stochastic = calc_S_stochastic(spot_normals[i]);
        S_cur = S_drift * S_stochastic;
        payoff_sum += _pOption->_pPayOff->operator()(S_cur);
    }

    return payoff_sum / static_cast<double>(_num_sims);
}

#endif