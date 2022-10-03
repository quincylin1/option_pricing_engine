#ifndef __MC_FDM_CPP
#define __MC_FDM_CPP

#include "mc_fdm_greeks.hpp"

#include <vector>

MonteCarloFDMGreeks::MonteCarloFDMGreeks(){}
MonteCarloFDMGreeks::MonteCarloFDMGreeks(unsigned long num_sims,
                                         Option *pOption, PayOff *pPayOff,
                                         StatisticalDistribution *sd)
: _num_sims(num_sims), _pOption(pOption), _pPayOff(pPayOff), _sd(sd) {}
MonteCarloFDMGreeks::~MonteCarloFDMGreeks(){}

double MonteCarloFDMGreeks::calc_S_drift(const double S_0, const double r, const double v, const double T) const
{
    return S_0 * exp((r - 0.5 * v * v) * T);
}

double MonteCarloFDMGreeks::calc_S_stochastic(const double v, const double T, const double spot_normal) const
{
    return exp(v * pow(T, 0.5) * spot_normal);
}

void MonteCarloFDMGreeks::monte_carlo_sim(double S_0, const double r,
                                          const double v, const double T,
                                          const double K,
                                          const double delta_S,
                                          double &price_Sp,
                                          double &price_S,
                                          double &price_Sm)
{
    double Sp_drift = calc_S_drift(S_0 + delta_S, r, v, T);
    double S_drift = calc_S_drift(S_0, r, v, T);
    double Sm_drift = calc_S_drift(S_0 - delta_S, r, v, T);

    double Sp_cur = 0.0;
    double S_cur = 0.0;
    double Sm_cur = 0.0;

    double payoff_sum_p = 0.0;
    double payoff_sum = 0.0;
    double payoff_sum_m = 0.0;

    std::vector<double> spot_uniforms(_num_sims, 0.0);
    std::vector<double> spot_normals(_num_sims, 0.0);

    for (int i = 0; i < spot_uniforms.size(); i++)
    {
        spot_uniforms[i] = rand() / static_cast<double>(RAND_MAX);
    }

    _sd->random_draws(spot_uniforms, spot_normals);

    for (int i = 0; i < _num_sims; i++){
        double S_stochastic = calc_S_stochastic(v, T, spot_normals[i]);
        Sp_cur = Sp_drift * S_stochastic;
        S_cur = S_drift * S_stochastic;
        Sm_cur = Sm_drift * S_stochastic;

        payoff_sum_p += _pPayOff->operator()(Sp_cur);
        payoff_sum += _pPayOff->operator()(S_cur);
        payoff_sum_m += _pPayOff->operator()(Sm_cur);
    }

    price_Sp = (payoff_sum_p / static_cast<double>(_num_sims)) * exp(-r * T);
    price_S = payoff_sum / static_cast<double>(_num_sims) * exp(-r * T);
    price_Sm = payoff_sum_m / static_cast<double>(_num_sims) * exp(-r * T);
}

double MonteCarloFDMGreeks::calc_delta(const double S_0, const double delta_S) {

    double price_Sp = 0.0;
    double price_S = 0.0;
    double price_Sm = 0.0;

    monte_carlo_sim(S_0, _pOption->_r,
                    _pOption->_v, _pOption->_T,
                    _pOption->_K, delta_S, price_Sp, price_S, price_Sm);

    return (price_Sp - price_S) / delta_S;
}

double MonteCarloFDMGreeks::calc_gamma(const double S_0, const double delta_S)
{

    double price_Sp = 0.0;
    double price_S = 0.0;
    double price_Sm = 0.0;

    monte_carlo_sim(S_0, _pOption->_r,
                    _pOption->_v, _pOption->_T,
                    _pOption->_K, delta_S,
                    price_Sp, price_S, price_Sm);

    return (price_Sp - 2 * price_S + price_Sm) / (delta_S * delta_S);
}

#endif 