#ifndef __HESTON_EULER_CPP
#define __HESTON_EULER_CPP

#include "heston_euler.hpp"
#include "../statistics/correlated_snd.h"
#include "../statistics/statistics.h"


HestonEuler::HestonEuler(){}
HestonEuler::HestonEuler(unsigned num_sims,
                         unsigned num_intervals,
                         Option* pOption,
                         double kappa, double theta,
                         double xi, double rho) : 
_num_sims(num_sims), _num_intervals(num_intervals), _pOption(pOption),
_kappa(kappa), _theta(theta), _xi(xi), _rho(rho) {}

HestonEuler::~HestonEuler(){}

void HestonEuler::generate_normal_correlated_paths(std::vector<double> &spot_normals,
                                                    std::vector<double> &cor_normals)
{

    StandardNormalDistribution snd;
    std::vector<double> snd_uniform_draws(_num_intervals, 0.0);

    // Simple random number generator based on RAND
    for (int i = 0; i < snd_uniform_draws.size(); i++)
    {
        snd_uniform_draws[i] = rand() / static_cast<double>(RAND_MAX);
    }

    // Create standard normal random draws
    snd.random_draws(snd_uniform_draws, spot_normals);

    CorrelatedSND csnd(_rho, &spot_normals);

    std::vector<double> csnd_uniform_draws(_num_intervals, 0.0);

    // Uniform generation of correlated SND
    for (int i = 0; i < csnd_uniform_draws.size(); i++)
    {
        csnd_uniform_draws[i] = rand() / static_cast<double>(RAND_MAX);
    }

    // Convert uniform draws to correlated SND
    csnd.random_draws(csnd_uniform_draws, cor_normals);
}

double HestonEuler::monte_carlo_sim(const double& S_0)
{

    std::vector<double> spot_draws(_num_intervals, 0.0);  // Vector of initial spot normal draws
    std::vector<double> vol_draws(_num_intervals, 0.0);   // Vector of initial correlated vol normal draws
    std::vector<double> spot_prices(_num_intervals, S_0); // Vector of initial spot prices
    std::vector<double> vol_prices(_num_intervals, _pOption->_v);  // Vector of initial vol prices

    double payoff_sum = 0.0;

    for (unsigned i = 0; i < _num_sims; i++)
    {
        generate_normal_correlated_paths(spot_draws, vol_draws);
        calculate_vol_path(vol_draws, vol_prices);
        calculate_spot_path(spot_draws, vol_prices, spot_prices);

        payoff_sum += _pOption->_pPayOff->operator()(spot_prices[_num_intervals - 1]);
    }

    double option_price = (payoff_sum / static_cast<double>(_num_sims)) * exp(-_pOption->_r * _pOption->_T);
}

HestonEulerReflection::HestonEulerReflection(unsigned num_sims,
                                             unsigned num_intervals,
                                             Option *pOption,
                                             double kappa, double theta,
                                             double xi, double rho)
 : HestonEuler(num_sims, num_intervals, pOption, kappa, theta, xi, rho) {}
HestonEulerReflection::~HestonEulerReflection() {}

void HestonEulerReflection::calculate_vol_path(const std::vector<double> &vol_draws,
                                                       std::vector<double> &vol_path)
{

    double dt = _pOption->_T / static_cast<double>(_num_intervals);

    for (int i = 1; i < _num_intervals; i++)
    {
        vol_path[i] = std::fabs(vol_path[i - 1]) + _kappa * dt * (_theta - std::fabs(vol_path[i - 1])) +
                      _xi * sqrt(std::fabs(vol_path[i - 1]) * dt) * vol_draws[i - 1];
    }
}

void HestonEulerReflection::calculate_spot_path(const std::vector<double> &spot_draws,
                                                        const std::vector<double> &vol_path,
                                                        std::vector<double> &spot_path)
{

    double dt = _pOption->_T / static_cast<double>(_num_intervals);

    for (int i = 1; i < _num_intervals; i++)
    {
        spot_path[i] = spot_path[i - 1] * exp((_pOption->_r - 0.5 * std::fabs(vol_path[i])) * dt +
                                              sqrt(std::fabs(vol_path[i]) * dt) * spot_draws[i - 1]);
    }
}

HestonEulerPartialTruncation::HestonEulerPartialTruncation(unsigned num_sims,
                                                           unsigned num_intervals,
                                                           Option *pOption,
                                                           double kappa, double theta,
                                                           double xi, double rho)
    : HestonEuler(num_sims, num_intervals, pOption, kappa, theta, xi, rho) {}
HestonEulerPartialTruncation::~HestonEulerPartialTruncation() {}

void HestonEulerPartialTruncation::calculate_vol_path(const std::vector<double> &vol_draws,
                                                              std::vector<double> &vol_path)
{

    double dt = _pOption->_T / static_cast<double>(_num_intervals);

    for (int i = 1; i < _num_intervals; i++)
    {
        vol_path[i] = vol_path[i - 1] + _kappa * dt * (_theta - vol_path[i - 1]) +
                      _xi * sqrt(vol_path[i - 1] * dt) * vol_draws[i - 1];
    }
}

void HestonEulerPartialTruncation::calculate_spot_path(const std::vector<double> &spot_draws,
                                                               const std::vector<double> &vol_path,
                                                               std::vector<double> &spot_path)
{

    double dt = _pOption->_T / static_cast<double>(_num_intervals);

    for (int i = 1; i < _num_intervals; i++)
    {
        spot_path[i] = spot_path[i - 1] * exp((_pOption->_r - 0.5 * vol_path[i]) * dt +
                                              sqrt(vol_path[i] * dt) * spot_draws[i - 1]);
    }
}

HestonEulerFullTruncation::HestonEulerFullTruncation(unsigned num_sims,
                                                     unsigned num_intervals,
                                                     Option *pOption,
                                                     double kappa, double theta,
                                                     double xi, double rho)
    : HestonEuler(num_sims, num_intervals, pOption, kappa, theta, xi, rho) {}
HestonEulerFullTruncation::~HestonEulerFullTruncation() {}

void HestonEulerFullTruncation::calculate_vol_path(const std::vector<double> &vol_draws,
                                                           std::vector<double> &vol_path)
{

    double dt = _pOption->_T / static_cast<double>(_num_intervals);

    for (int i = 1; i < _num_intervals; i++)
    {
        double v_max = std::max(vol_path[i - 1], 0.0);

        vol_path[i] = vol_path[i - 1] + _kappa * dt * (_theta - v_max) +
                      _xi * sqrt(v_max * dt) * vol_draws[i - 1];
    }
}

void HestonEulerFullTruncation::calculate_spot_path(const std::vector<double> &spot_draws,
                                                            const std::vector<double> &vol_path,
                                                            std::vector<double> &spot_path)
{

    double dt = _pOption->_T / static_cast<double>(_num_intervals);

    for (int i = 1; i < _num_intervals; i++)
    {
        double v_max = std::max(vol_path[i], 0.0);
        spot_path[i] = spot_path[i - 1] * exp((_pOption->_r - 0.5 * v_max) * dt +
                                              sqrt(v_max * dt) * spot_draws[i - 1]);
    }
}

#endif 