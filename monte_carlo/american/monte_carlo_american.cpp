#ifndef __MONTE_CARLO_AMERICAN_CPP
#define __MONTE_CARLO_AMERICAN_CPP

#include "monte_carlo_american.hpp"
#include "polyfit.hpp"
#include <iostream>

MonteCarloAmerican::MonteCarloAmerican(const unsigned long num_sims,
                                       const unsigned long num_intervals,
                                       const int order,
                                       Option *pOption,
                                       StatisticalDistribution *sd) 
: MonteCarloBase(num_sims, pOption, sd), _order(order), _num_intervals(num_intervals){}
MonteCarloAmerican::~MonteCarloAmerican(){}

void MonteCarloAmerican::generate_spot_paths(const double& S_0,
const double& dt,
    std::vector<std::vector<double> > &spot_paths,
std::vector<std::vector<double> > &spot_normals)
{
    for (int i = 0; i < _num_sims; i++){
        spot_paths[i][0] = S_0;
    }
    
    for (int i = 0; i < _num_sims; i++)
        {
        for (int j = 1; j < _num_intervals + 1; j++)
            {
                double S_drift = spot_paths[i][j - 1] * exp((_pOption->_r - 0.5 * (_pOption->_v) * (_pOption->_v)) * dt);
                double S_stochastic = exp(_pOption->_v * pow(dt, 0.5) * spot_normals[i][j]);
                spot_paths[i][j] = S_drift * S_stochastic;
            }
        }
}

double MonteCarloAmerican::monte_carlo_sim(const double& S_0){
    double payoff_sum = 0.0;
    double dt = _pOption->_T / static_cast<double>(_num_intervals);

    std::vector<std::vector<double> > spot_paths(_num_sims, std::vector<double>(_num_intervals + 1, 0.0));
    std::vector<std::vector<double> > spot_uniforms(_num_sims, std::vector<double>(_num_intervals + 1, 0.0));
    std::vector<std::vector<double> > spot_normals(_num_sims, std::vector<double>(_num_intervals + 1, 0.0));
    std::vector<std::vector<double> > cash_flow(_num_sims, std::vector<double>(_num_intervals + 1, 0.0));
    std::vector<std::vector<double> > value_matrix(_num_sims, std::vector<double>(_num_intervals + 1, 0.0));

    for (int i = 0; i < _num_sims; i++){
        for (int j = 0; j < _num_intervals + 1; j++){
            spot_uniforms[i][j] = rand() / static_cast<double>(RAND_MAX);
        }
    }

    for (int i = 0; i < _num_sims; i++){
        _sd->random_draws(spot_uniforms[i], spot_normals[i]);
    }

    generate_spot_paths(S_0, dt, spot_paths, spot_normals);

    for (int i = 0; i < _num_sims; i++){
        for (int j = 0; j < _num_intervals + 1; j++){
            cash_flow[i][j] = _pOption->_pPayOff->operator()(spot_paths[i][j]);
        }
        value_matrix[i][_num_intervals] = cash_flow[i][_num_intervals];
    }

    for (int j = _num_intervals - 1; j > 0; j--){
        std::vector<double> x;
        std::vector<double> y;
        std::vector<int> x_pos;
        for (int i = 0; i < _num_sims; i++)
        {
            if (cash_flow[i][j] > 0){
                x.push_back(spot_paths[i][j]);
                x_pos.push_back(i);
                y.push_back(value_matrix[i][j + 1] * exp(-_pOption->_r * dt));
            }
        }

        std::vector<double> coeff;
        polyfit(x, y, coeff, _order);

        std::vector<double> continuation(x.size(), 0.0);

        for (int k = 0; k < continuation.size(); k++){
            continuation[k] = coeff[0] + coeff[1] * pow(x[k], 1) + coeff[2] * pow(x[k], 2);

            if (cash_flow[x_pos[k]][j] > continuation[k])
            {
                value_matrix[x_pos[k]][j] = cash_flow[x_pos[k]][j];
                for (int jj = j + 1; jj < _num_intervals; jj++)
                {
                    value_matrix[x_pos[k]][jj] = 0;
                }
            }
        }

        for (int i = 0; i < _num_sims; i++){
            if (value_matrix[i][j] == 0){
                value_matrix[i][j] = value_matrix[i][j + 1] * exp(-_pOption->_r * dt);
            }
        }
    }

    for (int i = 0; i < _num_sims; i++){
        payoff_sum += value_matrix[i][1] * exp(-_pOption->_r * dt);
    }
    return payoff_sum / static_cast<double>(_num_sims);
}

#endif 