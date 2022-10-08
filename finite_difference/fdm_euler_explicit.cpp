#ifndef __FDM_EULER_EXPLICIT_CPP
#define __FDM_EULER_EXPLICIT_CPP

#include "fdm_euler_explicit.hpp"
#include <iostream>
#include <fstream>

FDMEulerExplicit::FDMEulerExplicit(double x_max, unsigned long J,
                                   double t_max, unsigned long N,
                                   ConvectionDiffusionPDE *pde)
    : FDMBase(x_max, J, t_max, N, pde)
{
    calc_step_sizes();
    set_initial_conditions();
}

void FDMEulerExplicit::calc_step_sizes()
{
    _dx = _x_max / static_cast<double>(_J - 1);
    _dt = _t_max / static_cast<double>(_N - 1);
}

void FDMEulerExplicit::set_initial_conditions()
{
    double cur_spot = 0.0;

    _old_result.resize(_J, 0.0);
    _new_result.resize(_J, 0.0);
    _x_values.resize(_J, 0.0);

    for (unsigned long j = 0; j < _J; j++)
    {
        cur_spot = static_cast<double>(j) * _dx;
        // Fill old result with initial spot
        _old_result[j] = _pde->init_cond(cur_spot);
        _x_values[j] = cur_spot;
    }

    // Temporal settings
    _prev_t = 0.0;
    _cur_t = 0.0;
}

void FDMEulerExplicit::calc_boundary_conditions()
{
    // Dirichlet conditions
    _new_result[0] = _pde->boundary_left(_prev_t, _x_values[0]);
    _new_result[_J - 1] = _pde->boundary_right(_prev_t, _x_values[_J - 1]);
}

void FDMEulerExplicit::calc_inner_domain()
{
    // Differencing coefficients
    double alpha, beta, gamma;

    // Only use inner result indices (1, ... J - 2)
    for (unsigned long j = 1; j < _J - 1; j++)
    {
        double dt_sig = _dt * (_pde->diff_coeff(_prev_t, _x_values[j]));
        double dt_sig_2 = _dt * _dx * 0.5 * (_pde->conv_coeff(_prev_t, _x_values[j]));

        alpha = dt_sig - dt_sig_2;
        beta = _dx * _dx - (2.0 * dt_sig) + (_dt * _dx * _dx * _pde->zero_coeff(_prev_t, _x_values[j]));
        gamma = dt_sig + dt_sig_2;

        _new_result[j] = ((alpha * _old_result[j - 1]) +
                          (beta * _old_result[j]) +
                          (gamma * _old_result[j + 1])) /
                             (_dx * _dx) -
                         (_dt * (_pde->source_coeff(_prev_t, _x_values[j])));
    }
}

void FDMEulerExplicit::step_march()
{
    std::ofstream fdm_out("fdm.csv");

    while (_cur_t < _t_max)
    {
        _cur_t = _prev_t + _dt;
        calc_boundary_conditions(); // Reapply BC
        calc_inner_domain();        // Calculate new inner solution domain
        for (int j = 0; j < _J; j++)
        {
            fdm_out << _x_values[j] << " " << _prev_t << " "
                    << _new_result[j] << std::endl;
        }

        // Set old result to new result
        _old_result = _new_result;
        // Update previous time
        _prev_t = _cur_t;
    }

    fdm_out.close();
}

#endif 