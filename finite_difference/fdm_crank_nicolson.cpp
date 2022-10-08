#ifndef __FDM_CRANK_NICOLSON_CPP
#define __FDM_CRANK_NICOLSON_CPP

#include "fdm_crank_nicolson.hpp"
#include <iostream>
#include <fstream>

FDMCrankNicolson::FDMCrankNicolson(double x_max, unsigned long J,
                                   double t_max, unsigned long N,
                                   ConvectionDiffusionPDE *pde)
    : FDMBase(x_max, J, t_max, N, pde)
{
    calc_step_sizes();
    set_initial_conditions();
}

void FDMCrankNicolson::calc_step_sizes()
{
    _dx = _x_max / static_cast<double>(_J - 1);
    _dt = _t_max / static_cast<double>(_N - 1);
    r_ = _dt / (_dx * _dx);
}

void FDMCrankNicolson::set_initial_conditions()
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

void FDMCrankNicolson::calc_boundary_conditions()
{
    // Dirichlet conditions
    _new_result[0] = _pde->boundary_left(_prev_t, _x_values[0]);
    _new_result[_J - 1] = _pde->boundary_right(_prev_t, _x_values[_J - 1]);
}

void FDMCrankNicolson::calc_inner_domain()
{
    // Only use inner result indices (1, ... J - 2)
    std::vector<double> a(_J - 1, -r_ / 2);
    std::vector<double> b(_J, 1 + r_);
    std::vector<double> c(_J, -r_ / 2);
    std::vector<double> d(_J, 0.0);

    for (int i = 0; i < _J; i++)
    {
        d[i] = 0.5 * r_ * _old_result[i + 1] + (1 - r_) * _old_result[i] + 0.5 * r_ * _old_result[i - 1];
    }

    thomas_algorithm(a, b, c, d, _new_result);
}

void FDMCrankNicolson::thomas_algorithm(const std::vector<double> &a,
                                        const std::vector<double> &b,
                                        const std::vector<double> &c,
                                        const std::vector<double> &d,
                                        std::vector<double> &f)
{
    std::vector<double> c_star(_J - 1, 0.0);
    std::vector<double> d_star(_J - 1, 0.0);

    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for (int i = 1; i < _J - 1; i++)
    {
        double m = 1.0 / (b[i] - a[i] * c_star[i - 1]);
        c_star[i] = c[i] * m;
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) * m;
    }

    for (int i = _J - 2; i > 0; i--)
    {
        f[i] = d_star[i] - c_star[i] * d[i + 1];
    }
}

void FDMCrankNicolson::step_march()
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