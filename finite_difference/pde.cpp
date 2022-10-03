#ifndef __PDE_CPP
#define __PDE_CPP

#include "pde.hpp"
#include <cmath>

BlackScholesPDE::BlackScholesPDE(Option* pOption) : _pOption(pOption) {}

// Diffusion coefficient
double BlackScholesPDE::diff_coeff(double t, double x) const
{
    double vol = _pOption->_r;
    return 0.5 * vol * vol * x * x;
}

// Convection coefficient
double BlackScholesPDE::conv_coeff(double t, double x) const
{
    return (_pOption->_r) * x; // rS
}

// Zero-term coefficient
double BlackScholesPDE::zero_coeff(double t, double x) const
{
    return -_pOption->_r; // -r
}

// Source-term coefficient
double BlackScholesPDE::source_coeff(double t, double x) const
{
    return 0.0;
}

// Left-boundary condition (vanilla call option)
double BlackScholesPDE::boundary_left(double t, double x) const
{
    return 0.0; // specifically for a CALL option
}

// Right-boundary condition
double BlackScholesPDE::boundary_right(double t, double x) const
{
    // This is via the Put-Call Parity and works for a call option
    return x - (_pOption->_K) * exp(-(_pOption->_r) * ((_pOption->_T) - t));
}

// Initial condition
double BlackScholesPDE::init_cond(double x) const
{
    return _pOption->_pPayOff->operator()(x);
}

#endif 