#ifndef __PDE_HPP
#define __PDE_HPP

#include "../options/option.hpp"


class ConvectionDiffusionPDE {
public:
    // PDE Coefficients
    virtual double diff_coeff(double t, double x) const = 0;
    virtual double conv_coeff(double t, double x) const = 0;
    virtual double zero_coeff(double t, double x) const = 0;
    virtual double source_coeff(double t, double x) const = 0;

    // Boundary and initial conditions
    virtual double boundary_left(double t, double x) const = 0;
    virtual double boundary_right(double t, double x) const = 0;

    virtual double init_cond(double x) const = 0;
};

class BlackScholesPDE : public ConvectionDiffusionPDE{
public:
    Option *_pOption;
    BlackScholesPDE(Option *pOption);
    virtual ~BlackScholesPDE(){};

    virtual double diff_coeff(double t, double x) const;
    virtual double conv_coeff(double t, double x) const;
    virtual double zero_coeff(double t, double x) const;
    virtual double source_coeff(double t, double x) const;

    virtual double boundary_left(double t, double x) const;
    virtual double boundary_right(double t, double x) const;

    virtual double init_cond(double x) const;
};

#endif 