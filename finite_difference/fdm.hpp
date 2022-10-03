#ifndef __FDM_HPP
#define __FDM_HPP

#include "pde.hpp"
#include <vector>

// Finite-difference method -- abstract base class
class FDMBase
{
protected:
    ConvectionDiffusionPDE *_pde;

    // Space discretisation
    double _x_max;                 // Spatial extent [0.0, x_max]
    unsigned long _J;              // Number of spatial differencing points
    double _dx;                    // Spatial step size (calculated from above)
    std::vector<double> _x_values; // Stores the coordinates of the x dimension

    // Time discretisation
    double _t_max;    // Temporal extent [0.0, t_max]
    unsigned long _N; // Number of temporal differencing points
    double _dt;       // Temporal step size (calculated from above)

    // Time-marching
    double _prev_t, _cur_t; // Previous and current time

    // Differencing coefficients
    double _alpha, _beta, _gamma;

    // Storage
    std::vector<double> _new_result; // New solution (becomes N + 1)
    std::vector<double> _old_result; // Old solution (becomes N)

    // Constructor
    FDMBase(double x_max, unsigned long J,
            double t_max, unsigned long N,
            ConvectionDiffusionPDE *pde);
    virtual ~FDMBase(){};

    // Override these virtual methods in derived classes for
    // specific FDM techniques , such as explicit Euler , Crank−Nicolson , etc .
    virtual void calc_step_sizes() = 0;
    virtual void set_initial_conditions() = 0;
    virtual void calc_boundary_conditions() = 0;
    virtual void calc_inner_domain() = 0;

public:
    // Carry out the actual time-stepping
    virtual void step_march() = 0;
};

class FDMEulerExplicit : public FDMBase
{
protected:
    void calc_step_sizes();
    void set_initial_conditions();
    void calc_boundary_conditions();
    void calc_inner_domain();

public:
    FDMEulerExplicit(double x_max, unsigned long J,
                     double t_max, unsigned long N,
                     ConvectionDiffusionPDE *pde);
    virtual ~FDMEulerExplicit(){};

    void step_march();
};

class FDMCrankNicolson : public FDMBase
{
protected:
    void calc_step_sizes();
    void set_initial_conditions();
    void calc_boundary_conditions();
    void calc_inner_domain();

public:
    FDMCrankNicolson(double x_max, unsigned long J,
                     double t_max, unsigned long N,
                     ConvectionDiffusionPDE *pde);
    virtual ~FDMCrankNicolson(){};

    void step_march();
};

#endif