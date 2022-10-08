#ifndef __FDM_EULER_EXPLICIT_HPP
#define __FDM_EULER_EXPLICIT_HPP

#include "fdm.hpp"

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

#endif 