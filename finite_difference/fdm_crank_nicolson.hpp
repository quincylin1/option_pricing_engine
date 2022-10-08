#ifndef __FDM_CRANK_NICOLSON_HPP
#define __FDM_CRANK_NICOLSON_HPP

#include "fdm.hpp"

class FDMCrankNicolson : public FDMBase
{
private:
    double r_;

protected:
    void thomas_algorithm(const std::vector<double> &a,
                          const std::vector<double> &b,
                          const std::vector<double> &c,
                          const std::vector<double> &d,
                          std::vector<double> &f);
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