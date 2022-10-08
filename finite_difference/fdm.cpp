#ifndef __FDM_CPP
#define __FDM_CPP

#include "fdm.hpp"

FDMBase::FDMBase(double x_max, unsigned long J,
                 double t_max, unsigned long N,
                 ConvectionDiffusionPDE *pde)
 : _x_max(x_max), _J(J), _t_max(t_max), _N(N), _pde(pde) {}

#endif 