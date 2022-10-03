#ifndef __MONTE_CARLO_CPP
#define __MONTE_CARLO_CPP

#include "monte_carlo.hpp"

MonteCarloBase::MonteCarloBase() {}
MonteCarloBase::MonteCarloBase(const unsigned long num_sims,
                               Option *pOption,
                               StatisticalDistribution *sd)
    : _num_sims(num_sims), _pOption(pOption), _sd(sd) {}
MonteCarloBase::~MonteCarloBase() {}

#endif