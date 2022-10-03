#ifndef __CORRELATED_SND_HPP
#define __CORRELATED_SND_HPP

#include "standard_normal.hpp"


class CorrelatedSND : public StandardNormalDistribution{
protected:
    double rho;
    const std::vector<double>* uncorr_draws;

    // Modify an uncorrelated set of random draws to be correlated
    virtual void correlation_calc(std::vector<double> &dist_draws);

public:
    CorrelatedSND(const double _rho,
                  const std::vector<double>* _uncorr_draws);
    virtual ~CorrelatedSND();

    // Obtain a sequence of correlated random numbers from another 
    // set of random draws
    virtual void random_draws(const std::vector<double> &uniform_draws,
                              std::vector<double> &dist_draws);
};

#endif
