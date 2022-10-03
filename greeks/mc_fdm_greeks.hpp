#ifndef __MC_FDM_H
#define __MC_FDM_H

#include "../options/option.hpp"
#include "../payoffs/payoff.hpp"
#include "../statistics/statistics.hpp"

class MonteCarloFDMGreeks {
private:

    unsigned long _num_sims;

    Option *_pOption;
    PayOff *_pPayOff;
    StatisticalDistribution *_sd;

public:
    MonteCarloFDMGreeks();
    MonteCarloFDMGreeks(unsigned long num_sims,
                        Option *pOption, PayOff *pPayOff,
                        StatisticalDistribution *sd);

    virtual ~MonteCarloFDMGreeks();

    double calc_S_drift(const double S_0, const double r,
                        const double v, const double T) const;

    double calc_S_stochastic(const double v, const double T, const double spot_normal) const;

    void monte_carlo_sim(const double S_0, const double r,
                         const double v, const double T,
                         const double K,
                         const double delta_S,
                         double &price_Sp,
                         double &price_S,
                         double &price_Sm);

    double calc_delta(const double S_0, const double delta_S);
    double calc_gamma(const double S_0, const double delta_S);
};

#endif 