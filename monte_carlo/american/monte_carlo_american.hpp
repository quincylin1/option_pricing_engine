#ifndef __MONTE_CARLO_AMERICAN_HPP
#define __MONTE_CARLO_AMERICAN_HPP

#include "../monte_carlo.hpp"

class MonteCarloAmerican : public MonteCarloBase {
private:
    unsigned long _num_intervals;
    int _order;

public:
    MonteCarloAmerican(const unsigned long num_sims,
                       const unsigned long num_intervals,
                       const int order,
                       Option *pOption,
                       StatisticalDistribution *sd);
    virtual ~MonteCarloAmerican();

    virtual void generate_spot_paths(const double &S_0,
                                     const double &dt,
                                     std::vector<std::vector<double>> &spot_paths,
                                     std::vector<std::vector<double>> &spot_normals);

    virtual double monte_carlo_sim(const double& S_0);
};





#endif 