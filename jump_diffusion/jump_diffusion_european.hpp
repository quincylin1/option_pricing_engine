#ifndef __JUMP_DIFFUSTION_EUROPEAN_HPP
#define __JUMP_DIFFUSTION_EUROPEAN_HPP

#include "../analytic/european/analytic_european.hpp"
#include "../statistics/statistics.hpp"
#include "../options/option.hpp"


class JumpDiffusionEuropean {
private:
    Option *_pOption;
    StatisticalDistribution *_sd;
    AnalyticEuropean *_pAnalytic;

    int _N;       
    double _m;
    double _lambda;
    double _nu;

public:
    JumpDiffusionEuropean(Option *pOption,
                          StatisticalDistribution *sd,
                          AnalyticEuropean *pAnalytic,
                          const int &N, const double &m,
                          const double &lambda, const double &nu);
    virtual ~JumpDiffusionEuropean();

    virtual double calc_call_price(const double& S_0) const;
    virtual double calc_put_price(const double& S_0) const;
};

#endif 