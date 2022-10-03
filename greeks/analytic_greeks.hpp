#ifndef __ANALYTIC_GREEKS_H
#define __ANALYTIC_GREEKS_H

#include "../options/option.hpp"
#include "../statistics/statistics.hpp"
#include "../analytic/european/analytic_european.hpp"


class AnalyticGreeks{
private:
    Option* _pOption;
    AnalyticEuropean *_pAnalytic;
    StatisticalDistribution *_sd;

public:
    AnalyticGreeks();
    AnalyticGreeks(Option *pOption,
                   AnalyticEuropean *pAnalytic,
                   StatisticalDistribution *sd);

    virtual ~AnalyticGreeks();

    double calc_call_delta(const double &S_0) const;
    double calc_call_gamma(const double &S_0) const;
    double calc_call_vega(const double &S_0) const;
    double calc_call_theta(const double &S_0) const;
    double calc_call_rho(const double &S_0) const;

    double calc_put_delta(const double &S_0) const;
    double calc_put_gamma(const double &S_0) const;
    double calc_put_vega(const double &S_0) const;
    double calc_put_theta(const double &S_0) const;
    double calc_put_rho(const double &S_0) const;
};

#endif