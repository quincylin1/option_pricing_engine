#ifndef __FDM_GREEKS_HPP
#define __FDM_GREEKS_HPP

#include "../options/option.hpp"
#include "../statistics/statistics.hpp"
#include "../analytic/european/analytic_european.hpp"

class FDMGreeks
{
private:
    Option *_pOption;
    AnalyticEuropean *_pAnalytic;
    StatisticalDistribution *_sd;

public:
    FDMGreeks();
    FDMGreeks(Option * _pOption,
              AnalyticEuropean *_pAnalytic,
              StatisticalDistribution *_sd);

    virtual ~FDMGreeks();

    double calc_call_delta(const double &S_0, const double &delta_S) const;
    double calc_call_gamma(const double &S_0, const double &delta_S) const;
    double calc_call_vega(const double &S_0, const double &delta_v) const;
    double calc_call_theta(const double &S_0, const double &delta_T) const;
    double calc_call_rho(const double &S_0, const double &delta_r) const;

    double calc_put_delta(const double &S_0, const double &delta_S) const;
    double calc_put_gamma(const double &S_0, const double &delta_S) const;
    double calc_put_vega(const double &S_0, const double &delta_v) const;
    double calc_put_theta(const double &S_0, const double &delta_T) const;
    double calc_put_rho(const double &S_0, const double &delta_r) const;
};

#endif