#ifndef __ANALYTIC_BARRIER_HPP
#define __ANALYTIC_BARRIER_HPP

#include "../../statistics/statistics.hpp"
#include "../european/analytic_european.hpp"

class AnalyticBarrier
{
private:
    StatisticalDistribution *_sd;
    AnalyticEuropean *_pAnalytic;

public:
    AnalyticBarrier(){};
    AnalyticBarrier(StatisticalDistribution *sd, AnalyticEuropean *pAnalytic);
    ~AnalyticBarrier(){};

    virtual double calc_gamma(const double &r, const double &q,
                              const double &v) const;

    virtual double calc_eta(const double &S, const double &B,
                            const double &v, const double &tau,
                            const double &K, const double &gamma) const;

    virtual double calc_nu(const double &S, const double &B,
                           const double &v, const double &tau,
                           const double &gamma) const;

    virtual double calc_lambda(const double &S, const double &B,
                               const double &v, const double &tau,
                               const double &gamma) const;

    virtual double calc_up_and_out_call_price(const double &S, const double &K, const double &r,
                                              const double &v, const double &T, const double &q,
                                              const double &B, const double &tau) const;
    virtual double calc_up_and_in_call_price(const double &S, const double &K, const double &r,
                                             const double &v, const double &T, const double &q,
                                             const double &B, const double &tau) const;

    virtual double calc_down_and_out_call_price(const double &S, const double &K, const double &r,
                                              const double &v, const double &T, const double &q,
                                              const double &B, const double &tau) const;
    virtual double calc_down_and_in_call_price(const double &S, const double &K, const double &r,
                                             const double &v, const double &T, const double &q,
                                             const double &B, const double &tau) const;

    virtual double calc_up_and_out_put_price(const double &S, const double &K, const double &r,
                                             const double &v, const double &T, const double &q,
                                             const double &B, const double &tau) const;
    virtual double calc_up_and_in_put_price(const double &S, const double &K, const double &r,
                                             const double &v, const double &T, const double &q,
                                             const double &B, const double &tau) const;

    virtual double calc_down_and_out_put_price(const double &S, const double &K, const double &r,
                                             const double &v, const double &T, const double &q,
                                             const double &B, const double &tau) const;
    virtual double calc_down_and_in_put_price(const double &S, const double &K, const double &r,
                                             const double &v, const double &T, const double &q,
                                             const double &B, const double &tau) const;
};

#endif