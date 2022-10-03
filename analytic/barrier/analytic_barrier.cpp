#ifndef __ANALYTIC_BARRIER_CPP
#define __ANALYTIC_BARRIER_CPP

#include "analytic_barrier.hpp"
#include <cmath>

AnalyticBarrier::AnalyticBarrier(StatisticalDistribution *sd, AnalyticEuropean* pAnalytic)
 : _sd(sd), _pAnalytic(pAnalytic) {}

double AnalyticBarrier::calc_gamma(const double &r, const double &q,
                                   const double &v) const
{
    return (r - q + 0.5 * v * v) / (v * v);
}

double AnalyticBarrier::calc_eta(const double &S, const double &B,
                                 const double &v, const double &tau,
                                 const double &K, const double &gamma) const
{
    return log((B * B) / (S * K)) / (v * sqrt(tau)) + gamma * v * sqrt(tau);
}

double AnalyticBarrier::calc_mu(const double &S, const double &B,
                                const double &v, const double &tau,
                                const double &gamma) const
{
    return log(S / B) / (v * sqrt(tau)) + gamma * v * sqrt(tau);
}

double AnalyticBarrier::calc_lambda(const double &S, const double &B,
                                    const double &v, const double &tau,
                                    const double &gamma) const
{
    return log(B / S) / (v * sqrt(tau)) + gamma * v * sqrt(tau);
}

double AnalyticBarrier::calc_up_and_out_call_price(const double &S, const double &K, const double &r,
                                                   const double &v, const double &T, const double &q,
                                                   const double &B, const double &tau) const
{

    if (B < K)
    {
        return 0.0;
    }

    double european_call_price = _pAnalytic->calc_call_price(S, K, r, v, T);
    double up_and_in_call_price = calc_up_and_in_call_price(S, K, r, v, T, q, B, tau);
    return european_call_price - up_and_in_call_price;
}

double AnalyticBarrier::calc_up_and_in_call_price(const double &S, const double &K, const double &r,
                                                      const double &v, const double &T, const double &q,
                                                      const double &B, const double &tau) const
{
    if (B < K) {
        return _pAnalytic->calc_call_price(S, K, r, v, T);
    }

    double gamma = calc_gamma(r, q, v);
    double mu = calc_mu(S, B, v, tau, gamma);
    double lambda = calc_lambda(S, B, v, tau, gamma);
    double eta = calc_eta(S, B, v, tau, K, gamma);

    return S * exp(-q * tau) * _sd->cdf(mu) 
      - K * exp(-r * tau) * _sd->cdf(mu - v * sqrt(tau)) 
      - S * exp(-q * tau) * pow((B / S), 2.0 * gamma) * (_sd->cdf(-eta) - _sd->cdf(-lambda)) 
      + K * exp(-r * tau) * pow((B / S), 2.0 * gamma - 2.0) * (_sd->cdf(-eta + v * sqrt(tau)) - _sd->cdf(-lambda + v * sqrt(tau)));
}

double AnalyticBarrier::calc_down_and_out_call_price(const double &S, const double &K, const double &r,
                                                   const double &v, const double &T, const double &q,
                                                   const double &B, const double &tau) const
{

    if (B < K)
    {
        double european_call_price = _pAnalytic->calc_call_price(S, K, r, v, T);
        double down_and_in_call_price = calc_down_and_in_call_price(S, K, r, v, T, q, B, tau);
        return european_call_price - down_and_in_call_price;
    }

    double gamma = calc_gamma(r, q, v);
    double mu = calc_mu(S, B, v, tau, gamma);
    double lambda = calc_lambda(S, B, v, tau, gamma);

    return S * exp(-q * tau) * _sd->cdf(mu) 
      - K * exp(-r * tau) * _sd->cdf(mu - v * sqrt(tau)) 
      - S * exp(-q * tau) * pow((B / S), 2.0 * gamma) * _sd->cdf(lambda) 
      + K * exp(-r * tau) * pow((B / S), 2.0 * gamma - 2.0) * _sd->cdf(lambda - v * sqrt(tau));
}

double AnalyticBarrier::calc_down_and_in_call_price(const double &S, const double &K, const double &r,
                                                   const double &v, const double &T, const double &q,
                                                   const double &B, const double &tau) const
{
    if (B < K){
        double gamma = calc_gamma(r, q, v);
        double eta = calc_eta(S, B, v, tau, K, gamma);

        return S * exp(-q * tau) * pow((B / S), 2.0 * gamma) * _sd->cdf(eta) 
          - K * exp(-r * tau) * pow((B / S), 2.0 * gamma - 2.0) * _sd->cdf(eta - v * sqrt(tau));
    }

    double european_call_price = _pAnalytic->calc_call_price(S, K, r, v, T);
    double down_and_out_call_price = calc_down_and_out_call_price(S, K, r, v, T, q, B, tau);
    return european_call_price - down_and_out_call_price;
}

double AnalyticBarrier::calc_up_and_out_put_price(const double &S, const double &K, const double &r,
                                                   const double &v, const double &T, const double &q,
                                                   const double &B, const double &tau) const
{

    if (B < K)
    {
        double gamma = calc_gamma(r, q, v);
        double mu = calc_mu(S, B, v, tau, gamma);
        double lambda = calc_lambda(S, B, v, tau, gamma);

        return -S * exp(-q * tau) * _sd->cdf(-mu) 
          + K * exp(-r * tau) * _sd->cdf(-mu + v * sqrt(tau)) 
          + S * exp(-q * tau) * pow((B / S), 2.0 * gamma) * _sd->cdf(-lambda) 
          - K * exp(-r * tau) * pow((B / S), 2.0 * gamma - 2.0) * _sd->cdf(-lambda + v * sqrt(tau));
    }

    double european_put_price = _pAnalytic->calc_put_price(S, K, r, v, T);
    double up_and_in_put_price = calc_up_and_in_put_price(S, K, r, v, T, q, B, tau);
    return european_put_price - up_and_in_put_price;
}

double AnalyticBarrier::calc_up_and_in_put_price(const double &S, const double &K, const double &r,
                                                  const double &v, const double &T, const double &q,
                                                  const double &B, const double &tau) const
{
    if (B < K)
    {
        double european_put_price = _pAnalytic->calc_put_price(S, K, r, v, T);
        double up_and_out_put_price = calc_up_and_out_put_price(S, K, r, v, T, q, B, tau);
        return european_put_price - up_and_out_put_price;
    }

    double gamma = calc_gamma(r, q, v);
    double eta = calc_eta(S, B, v, tau, K, gamma);

    return -S * exp(-q * tau) * pow((B / S), 2.0 * gamma) * _sd->cdf(-eta) 
      + K * exp(-r * tau) * pow((B / S), 2.0 * gamma - 2.0) * _sd->cdf(-eta + v * sqrt(tau));
}

double AnalyticBarrier::calc_down_and_out_put_price(const double &S, const double &K, const double &r,
                                                     const double &v, const double &T, const double &q,
                                                     const double &B, const double &tau) const
{

    if (B < K)
    {
        double european_put_price = _pAnalytic->calc_put_price(S, K, r, v, T);
        double down_and_in_put_price = calc_down_and_in_put_price(S, K, r, v, T, q, B, tau);
        return european_put_price - down_and_in_put_price;
    }

    return 0.0;
}

double AnalyticBarrier::calc_down_and_in_put_price(const double &S, const double &K, const double &r,
                                                    const double &v, const double &T, const double &q,
                                                    const double &B, const double &tau) const
{
    if (B < K)
    {
        double gamma = calc_gamma(r, q, v);
        double mu = calc_mu(S, B, v, tau, gamma);
        double eta = calc_eta(S, B, v, tau, K, gamma);
        double lambda = calc_lambda(S, B, v, tau, gamma);

        return -S * exp(-q * tau) * _sd->cdf(-mu) 
          + K * exp(-r * tau) * _sd->cdf(-mu + v * sqrt(tau)) 
          + S * exp(-q * tau) * pow((B / S), 2.0 * gamma) * (_sd->cdf(eta) - _sd->cdf(lambda)) 
          - K * exp(-r * tau) * pow((B / S), 2.0 * gamma - 2.0) * (_sd->cdf(eta - v * sqrt(tau) - _sd->cdf(lambda - v * sqrt(lambda))));
    }

    return _pAnalytic->calc_put_price(S, K, r, v, T);
}
#endif