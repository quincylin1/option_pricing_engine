#ifndef __JUMP_DIFFUSION_EUROPEAN_CPP
#define __JUMP_DIFFUSION_EUROPEAN_CPP

#include "jump_diffusion_european.hpp"

JumpDiffusionEuropean::JumpDiffusionEuropean(Option *pOption,
                                             StatisticalDistribution *sd,
                                             AnalyticEuropean *pAnalytic,
                                             const int &N, const double &m,
                                             const double &lambda, const double &nu) :
_pOption(pOption), _sd(sd), _pAnalytic(pAnalytic),
_N(N), _m(m), _lambda(lambda), _nu(nu) {}
JumpDiffusionEuropean::~JumpDiffusionEuropean(){}

double JumpDiffusionEuropean::calc_call_price(const double& S_0) const {
    double call_price = 0.0;
    double factorial = 1.0;

    double lambda_p = _lambda * _m;
    double lambda_p_T = lambda_p * _pOption->_T;

    // Calculate the finite sum over N terms
    for (int n = 0; n < _N; n++)
    {
        double sigma_n = sqrt(_pOption->_v * _pOption->_v + n * _nu * _nu / _pOption->_T);
        double r_n = _pOption->_r - _lambda * (_m - 1) + n * log(_m) / _pOption->_T;

        if (n == 0)
        {
            factorial *= 1.0;
        }
        else
        {
            factorial *= n;
        }

        call_price += (exp(-lambda_p_T) * pow(lambda_p_T, n) / factorial) *
                 _pAnalytic->calc_call_price(S_0, _pOption->_K, r_n, sigma_n, _pOption->_T);
    }
    return call_price;
}

double JumpDiffusionEuropean::calc_put_price(const double& S_0) const
{
    double put_price = 0.0;
    double factorial = 1.0;

    double lambda_p = _lambda * _m;
    double lambda_p_T = lambda_p * _pOption->_T;

    for (int n = 0; n < _N; n++)
    {
        double sigma_n = sqrt(_pOption->_v * _pOption->_v + n * _nu * _nu / _pOption->_T);
        double r_n = _pOption->_r - _lambda * (_m - 1) + n * log(_m) / _pOption->_T;

        if (n == 0)
        {
            factorial *= 1.0;
        }
        else
        {
            factorial *= n;
        }

        put_price += (exp(-lambda_p_T) * pow(lambda_p_T, n) / factorial) *
                      _pAnalytic->calc_put_price(S_0, _pOption->_K, r_n, sigma_n, _pOption->_T);
    }
    return put_price;
}

#endif 