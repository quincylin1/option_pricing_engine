#ifndef __BINOMIAL_TREE_EUROPEAN_CPP
#define __BINOMIAL_TREE_EUROPEAN_CPP

#include "binomial_tree_european.hpp"
#include <vector>
#include <cmath>
#include <iostream>

BinomialTreeEuropean::BinomialTreeEuropean(int num_steps, Option *pOption) 
: BinomialTreeBase(num_steps, pOption){}
BinomialTreeEuropean::~BinomialTreeEuropean(){}

int BinomialTreeEuropean::factorial(const int& n) const {
    int factorial = 1;
    if (n == 0)
        return 1;
    for (int i = 1; i <= n; i++)
    {
        factorial *= i;
    }
    return factorial;
}

double BinomialTreeEuropean::calc_binomial_option_price(const double& S_0)
{
    std::vector<std::vector<double> > stockprice_tree(_num_steps + 1, std::vector<double>(_num_steps + 1, 0.0));
    std::vector<std::vector<double> > probability_tree(_num_steps + 1, std::vector<double>(_num_steps + 1, 0.0));
    std::vector<std::vector<double> > payoff_tree(_num_steps + 1, std::vector<double>(_num_steps + 1, 0.0));

    calc_binomial_parameters(S_0);
    generate_stockprice_tree(S_0, stockprice_tree);

    generate_probability_tree(probability_tree);
    generate_payoff_tree(stockprice_tree, payoff_tree);

    return backwards_discounted_price(probability_tree, payoff_tree);
}

void BinomialTreeEuropean::generate_stockprice_tree(const double& S_0, std::vector<std::vector<double> > &stockprice_tree)
{
    for (int j = 0; j < stockprice_tree[0].size(); j++)
    {
        for (int i = 0; i < j + 1; i++)
        {
            stockprice_tree[i][j] = S_0 * pow(_u, i) * pow(_d, j - i);
        }
    }
}

void BinomialTreeEuropean::generate_probability_tree(std::vector<std::vector<double> > &probability_tree)
{
    for (int j = 0; j < probability_tree[0].size(); j++)
    {
        for (int i = 0; i < j + 1; i++)
        {
            probability_tree[i][j] = factorial(j) / (factorial(i) * factorial(j - i)) *
                                     pow(_pu, i) * pow(_pd, j - i);
        }
    }
}

void BinomialTreeEuropean::generate_payoff_tree(const std::vector<std::vector<double> > &stockprice_tree,
                                                std::vector<std::vector<double> > &payoff_tree)
{
    for (int i = 0; i < _num_steps + 1; i++)
    {
        payoff_tree[i][_num_steps] = _pOption->_pPayOff->operator()(stockprice_tree[i][_num_steps]);
    }
}

double BinomialTreeEuropean::backwards_discounted_price(const std::vector<std::vector<double> > &probability_tree,
                                                        const std::vector<std::vector<double> > &payoff_tree)
{
    double payoff = 0.0;
    for (int i = 0; i < _num_steps + 1; i++)
    {
        payoff += payoff_tree[i][_num_steps] * probability_tree[i][_num_steps];
    }
    return payoff * exp(-_pOption->_r * _pOption->_T);
}

BinomialTreeEuropeanCoxRossRubinstein::BinomialTreeEuropeanCoxRossRubinstein(int num_steps, Option *pOption)
: BinomialTreeEuropean(num_steps, pOption){}
BinomialTreeEuropeanCoxRossRubinstein::~BinomialTreeEuropeanCoxRossRubinstein(){}

void BinomialTreeEuropeanCoxRossRubinstein::calc_binomial_parameters(const double &S_0)
{
    double dt = _pOption->_T / static_cast<double>(_num_steps);
    _u = exp(_pOption->_v * sqrt(dt));
    _d = 1.0 / _u;
    _pu = (exp((_pOption->_r - _pOption->_q) * dt) - _d) / (_u - _d);
    _pd = 1.0 - _pu;
}

BinomialTreeEuropeanJarrowRudd::BinomialTreeEuropeanJarrowRudd(int num_steps, Option *pOption)
    : BinomialTreeEuropean(num_steps, pOption) {}
BinomialTreeEuropeanJarrowRudd::~BinomialTreeEuropeanJarrowRudd() {}

void BinomialTreeEuropeanJarrowRudd::calc_binomial_parameters(const double &S_0)
{
    double dt = _pOption->_T / static_cast<double>(_num_steps);
    _u = exp((_pOption->_r - _pOption->_q - 0.5 * pow(_pOption->_v, 2)) * dt + _pOption->_v * sqrt(dt));
    _d = exp((_pOption->_r - _pOption->_q - 0.5 * pow(_pOption->_v, 2)) * dt - _pOption->_v * sqrt(dt));
    _pu = 0.5;
    _pd = 1.0 - _pu;
}

BinomialTreeEuropeanLeisenReimer::BinomialTreeEuropeanLeisenReimer(int num_steps, 
Option *pOption, AnalyticEuropean *pAnalytic)
: BinomialTreeEuropean(num_steps, pOption), _pAnalytic(pAnalytic) {}
BinomialTreeEuropeanLeisenReimer::~BinomialTreeEuropeanLeisenReimer(){}

double BinomialTreeEuropeanLeisenReimer::h_inv(const double& z) {
    int sign_z = 1;
    if (z < 0)
        sign_z = -1;
    double h_exp = pow(z / (_num_steps + 1 / 3 + (0.1 / (_num_steps + 1))), 2) * (_num_steps + 1 / 6);

    return 0.5 + 0.5 * sign_z * sqrt(1 - exp(-h_exp));
}
void BinomialTreeEuropeanLeisenReimer::calc_binomial_parameters(const double& S_0)
{
    double dt = _pOption->_T / static_cast<double>(_num_steps);
    double d_1 = _pAnalytic->d_j(1, S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T);
    double d_2 = _pAnalytic->d_j(2, S_0, _pOption->_K, _pOption->_r, _pOption->_v, _pOption->_T);

    _pu = h_inv(d_2);
    _pd = 1.0 - _pu;

    double pu_prime = h_inv(d_1);

    _u = exp((_pOption->_r - _pOption->_q) * dt) * pu_prime / _pu;
    _d = exp((_pOption->_r - _pOption->_q) * dt) * (1.0 - pu_prime) / (1.0 - _pu);
}

#endif 