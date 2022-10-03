#ifndef __BINOMIAL_TREE_HPP
#define __BINOMIAL_TREE_HPP

#include <vector>
#include "../options/option.hpp"


class BinomialTreeBase {
protected:
    double _u;
    double _d;
    double _pu;
    double _pd;

    Option *_pOption;

    int _num_steps;

public:
    BinomialTreeBase(){};
    BinomialTreeBase(int num_steps, Option *pOption);
    virtual ~BinomialTreeBase();

    virtual int factorial(const int &n) const = 0;

    virtual void calc_binomial_parameters(const double& S_0) = 0;

    virtual double calc_binomial_option_price(const double& S_0) = 0;

    virtual void generate_stockprice_tree(const double &S_0,
                                          std::vector<std::vector<double>> &stockprice_tree) = 0;
    virtual void generate_probability_tree(std::vector<std::vector<double>> &probability_tree) = 0;

    virtual void generate_payoff_tree(const std::vector<std::vector<double> > &stockprice_tree,
                                      std::vector<std::vector<double> > &payoff_tree) = 0;
    virtual double backwards_discounted_price(const std::vector<std::vector<double> > &probability_tree,
                                              const std::vector<std::vector<double> > &payoff_tree) = 0;
};

// class BinomialTreeJarrowRudd : public BinomialTreeBase
// {
// public:
//     BinomialTreeJarrowRudd(int _num_step, Option *_pOption);
//     virtual ~BinomialTreeJarrowRudd();

//     virtual void calc_binomial_parameters();
// };

// class BinomialTreeLeisenReimer : public BinomialTreeBase
// {
// public:
//     BinomialTreeLeisenReimer(int _num_step, Option *_pOption);
//     virtual ~BinomialTreeLeisenReimer();

//     virtual void calc_binomial_parameters();
// };

#endif