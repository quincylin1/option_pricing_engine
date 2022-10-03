#ifndef __BINOMIAL_TREE_AMERICAN_HPP
#define __BINOMIAL_TREE_AMERICAN_HPP

#include "../binomial_tree.hpp"
#include "../../analytic/european/analytic_european.hpp"
#include "../../options/option.hpp"

class BinomialTreeAmerican : public BinomialTreeBase
{
public:
    BinomialTreeAmerican(int num_steps, Option *pOption);
    virtual ~BinomialTreeAmerican();

    virtual int factorial(const int &n) const;

    virtual double calc_binomial_option_price(const double& S_0);

    virtual void generate_stockprice_tree(const double &S_0,
                                          std::vector<std::vector<double>> &stockprice_tree);
    virtual void generate_probability_tree(std::vector<std::vector<double>> &probability_tree);

    virtual void generate_payoff_tree(const std::vector<std::vector<double> > &stockprice_tree,
                                      std::vector<std::vector<double> > &payoff_tree);

    virtual double backwards_discounted_price(const std::vector<std::vector<double> > &probability_tree,
                                              const std::vector<std::vector<double> > &payoff_tree);
};

class BinomialTreeAmericanCoxRossRubinstein : public BinomialTreeAmerican
{
public:
    BinomialTreeAmericanCoxRossRubinstein(int num_steps, Option *pOption);
    virtual ~BinomialTreeAmericanCoxRossRubinstein();

    virtual void calc_binomial_parameters(const double& S_0);
};

class BinomialTreeAmericanJarrowRudd : public BinomialTreeAmerican
{
public:
    BinomialTreeAmericanJarrowRudd(int num_steps, Option *pOption);
    virtual ~BinomialTreeAmericanJarrowRudd();

    virtual void calc_binomial_parameters(const double &S_0);
};

class BinomialTreeAmericanLeisenReimer : public BinomialTreeAmerican
{
private:
    AnalyticEuropean *_pAnalytic;

public:
    BinomialTreeAmericanLeisenReimer(int num_steps,
                                     Option *pOption,
                                     AnalyticEuropean *pAnalytic);
    virtual ~BinomialTreeAmericanLeisenReimer();

    virtual void calc_binomial_parameters(const double &S_0);
    virtual double h_inv(const double &z);
};

#endif