#ifndef __BINOMIAL_TREE_EUROPEAN_HPP
#define __BINOMIAL_TREE_EUROPEAN_HPP

#include "../../options/option.hpp"
#include "../binomial_tree.hpp"
#include "../../analytic/european/analytic_european.hpp"

class BinomialTreeEuropean : public BinomialTreeBase {
public:
    BinomialTreeEuropean(int num_steps, Option *pOption);
    virtual ~BinomialTreeEuropean();

    virtual int factorial(const int &n) const;

    virtual double calc_binomial_option_price(const double& S_0);

    virtual void generate_stockprice_tree(
        const double& S_0,
        std::vector<std::vector<double> > &stockprice_tree);
    virtual void generate_probability_tree(std::vector<std::vector<double> > &probability_tree);

    virtual void generate_payoff_tree(const std::vector<std::vector<double> > &stockprice_tree,
                                      std::vector<std::vector<double> > &payoff_tree);

    virtual double backwards_discounted_price(const std::vector<std::vector<double> > &probability_tree,
                                              const std::vector<std::vector<double> > &payoff_tree);
};

class BinomialTreeEuropeanCoxRossRubinstein : public BinomialTreeEuropean
{
public:
    BinomialTreeEuropeanCoxRossRubinstein(int num_steps, Option *pOption);
    virtual ~BinomialTreeEuropeanCoxRossRubinstein();

    virtual void calc_binomial_parameters(const double &S_0);
};

class BinomialTreeEuropeanJarrowRudd : public BinomialTreeEuropean
{
public:
    BinomialTreeEuropeanJarrowRudd(int num_steps, Option *pOption);
    virtual ~BinomialTreeEuropeanJarrowRudd();

    virtual void calc_binomial_parameters(const double &S_0);
};

class BinomialTreeEuropeanLeisenReimer : public BinomialTreeEuropean 
{
private:
    AnalyticEuropean *_pAnalytic;

public:
    BinomialTreeEuropeanLeisenReimer(int num_steps,
                                     Option *pOption,
                                     AnalyticEuropean *pAnalytic);
    virtual ~BinomialTreeEuropeanLeisenReimer();

    virtual void calc_binomial_parameters(const double& S_0);
    virtual double h_inv(const double &z);
};


#endif