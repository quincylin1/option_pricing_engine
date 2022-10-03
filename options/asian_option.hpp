#ifndef __ASIAN_OPTION_HPP
#define __ASIAN_OPTION_HPP

#include "../payoffs/payoff.h"
#include "option.hpp"
#include <vector>

class AsianOption : public Option
{
public:
    AsianOption(double K,
                double r, double T,
                double v, double q,
                PayOff *pPayOff);
    virtual ~AsianOption(){};

    // Pure virtual payoff operator (this will determine arithmetic
    // or geometric)
    virtual double pay_off_price(const std::vector<double> &spot_prices) const = 0;
};

class AsianOptionArithmetic : public AsianOption
{
public:
    AsianOptionArithmetic(double K,
                          double r, double T,
                          double v, double q,
                          PayOff *pPayOff);
    virtual ~AsianOptionArithmetic(){};

    virtual double pay_off_price(const std::vector<double> &spot_prices) const;
};

class AsianOptionGeometric : public AsianOption
{
public:
    AsianOptionGeometric(double K,
                         double r, double T,
                         double v, double q,
                         PayOff *pPayOff);
    virtual ~AsianOptionGeometric(){};

    virtual double pay_off_price(const std::vector<double> &spot_prices) const;
};

#endif