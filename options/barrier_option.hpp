#ifndef __BARRIER_OPTION_HPP
#define __BARRIER_OPTION_HPP

#include "option.hpp"

class BarrierOption : public Option {
protected:
    double _B;

public:
    BarrierOption(double K,
                  double r, double T,
                  double v, double q,
                  double B,
                  PayOff *pPayOff);
    virtual ~BarrierOption();

    virtual double calc_payoff(const std::vector<double> &spot_prices) const = 0;
};

class BarrierOptionUpAndOut : public BarrierOption {
public:
    BarrierOptionUpAndOut(double K,
                          double r, double T,
                          double v, double q,
                          double B,
                          PayOff *pPayOff);
    virtual ~BarrierOptionUpAndOut();

    virtual double calc_payoff(const std::vector<double> &spot_prices) const;
};

class BarrierOptionDownAndOut : public BarrierOption
{
public:
    BarrierOptionDownAndOut(double K,
                            double r, double T,
                            double v, double q,
                            double B,
                            PayOff *pPayOff);
    virtual ~BarrierOptionDownAndOut();

    virtual double calc_payoff(const std::vector<double> &spot_prices) const;
};

class BarrierOptionUpAndIn : public BarrierOption
{
public:
    BarrierOptionUpAndIn(double K,
                         double r, double T,
                         double v, double q,
                         double B,
                         PayOff *pPayOff);
    virtual ~BarrierOptionUpAndIn();

    virtual double calc_payoff(const std::vector<double> &spot_prices) const;
};

class BarrierOptionDownAndIn : public BarrierOption
{
public:
    BarrierOptionDownAndIn(double K,
                           double r, double T,
                           double v, double q,
                           double B,
                           PayOff *pPayOff);
    virtual ~BarrierOptionDownAndIn();

    virtual double calc_payoff(const std::vector<double> &spot_prices) const;
};

#endif 