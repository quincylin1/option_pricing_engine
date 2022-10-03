#ifndef __BARRIER_OPTION_CPP
#define __BARRIER_OPTION_CPP

#include "barrier_option.hpp"
#include <cmath>
#include <vector>

BarrierOption::BarrierOption(double K,
                             double r, double T,
                             double v, double q,
                             double B,
                             PayOff *pPayOff)
    : Option(K, r, T, v, q, pPayOff), _B(B) {}
BarrierOption::~BarrierOption(){}

BarrierOptionUpAndOut::BarrierOptionUpAndOut(double _S_0, double _K,
                                             double _r, double _T,
                                             double _v, double _q,
                                             double _B,
                                             PayOff *_pPayOff)
: BarrierOption(_S_0, _K, _r, _T, _v, _q, _B, _pPayOff) {}
BarrierOptionUpAndOut::~BarrierOptionUpAndOut(){}

double BarrierOptionUpAndOut::calc_payoff(const std::vector<double> &spot_prices) const {
    double M_max = *max_element(spot_prices.begin(), spot_prices.end());
    return pay_off->operator()(spot_prices[spot_prices.size() - 1]) * (M_max <= B);
}

BarrierOptionDownAndOut::BarrierOptionDownAndOut(double _S_0, double _K,
                                             double _r, double _T,
                                             double _v, double _q,
                                             double _B,
                                             PayOff *_pay_off)
: BarrierOption(_S_0, _K, _r, _T, _v, _q, _B, _pay_off) {}
BarrierOptionDownAndOut::~BarrierOptionDownAndOut() {}

double BarrierOptionDownAndOut::calc_payoff(const std::vector<double> &spot_prices) const
{
    double M_min = *min_element(spot_prices.begin(), spot_prices.end());
    return pay_off->operator()(spot_prices[spot_prices.size() - 1]) * (M_min > B);
}

BarrierOptionUpAndIn::BarrierOptionUpAndIn(double _S_0, double _K,
                                           double _r, double _T,
                                           double _v, double _q,
                                           double _B,
                                           PayOff *_pay_off)
    : BarrierOption(_S_0, _K, _r, _T, _v, _q, _B, _pay_off) {}
BarrierOptionUpAndIn::~BarrierOptionUpAndIn() {}

double BarrierOptionUpAndIn::calc_payoff(const std::vector<double> &spot_prices) const
{
    double M_max = *min_element(spot_prices.begin(), spot_prices.end());
    return pay_off->operator()(spot_prices[spot_prices.size() - 1]) * (M_max > B);
}

BarrierOptionDownAndIn::BarrierOptionDownAndIn(double _S_0, double _K,
                                             double _r, double _T,
                                             double _v, double _q,
                                             double _B,
                                             PayOff *_pay_off)
    : BarrierOption(_S_0, _K, _r, _T, _v, _q, _B, _pay_off) {}
BarrierOptionDownAndIn::~BarrierOptionDownAndIn() {}

double BarrierOptionDownAndIn::calc_payoff(const std::vector<double> &spot_prices) const
{
    double M_min = *min_element(spot_prices.begin(), spot_prices.end());
    return pay_off->operator()(spot_prices[spot_prices.size() - 1]) * (M_min <= B);
}

#endif 