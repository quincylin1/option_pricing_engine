#ifndef __PAY_OFF_CPP
#define __PAY_OFF_CPP

#include "payoff.hpp"

PayOff::PayOff(){}

// ==========
// PayOffCall
// ==========

// Constructor with a single parameter
PayOffCall::PayOffCall(const double &K) { _K = K; }

// Overriden operator() method, which turns PayOffCall into a function object
double PayOffCall::operator()(const double& S) const {
    return std::max(S - _K, 0.0);
}

// =========
// PayOffPut
// =========
PayOffPut::PayOffPut(const double &K) { _K = K; }

double PayOffPut::operator()(const double& S) const {
    return std::max(_K - S, 0.0);
}

PayOffDigitalCall::PayOffDigitalCall(const double &K) { _K = K; }

double PayOffDigitalCall::operator()(const double& S) const {
    if (S >= _K) {
        return 1.0;
    }
    return 0.0;
}

PayOffDigitalPut::PayOffDigitalPut(const double &K) { _K = K; }

double PayOffDigitalPut::operator()(const double &S) const
{
    if (S <= _K)
    {
        return 1.0;
    }
    return 0.0;
}

// Constructor with two strike parameters, upper and lower barrier
PayOffDoubleDigital::PayOffDoubleDigital(const double U, const double D)
{
    _U = U;
    _D = D;
}

// Destructor (no need to use virtual keyword in source file)
PayOffDoubleDigital::~PayOffDoubleDigital() {}

// Overriden operator() method, which turns PayOffDoubleDigital into a function object
double PayOffDoubleDigital::operator()(const double S) const
{
    if (S >= _D && S <= _U)
    {
        return 1.0;
    }
    return 0.0;
}

#endif