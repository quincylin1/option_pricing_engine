#ifndef __PAY_OFF_HPP
#define __PAY_OFF_HPP

#include <algorithm>


class PayOff{
  public:
      PayOff();
      virtual ~PayOff() {};

      virtual double operator()(const double &S) const = 0;
};

class PayOffCall : public PayOff{
private:
    double _K; // Strike price

public:
    PayOffCall(const double &K);
    virtual ~PayOffCall(){};

    // Override PayOff virtual operator() function 
    virtual double operator()(const double &S) const;
};

class PayOffPut : public PayOff{
private:
    double _K;

public:
    PayOffPut(const double &K);
    virtual ~PayOffPut(){};

    virtual double operator()(const double &S) const;
};

class PayOffDigitalCall : public PayOff {
private:
    double _K;

public:
    PayOffDigitalCall(const double &K);
    virtual ~PayOffDigitalCall(){};

    virtual double operator()(const double &S) const;
};

class PayOffDigitalPut : public PayOff
{
private:
    double _K;

public:
    PayOffDigitalPut(const double &K);
    virtual ~PayOffDigitalPut(){};

    virtual double operator()(const double &S) const;
};

class PayOffDoubleDigital : public PayOff
{
private:
    double _U; // Upper strike price
    double _D; // Lower strike price

public:
    // Two strike parameters
    PayOffDoubleDigital(const double U, const double D);
    virtual ~PayOffDoubleDigital();

    // Pay-off is one if spot within strike barrier, 0 otherwise
    virtual double operator()(const double S) const;
};
#endif
