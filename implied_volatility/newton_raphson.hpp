#ifndef __NEWTON_RAPHSON_HPP
#define __NEWTON_RAPHSON_HPP

#include <cmath>

template <typename T,
          double (T::*g)(double) const,
          double (T::*g_prime)(double) const>
double newton_raphson(double y_target,
                      double x_0,
                      double epsilon,
                      const T &root_func){

    double y = (root_func.*g)(x_0);
    double x = x_0;

    while (fabs(y_target - y) > epsilon) {
        double d_x = (root_func.*g_prime)(x);
        x += (y_target - y) / d_x;
        y = (root_func.*g)(x);
    }

    return x;
}

#endif