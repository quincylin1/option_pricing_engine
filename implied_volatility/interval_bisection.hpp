#ifndef __INTERVAL_BISECTION_HPP
#define __INTERVAL_BISECTION_HPP

#include <cmath>

template <typename T>
double interval_bisection(double y_target,
                          double left,
                          double right,
                          double epsilon,
                          T g) {

    double x = 0.5 * (left + right);
    double y = g(x);

    do {
        if (y < y_target) {
            left = x;
        }
        if (y > y_target) {
            right = x;
        }

        x = 0.5 * (left + right);
        y = g(x);
    } while (fabs(y - y_target) > epsilon);

    return x;
}

#endif 