#ifndef __POLYFIT_HPP
#define __POLYFIT_HPP

#include <Eigen/Dense>
#include <Eigen/QR>

void polyfit(const std::vector<double> &x,
             const std::vector<double> &y,
             std::vector<double> &coeff,
             int order)
{
    // Create Matrix Placeholder of size n x k, n= number of datapoints, k = order of polynomial, for exame k = 3 for cubic polynomial
    Eigen::MatrixXd X(x.size(), order + 1);
    Eigen::VectorXd Y = Eigen::VectorXd::Map(&y.front(), y.size());
    Eigen::VectorXd Y_predicted;

    for (int i = 0; i < x.size(); i++)
    {
        for (int n = 0; n < order + 1; n++)
        {
            X(i, n) = pow(x.at(i), n);
        }
    }
    // std::cout << X << std::endl;

    // Solve for linear least square fit
    Y_predicted = X.householderQr().solve(Y);

    coeff.resize(order + 1);
    for (int n = 0; n < order + 1; n++)
    {
        coeff[n] = Y_predicted[n];
    }
};

#endif 