#ifndef __STANDARD_GAMMA_CPP
#define __STANDARD_GAMMA_CPP

#include "standard_gamma.hpp"
#include <iostream>

std::vector<double> p_values = {676.5203681218851,
                         -1259.1392167224028,
                         771.32342877765313,
                         -176.61502916214059,
                         12.507343278686905,
                         -0.13857109526572012,
                         9.9843695780195716e-6,
                         1.5056327351493116e-7};

StandardGammaDistribution::StandardGammaDistribution(double alpha) : alpha_(alpha) {}
StandardGammaDistribution::~StandardGammaDistribution() {}

double StandardGammaDistribution::gamma_function(double z)
{
    double y;
    if (z < 0.5)
    {
        y = M_PI / (sin(M_PI * z) * gamma_function(1.0 - z));
    }
    else{
        z -= 1.0;
        double x = 0.99999999999980993;

        for (int i = 0; i < p_values.size(); i++) {
            x += p_values[i] / (z + static_cast<double>(i) + 1.0);
        }
        double t = z + static_cast<double>(p_values.size()) - 0.5;
        y = sqrt(2.0 * M_PI) * pow(t, (z + 0.5)) * exp(-t) * x;
    }
    return y;
}

double StandardGammaDistribution::pdf(const double &x)
{
    if (x < 0) {
        std::cout << "x must be greater than or equal to 0" << std::endl;
    }

    double z = x;
    double gamma = gamma_function(z);
    return pow(x, alpha_ - 1.0) * exp(-x) / gamma;
}

double StandardGammaDistribution::cdf(const double &x)
{
}

double StandardGammaDistribution::inv_cdf(const double &x)
{
}

double StandardGammaDistribution::mean() const { return alpha_; }

double StandardGammaDistribution::var() const { return alpha_; }

double StandardGammaDistribution::stdev() const { return sqrt(alpha_); }

void StandardGammaDistribution::random_draws(const std::vector<double> &uniform_draws,
                                         std::vector<double> &dist_draws)
{
}

#endif