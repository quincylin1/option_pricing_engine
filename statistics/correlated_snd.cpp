#ifndef __CORRELATED_SND_CPP
#define __CORRELATED_SND_CPP

#include "correlated_snd.h"
#include <iostream>
#include <cmath>

CorrelatedSND::CorrelatedSND(const double _rho,
                             const std::vector<double> *_uncorr_draws)
    : rho(_rho), uncorr_draws(_uncorr_draws){};

CorrelatedSND::~CorrelatedSND(){};

// This carries out the actual correlation modification. It is easy to see that if
// rho = 0.0 , then dist draws is unmodified , whereas if rho = 1.0 , then dist draws
// is simply set equal to uncorr draws. Thus with 0 < rho < 1 we have a 
// weighted average of each set .
void CorrelatedSND::correlation_calc(std::vector<double>& dist_draws){
    for (int i = 0; i < dist_draws.size(); i++){
        dist_draws[i] = rho * (*uncorr_draws)[i] + (sqrt(1 - rho * rho)) * dist_draws[i];
    }
}

void CorrelatedSND::random_draws(const std::vector<double>& uniform_draws,
                                 std::vector<double>& dist_draws){
// First check if the two vectors are of equal size
if (uniform_draws.size() != dist_draws.size()){
    std::cout << "The draws vectors are of unequal size" << std::endl;
    return;
}

// Check if size is even
if (uniform_draws.size() % 2 != 0){
    std::cout << "The size of draws vectors must be even" << std::endl;
    return;
}

// Box-Muller to generate gaussian random numbers from uniform random numbers
for (int i = 0; i < uniform_draws.size(); i++){
    dist_draws[2 * i] = sqrt(-2.0 * log(uniform_draws[2 * i])) *
                        sin(2 * M_PI * uniform_draws[2 * i + 1]);
    dist_draws[2 * i + 1] = sqrt(-2.0 * log(uniform_draws[2 * i])) *
                            cos(2 * M_PI * uniform_draws[2 * i + 1]);
}

// Convert gaussian dist_draws to correlated dist_draws
correlation_calc(dist_draws);

return;
}

#endif 