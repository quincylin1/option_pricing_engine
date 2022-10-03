#ifndef __HESTON_EULER_HPP
#define __HESTON_EULER_HPP

#include "../options/option.hpp"
#include "../statistics/statistics.hpp"


class HestonEuler {
protected:
    unsigned _num_sims;
    unsigned _num_intervals;

    Option *_pOption;

    double _kappa;
    double _theta;
    double _xi;
    double _rho;

public:
    HestonEuler();
    HestonEuler(unsigned num_sims,
                unsigned num_intervals,
                Option* pOption,
                double kappa, double theta,
                double xi, double rho);
    virtual ~HestonEuler();

    virtual void generate_normal_correlated_paths(std::vector<double> &spot_normals,
                                                  std::vector<double> &cor_normals);

    virtual double monte_carlo_sim(const double& S_0) = 0;

    // Calculate the volatility path
    virtual void calculate_vol_path(const std::vector<double> &vol_draws,
                       std::vector<double> &vol_path) = 0;

    // Calculate the asset price path
    virtual void calculate_spot_path(const std::vector<double> &spot_draws,
                        const std::vector<double> &vol_path,
                        std::vector<double> &spot_path) = 0;
};

class HestonEulerReflection : public HestonEuler
{
public:
    HestonEulerReflection(unsigned num_sims,
                          unsigned num_intervals,
                          Option *pOption,
                          double kappa, double theta,
                          double xi, double rho);
    virtual ~HestonEulerReflection();

    // Calculate the volatility path
    virtual void calculate_vol_path(const std::vector<double> &vol_draws,
                                    std::vector<double> &vol_path);

    // Calculate the asset price path
    virtual void calculate_spot_path(const std::vector<double> &spot_draws,
                                     const std::vector<double> &vol_path,
                                     std::vector<double> &spot_path);
};

class HestonEulerPartialTruncation : public HestonEuler
{
public:
    HestonEulerPartialTruncation(unsigned num_sims,
                                 unsigned num_intervals,
                                 Option *pOption,
                                 double kappa, double theta,
                                 double xi, double rho);
    virtual ~HestonEulerPartialTruncation();

    // Calculate the volatility path
    virtual void calculate_vol_path(const std::vector<double> &vol_draws,
                                    std::vector<double> &vol_path);

    // Calculate the asset price path
    virtual void calculate_spot_path(const std::vector<double> &spot_draws,
                                     const std::vector<double> &vol_path,
                                     std::vector<double> &spot_path);
};

class HestonEulerFullTruncation : public HestonEuler
{
public:
    HestonEulerFullTruncation(unsigned num_sims,
                              unsigned num_intervals,
                              Option *pOption,
                              double kappa, double theta,
                              double xi, double rho);
    virtual ~HestonEulerFullTruncation();

    // Calculate the volatility path
    virtual void calculate_vol_path(const std::vector<double> &vol_draws,
                                    std::vector<double> &vol_path);

    // Calculate the asset price path
    virtual void calculate_spot_path(const std::vector<double> &spot_draws,
                                     const std::vector<double> &vol_path,
                                     std::vector<double> &spot_path);
};

#endif