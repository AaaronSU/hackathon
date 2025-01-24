// Purpose: To compare the performance of different random number generators in the context of the Black-Scholes-Merton (BSM) model.
// Libs: boost, trng, gsl, armpl
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <trng/yarn2.hpp>
#include <trng/normal_dist.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <armpl.h>

// Standard C++ headers
#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <limits>
#include <algorithm>
#include <iomanip>  
#include <sys/time.h>
#include <functional>
#include <omp.h>

using ui64 = uint64_t;

//
double dml_micros() {
    static struct timeval tv;
    gettimeofday(&tv, nullptr);
    return((tv.tv_sec*1000000.0)+tv.tv_usec);
}

// Gaussian random number generator using C++11 random library with mt19937
double gaussian_box_muller_std() {
    static thread_local std::mt19937 generator(std::random_device{}());
    static thread_local std::normal_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
}

// Gaussian random number generator using Boost library with mt19937
double gaussian_box_muller_boost() {
    static thread_local boost::random::mt19937 generator(std::random_device{}());
    static thread_local boost::random::normal_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
}

// Gaussian random number generator using TRNG library with mt19937
double gaussian_box_muller_trng() {
    static thread_local std::mt19937 generator(std::random_device{}());
    static thread_local trng::normal_dist<double> distribution(0.0, 1.0);
    return distribution(generator);
}

// Gaussian random number generator using GSL library with mt19937
double gaussian_box_muller_gsl() {
    static thread_local gsl_rng* generator = gsl_rng_alloc(gsl_rng_mt19937);

    static thread_local bool is_seeded = []() {
        gsl_rng_set(generator, std::random_device{}() + omp_get_thread_num());
        return true;
    }();

    return gsl_ran_gaussian(generator, 1.0);
}

// Gaussian random number generator using Armpl library with mt19937
double gaussian_box_muller_armpl() {
    static thread_local VSLStreamStatePtr stream;
    static thread_local bool is_seeded = []() {
        int errcode = vslNewStream(&stream, VSL_BRNG_MT19937, std::random_device{}());
        if (errcode != VSL_STATUS_OK) {
            std::cerr << "Error creating random number stream!" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        return true;
    }();

    double result;
    vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, 1, &result, 0.0, 1.0);
    return result;
}

// Function to calculate the Black-Scholes call option price using Monte Carlo method
double black_scholes_monte_carlo(double f1, double f2, double f3, ui64 K, ui64 num_simulations, std::function<double()> gaussian_generator) {
    double sum_payoffs = 0.0;
    
    double Z, ST, payoff;
    #pragma omp parallel for reduction(+:sum_payoffs)
    for (ui64 i = 0; i < num_simulations; ++i) {
        Z = gaussian_generator();
        ST = f1 * exp(f2 * Z);
        payoff = (ST > K) * (ST - K);
        sum_payoffs += payoff;
    }
    return f3 * (sum_payoffs / num_simulations);
}

void run_simulation(const std::string& generator_name, std::function<double()> gaussian_generator, double factor1, double factor2, double factor3, ui64 K, ui64 num_simulations, ui64 num_runs) {
    std::vector<double> bms;
    bms.reserve(num_runs);
    double t1 = dml_micros();

    #pragma omp parallel for 
    for (ui64 run = 0; run < num_runs; ++run) {
        double result = black_scholes_monte_carlo(factor1, factor2, factor3, K, num_simulations, gaussian_generator);
        bms.push_back(result);
    }
    double t2 = dml_micros();

    double min_val = *std::min_element(bms.begin(), bms.end());
    double max_val = *std::max_element(bms.begin(), bms.end());
    double mean = 0.0;
    double variance = 0.0;
    for (ui64 i = 0; i < bms.size(); ++i) {
        double delta = bms[i] - mean;
        mean += delta / (i + 1);
        variance += delta * (bms[i] - mean);
    }
    variance /= (bms.size() - 1);
    double stddev = std::sqrt(variance);

    double dev_percent = (mean != 0.0) ? (stddev * 100.0) / mean : 0.0;

    std::cout << std::setw(20) << generator_name
          << std::setw(20) << min_val
          << std::setw(20) << max_val
          << std::setw(20) << mean
          << std::setw(20) << stddev
          << std::setw(20) << dev_percent
          << std::setw(20) << (t2 - t1) / 1000000.0 << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <num_simulations> <num_runs>" << std::endl;
        return 1;
    }

    ui64 num_simulations = std::stoull(argv[1]);
    ui64 num_runs        = std::stoull(argv[2]);

    constexpr ui64 S0 = 100;
    constexpr ui64 K = 110;
    constexpr double T = 1.0;
    constexpr double r = 0.06;
    constexpr double sigma = 0.2;
    constexpr double q = 0.03;

    constexpr double factor1 = S0 * exp((r - q - 0.5 * sigma * sigma) * T);
    constexpr double factor2 = sigma * sqrt(T);
    constexpr double factor3 = exp(-r * T);
    constexpr double factor4 = log(K / factor1) / factor2;

    std::random_device rd;
    unsigned long long global_seed = rd();

    std::cout << std::fixed << std::setprecision(6);

    // Display the header
    std::cout << "------------------------------------------" << std::endl; 
    std::cout << "Global initial seed: " << global_seed << "      argv[1]= " << argv[1] << "     argv[2]= " << argv[2] <<  std::endl;
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "Statistics Results:" << std::endl;
    std::cout << "------------------------------------------" << std::endl;
    std::cout << std::setw(20) << "Generator"
              << std::setw(20) << "Min"
              << std::setw(20) << "Max"
              << std::setw(20) << "Mean"
              << std::setw(20) << "Stddev"
              << std::setw(20) << "COV(%)"
              << std::setw(20) << "Time Taken (s)" << std::endl;
    
    std::vector<std::pair<std::string, std::function<double()>>> generators = {
        {"STD", gaussian_box_muller_std},
        {"TRNG", gaussian_box_muller_trng},
        {"Boost", gaussian_box_muller_boost},
        {"ARMPL", gaussian_box_muller_armpl},
        {"GSL", gaussian_box_muller_gsl}

    };

    for (const auto& [name, generator] : generators) {
        run_simulation(name, generator, factor1, factor2, factor3, K, num_simulations, num_runs);
    }

    return 0;
}
