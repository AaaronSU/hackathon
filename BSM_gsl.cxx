#include <iostream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iomanip>   // For setting precision
#include <sys/time.h>
#include <random>

#define ui64 u_int64_t

// Function to get current time in microseconds
double dml_micros() {
    static struct timezone tz;
    static struct timeval tv;
    gettimeofday(&tv, &tz);
    return ((tv.tv_sec * 1000000.0) + tv.tv_usec);
}

// Function to generate Gaussian noise using Box-Muller transform
double gaussian_box_muller(gsl_rng* rng) {
    return gsl_ran_gaussian(rng, 1.0); // 1.0 is the standard deviation (mean=0)
}

// Function to calculate the Black-Scholes call option price using Monte Carlo method
double black_scholes_monte_carlo(double f1, double f2, double f3, ui64 K, ui64 num_simulations, gsl_rng* rng) {
    double sum_payoffs = 0.0;
    for (ui64 i = 0; i < num_simulations; ++i) {
        double Z = gaussian_box_muller(rng);
        double ST = f1 * exp(f2 * Z);
        double payoff = std::max(ST - K, 0.0); // Call option payoff
        sum_payoffs += payoff;
    }
    return f3 * (sum_payoffs / num_simulations);
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <num_simulations> <num_runs>" << std::endl;
        return 1;
    }

    ui64 num_simulations = std::stoull(argv[1]);
    ui64 num_runs = std::stoull(argv[2]);

    // Input parameters
    ui64 S0 = 100;                   // Initial stock price
    ui64 K = 110;                    // Strike price
    double T = 1.0;                   // Time to maturity (1 year)
    double r = 0.06;                  // Risk-free interest rate
    double sigma = 0.2;               // Volatility
    double q = 0.03;                  // Dividend yield

    double factor1 = S0 * exp((r - q - 0.5 * sigma * sigma) * T);
    double factor2 = sigma * sqrt(T);
    double factor3 = exp(-r * T);

    // Initialize GSL random number generator
    const gsl_rng_type* T_rng;
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937); // Mersenne Twister RNG
    gsl_rng_set(rng, std::random_device{}()); // Set random seed from random_device

    std::vector<double> bms;
    double t1 = dml_micros();
    for (ui64 run = 0; run < num_runs; ++run) {
        double result = black_scholes_monte_carlo(factor1, factor2, factor3, K, num_simulations, rng);
        bms.push_back(result);
    }
    double t2 = dml_micros();

    // Calculate statistics
    double min_val = *std::min_element(bms.begin(), bms.end());
    double max_val = *std::max_element(bms.begin(), bms.end());
    double mean = std::accumulate(bms.begin(), bms.end(), 0.0) / bms.size();

    double variance = 0.0;
    for (double value : bms) {
        variance += (value - mean) * (value - mean);
    }
    variance /= bms.size();
    double stddev = std::sqrt(variance);

    double dev_percent = (stddev * 100.0) / mean;

    std::cout << std::fixed << std::setprecision(6)
              << "Min: " << min_val << "\n"
              << "Max: " << max_val << "\n"
              << "Mean: " << mean << "\n"
              << "Standard Deviation: " << stddev << "\n"
              << "Standard Deviation (% of Mean): " << dev_percent << "%" << "\n"
              << "Time taken: " << (t2 - t1) / 1000000.0 << " seconds" << std::endl;

    // Free GSL RNG memory
    gsl_rng_free(rng);

    return 0;
}
