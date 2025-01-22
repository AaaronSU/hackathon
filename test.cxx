#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <omp.h>
#include <sys/time.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using ui64 = unsigned long long;
gsl_rng* rng;

// Function to get the current time in microseconds
double
dml_micros()
{
        static struct timezone tz;
        static struct timeval  tv;
        gettimeofday(&tv,&tz);
        return((tv.tv_sec*1000000.0)+tv.tv_usec);
}

// Function to generate Gaussian noise using Box-Muller transform
double gaussian_box_muller() {
    return gsl_ran_gaussian(rng, 1.0);
}

// Function to calculate the exp values table for Z in range [-2.576, 2.576] (> 99%)
// with better accuracy for Z in range [-1.281, 1.281] (> 80%)
void exp_table(double f2, int n1, int n2, std::vector<double>& table_x, std::vector<double>& table_y, double h1, double h2)
{
    #pragma omp parallel for
    for (int i = 0; i < n1 + 1; i++) 
    {
        table_x[i] = -2.576 + i * h1;
        table_y[i] = exp((-2.576 + i * h1) * f2);
        table_x[i + n1 + n2] = 1.281 + i * h1;
        table_y[i + n1 + n2] = exp((1.281 + i * h1) * f2);
    }

    #pragma omp parallel for 
    for (int i = 1; i < n2; i++) 
    {
        table_x[i + n1] = -1.281 + i * h2;
        table_y[i + n1] = exp((-1.281 + i * h2) * f2);
    }
}

// Function to check approximate exp value in the table
double exp_fast(double f2, double Z, int n1, int n2, const std::vector<double>& table_x, const std::vector<double>& table_y, double inv_h1, double inv_h2)
{
    int i = 0;

    // For Z in range [-2.576, -1.281], use h1
    // For Z in range [1.281, 2.576], use h1
    // For Z in range [-1.281, 1.281], use h2
    i = (Z <= -1.281) ? (Z + 2.576) * inv_h1 :
        (Z >= 1.281) ? n1 + n2 + (Z - 1.281) * inv_h1 :
        n1 + (Z + 1.281) * inv_h2;

    if (Z >= -2.576 && Z <= 2.576) 
    {
        double a1 = table_x[i];
        double a2 = table_x[i + 1];
        double b1 = table_y[i];
        double b2 = table_y[i + 1];
        double tmp = b1 * f2 * (Z - a1) + b1 + b2 * f2 * (Z - a2) + b2;
        tmp *= 0.5;
        return tmp;
    }
    else 
    {
        return exp(Z * f2);
    }
}

// Function to calculate the Black-Scholes call option price using Monte Carlo method
double black_scholes_monte_carlo(double f1, double f2, double f3, ui64 K, ui64 num_simulations, int n1, int n2, const std::vector<double>& table_x, const std::vector<double>& table_y, double inv_h1, double inv_h2) {
    double sum_payoffs = 0.0;

    #pragma omp parallel for reduction(+:sum_payoffs)
    for (ui64 i = 0; i < num_simulations; ++i) 
    {
        double Z = gaussian_box_muller();
        double exp_output = exp_fast(f2, Z, n1, n2, table_x, table_y, inv_h1, inv_h2);
        
        double ST = exp_output * f1;
        double payoff = (ST > K) * (ST - K);
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
    ui64 num_runs        = std::stoull(argv[2]);

    // Input parameters
    ui64 S0      = 100;                   // Initial stock price
    ui64 K       = 110;                   // Strike price
    double T     = 1.0;                   // Time to maturity (1 year)
    double r     = 0.06;                  // Risk-free interest rate
    double sigma = 0.2;                   // Volatility
    double q     = 0.03;                  // Dividend yield

    // Initialize GSL random number generator
    //const gsl_rng_type* T_rng;
    rng = gsl_rng_alloc(gsl_rng_mt19937); // Mersenne Twister RNG
    gsl_rng_set(rng, std::random_device{}());

    std::cout << "Global initial seed: " << gsl_rng_default_seed << "      argv[1]= " << argv[1] << "     argv[2]= " << argv[2] <<  std::endl;

    // Precompute the factors
    double factor1 = S0 * exp((r - q - 0.5 * sigma * sigma) * T);
    double factor2 = sigma * sqrt(T);
    double factor3 = exp(-r * T);

    // Precompute the exp table
    double h1 = 1e-6;
    double h2 = 1e-8;
    double inv_h1 = 1.0 / h1;
    double inv_h2 = 1.0 / h2;
    int n1 = (2.576 - 1.281) / h1;
    int n2 = (1.281 + 1.281) / h2;
    int size = n1 + n2 + n1 + 1;
    std::vector<double> table_x(size);
    std::vector<double> table_y(size);
    exp_table(factor2, n1, n2, table_x, table_y, h1, h2);

    std::vector<double> bms;
    //double sum=0.0;
    double t1=dml_micros();
    for (ui64 run = 0; run < num_runs; ++run) {
        double value = black_scholes_monte_carlo(factor1, factor2, factor3, K, num_simulations, n1, n2, table_x, table_y, inv_h1, inv_h2);
        bms.push_back(value);
    }
    double t2=dml_micros();

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

    //std::cout << std::fixed << std::setprecision(6) << " value= " << sum/num_runs << " in " << (t2-t1)/1000000.0 << " seconds" << std::endl;

    return 0;
}
