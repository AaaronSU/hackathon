#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
//#include <trng/yarn2.hpp>
//#include <trng/normal_dist.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <armpl.h>
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

double dml_micros() {
    static struct timeval tv;
    gettimeofday(&tv, nullptr);
    return (tv.tv_sec * 1000000.0) + tv.tv_usec;
}

// Standard normal random number generators
std::vector<double> gaussian_box_muller_std(ui64 num_simulations) {
    std::vector<double> random_numbers(num_simulations);
    #pragma omp parallel
    {
        std::mt19937 generator(std::random_device{}() + omp_get_thread_num());
        std::normal_distribution<double> distribution(0.0, 1.0);
        #pragma omp for
        for (ui64 i = 0; i < num_simulations; ++i) {
            random_numbers[i] = distribution(generator);
        }
    }
    return random_numbers;
}

std::vector<double> gaussian_box_muller_boost(ui64 num_simulations) {
    std::vector<double> random_numbers(num_simulations);
    #pragma omp parallel
    {
        boost::random::mt19937 generator(std::random_device{}() + omp_get_thread_num());
        boost::random::normal_distribution<double> distribution(0.0, 1.0);
        #pragma omp for
        for (ui64 i = 0; i < num_simulations; ++i) {
            random_numbers[i] = distribution(generator);
        }
    }
    return random_numbers;
}

/* std::vector<double> gaussian_box_muller_trng(ui64 num_simulations) {
    std::vector<double> random_numbers(num_simulations);
    #pragma omp parallel
    {
        std::mt19937 generator(std::random_device{}() + omp_get_thread_num());
        trng::normal_dist<double> distribution(0.0, 1.0);
        #pragma omp for
        for (ui64 i = 0; i < num_simulations; ++i) {
            random_numbers[i] = distribution(generator);
        }
    }
    return random_numbers;
} */

std::vector<double> gaussian_box_muller_gsl(ui64 num_simulations) {
    std::vector<double> random_numbers(num_simulations);
    #pragma omp parallel
    {
        gsl_rng* generator = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(generator, std::random_device{}() + omp_get_thread_num());
        #pragma omp for
        for (ui64 i = 0; i < num_simulations; ++i) {
            random_numbers[i] = gsl_ran_gaussian(generator, 1.0);
        }
        gsl_rng_free(generator);
    }
    return random_numbers;
}

std::vector<double> gaussian_box_muller_armpl(ui64 num_simulations) {
    std::vector<double> random_numbers(num_simulations);

    #pragma omp parallel
    {
        VSLStreamStatePtr stream;
        vslNewStream(&stream, VSL_BRNG_MT19937, std::random_device{}());

        #pragma omp for
        for (ui64 i = 0; i < num_simulations; ++i) {
            vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, 1, &random_numbers[i], 0.0, 1.0);
        }

        vslDeleteStream(&stream);
    }

    return random_numbers;
}

double black_scholes_monte_carlo(double f1, double f2, double f3, ui64 K, ui64 num_simulations, std::vector<double>& r) {
    double sum_payoffs = 0.0;
    #pragma omp parallel for reduction(+:sum_payoffs)
    for (ui64 i = 0; i < num_simulations; ++i) {
        double Z = r[i];
        double ST = f1 * exp(f2 * Z);
        double payoff = (ST > K) * (ST - K);
        sum_payoffs += payoff;
    }
    return f3 * (sum_payoffs / num_simulations);
}

void run_simulation(const std::string& generator_name, std::function<std::vector<double>(ui64)> gaussian_generator, double factor1, double factor2, double factor3, ui64 K, ui64 num_simulations, ui64 num_runs) {

    std::vector<double> bms(num_runs);

    // Pre-generate all random numbers outside the simulation loop
    std::vector<std::vector<double>> all_random_numbers(num_runs);
    for (ui64 run = 0; run < num_runs; ++run) {
        all_random_numbers[run] = gaussian_generator(num_simulations);
    }

    double t1 = dml_micros();
    #pragma omp parallel
    {
        std::vector<double> r;
        #pragma omp for schedule(dynamic)
        for (ui64 run = 0; run < num_runs; ++run) {
            // Use pre-generated random numbers
            r = all_random_numbers[run];
            double result = black_scholes_monte_carlo(factor1, factor2, factor3, K, num_simulations, r);
            bms[run] = result;
        }
    }
    double t2 = dml_micros();

    double min_val = *std::min_element(bms.begin(), bms.end());
    double max_val = *std::max_element(bms.begin(), bms.end());
    double mean = std::accumulate(bms.begin(), bms.end(), 0.0) / bms.size();
    double variance = 0.0;
    for (double x : bms) {
        variance += (x - mean) * (x - mean);
    }
    variance /= bms.size();
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
    ui64 num_runs = std::stoull(argv[2]);

    constexpr ui64 S0 = 100;
    constexpr ui64 K = 110;
    constexpr double T = 1.0;
    constexpr double r = 0.06;
    constexpr double sigma = 0.2;
    constexpr double q = 0.03;

    const double factor1 = S0 * exp((r - q - 0.5 * sigma * sigma) * T);
    const double factor2 = sigma * sqrt(T);
    const double factor3 = exp(-r * T);

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

    std::vector<std::pair<std::string, std::function<std::vector<double>(ui64)>>> generators = {
        {"STD", gaussian_box_muller_std},
        //{"TRNG", gaussian_box_muller_trng},
        {"Boost", gaussian_box_muller_boost},
        {"GSL", gaussian_box_muller_gsl},
        {"ARMPL", gaussian_box_muller_armpl}
    };

    for (const auto& [name, generator] : generators) {
        run_simulation(name, generator, factor1, factor2, factor3, K, num_simulations, num_runs);
    }

    return 0;
}