
#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <numeric>
// #include <arm_neon.h>
#include <sys/time.h>

#define ui64 uint64_t

inline double dml_micros() {
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return (tv.tv_sec * 1e6) + tv.tv_usec;
}

class FastRNG {
public:
    FastRNG() : generator(std::random_device{}()), distribution(0.0, 1.0) {}

    inline double next() {
        return distribution(generator);
    }

private:
    std::mt19937_64 generator;
    std::normal_distribution<double> distribution;
};

double black_scholes_monte_carlo(const double f1, const double f2, const double f3, const double f4, const ui64 K, const ui64 num_simulations) {
    double sum_payoffs = 0.0;
    #pragma omp parallel
    {
        FastRNG rng_local;

        double local_sum = 0.0;

        #pragma omp for
        for (ui64 i = 0; i < num_simulations; ++i) {
            double Z  = rng_local.next();
            if (Z <= f4) continue;
            double ST = f1 * std::exp(f2 * Z) - K;
            local_sum += ST;
        }

        #pragma omp atomic
        sum_payoffs += local_sum;
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

    std::cout << "Global initial seed: " << global_seed << "  "
              << "argv[1]= " << argv[1] << "  "
              << "argv[2]= " << argv[2] << std::endl;

    alignas(64) std::vector<double> bms(num_runs);
    double t1 = dml_micros();

    for (ui64 run = 0; run < num_runs; ++run) {
        bms[run] = black_scholes_monte_carlo(factor1, factor2, factor3, factor4, K, num_simulations);
    }

    double t2 = dml_micros();

    double min_val = *std::min_element(bms.begin(), bms.end());
    double max_val = *std::max_element(bms.begin(), bms.end());
    double mean = std::accumulate(bms.begin(), bms.end(), 0.0) / bms.size();

    double variance = std::accumulate(bms.begin(), bms.end(), 0.0, 
        [mean](double acc, double x) { return acc + (x - mean) * (x - mean); }) / bms.size();
    double stddev = std::sqrt(variance);
    double dev_percent = (stddev * 100.0) / mean;

    std::cout << std::fixed << std::setprecision(6)
              << "Min: " << min_val << "\n"
              << "Max: " << max_val << "\n"
              << "Mean: " << mean << "\n"
              << "Standard Deviation: " << stddev << "\n"
              << "Standard Deviation (% of Mean): " << dev_percent << "%" << "\n"
              << "Time taken: " << (t2 - t1) / 1e6 << " seconds" << std::endl;

    return 0;
}