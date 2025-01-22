#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <omp.h>
#include <arm_sve.h>
#include <sleef.h>

#define ui64 u_int64_t

#include <sys/time.h>
double dml_micros() {
    static struct timezone tz;
    static struct timeval tv;
    gettimeofday(&tv, &tz);
    return ((tv.tv_sec * 1000000.0) + tv.tv_usec);
}

// Vectorized Gaussian random number generation using Box-Muller transform
svfloat64_t gaussian_box_muller_sve(svbool_t pg, svfloat64_t u1, svfloat64_t u2) {
    svfloat64_t log_u1 = Sleef_finz_logdx_u10sve(u1);
    // svfloat64_t sqrt_term = Sleef_finz_powdx_u10sve(svmul_f64_z(pg, svdup_f64(-2.0), log_u1), svdup_f64(0.5));
    svfloat64_t sqrt_term = svsqrt_f64_z(pg, svmul_f64_z(pg, svdup_f64(-2.0), log_u1));
    svfloat64_t theta = svmul_f64_z(pg, svdup_f64(2.0 * M_PI), u2);
    return svmul_f64_z(pg, sqrt_term, Sleef_finz_cosdx_u10sve(theta));
}

// Vectorized Monte Carlo Black-Scholes simulation
double black_scholes_monte_carlo_sve(double f1, double f2, double f3, double f4, ui64 K, ui64 num_simulations) {
    double sum_payoffs = 0.0;
    ui64 i = 0;
    svbool_t pg = svwhilelt_b64(i, num_simulations);
    

    std::mt19937_64 rng(std::random_device{}());
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

    while (svptest_any(svptrue_b64(), pg)) {
        // Generate uniform random numbers
        double u1[svcntd()];
        double u2[svcntd()];
        for (int j = 0; j < svcntd(); ++j) {
            u1[j] = uniform_dist(rng);
            u2[j] = uniform_dist(rng);
        }

        svfloat64_t sv_u1 = svld1(pg, u1);
        svfloat64_t sv_u2 = svld1(pg, u2);

        // Convert to Gaussian random numbers
        svfloat64_t Z = gaussian_box_muller_sve(pg, sv_u1, sv_u2);

        // Compute ST = f1 * exp(f2 * Z)
        //svfloat64_t ST = svmul_f64_z(pg, svdup_f64(f1), Sleef_finz_expdx_u10sve(svmul_f64_z(pg, svdup_f64(f2), Z)));

        // Compute payoff = max(ST - K, 0.0)
        //svfloat64_t payoff = svmax_f64_z(pg, svsub_f64_z(pg, ST, svdup_f64(K)), svdup_f64(0.0));

       
        // double payoff = (Z > f4) * (f1 * exp(f2 * Z) - K);
        /*
        svfloat64_t payoff = svmul_f64_z(
            svcmpgt_f64(pg, Z, svdup_f64(f4)), // predicate: Z <= f4 then return 0 
            svdup_f64(1.0),
            svsub_f64_z(pg, svmul_f64_z(pg, svdup_f64(f1), Sleef_finz_expdx_u10sve(svmul_f64_z(pg, svdup_f64(f2), Z))), svdup_f64(K)));
        */

        // Accumulate payoffs
        // sum_payoffs += svaddv_f64(pg, payoff);
        svbool_t pd = svcmpgt_f64(pg, Z, svdup_f64(f4)); // predicate: Z <= f4 then return 0
        sum_payoffs += svaddv_f64(pd, svsub_f64_z(pg, svmul_f64_z(pg, svdup_f64(f1), Sleef_finz_expdx_u10sve(svmul_f64_z(pg, svdup_f64(f2), Z))), svdup_f64(K)));


        // Move to the next batch
        i += svcntd();
        pg = svwhilelt_b64(i, num_simulations);
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

    ui64 S0 = 100;
    ui64 K = 110;
    double T = 1.0;
    double r = 0.06;
    double sigma = 0.2;
    double q = 0.03;

    double factor1 = S0 * exp((r - q - 0.5 * sigma * sigma) * T);
    double factor2 = sigma * sqrt(T);
    double factor3 = exp(-r * T);
    double factor4 = log(K / factor1) / factor2;

    double sum = 0.0;
    double t1 = dml_micros();
    #pragma omp parallel for reduction(+:sum)
    std::vector<double> bms;
    for (ui64 run = 0; run < num_runs; ++run) {
        double result = black_scholes_monte_carlo(factor1, factor2, factor3, factor4, K, num_simulations);
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

    return 0;
}
