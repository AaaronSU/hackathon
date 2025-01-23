/* 

    Monte Carlo Hackathon created by Hafsa Demnati and Patrick Demichel @ Viridien 2024
    The code compute a Call Option with a Monte Carlo method and compare the result with the analytical equation of Black-Scholes Merton : more details in the documentation

    Compilation : g++ -O BSM.cxx -o BSM

    Exemple of run: ./BSM #simulations #runs

./BSM 100 1000000
Global initial seed: 21852687      argv[1]= 100     argv[2]= 1000000
 value= 5.136359 in 10.191287 seconds

./BSM 100 1000000
Global initial seed: 4208275479      argv[1]= 100     argv[2]= 1000000
 value= 5.138515 in 10.223189 seconds
 
   We want the performance and value for largest # of simulations as it will define a more precise pricing
   If you run multiple runs you will see that the value fluctuate as expected
   The large number of runs will generate a more precise value then you will converge but it require a large computation

   give values for ./BSM 100000 1000000        
               for ./BSM 1000000 1000000
               for ./BSM 10000000 1000000
               for ./BSM 100000000 1000000
  
   We give points for best performance for each group of runs 

   You need to tune and parallelize the code to run for large # of simulations

*/

#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <limits>
#include <algorithm>
#include <iomanip>   // For setting precision
#include <armpl.h>

#define ui64 u_int64_t

#include <sys/time.h>
double dml_micros() {
    static struct timezone tz;
    static struct timeval tv;
    gettimeofday(&tv, &tz);
    return ((tv.tv_sec * 1000000.0) + tv.tv_usec);
}

// Fonction de simulation de Monte Carlo avec Black-Scholes utilisant MKL
double black_scholes_monte_carlo_mkl(double f1, double f2, double f3, ui64 K, ui64 num_simulations) {
    double sum_payoffs = 0.0;

    VSLStreamStatePtr stream;
    int errcode = vslNewStream(&stream, VSL_BRNG_MT19937, std::random_device{}());
    
    if (errcode != VSL_STATUS_OK) {
        std::cerr << "Erreur lors de la création du flux de nombres aléatoires!" << std::endl;
        return -1.0;
    }

    std::vector<double> r(num_simulations);
    vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, num_simulations, r.data(), 0.0, 1.0);
    
    // Calcul des payoffs
    //#pragma omp parallel for reduction(+:sum_payoffs)
    //#pragma GCC unroll 8
    for (ui64 i = 0; i < num_simulations; ++i) {
        double Z = r[i]; 
        //std::cout << Z << std::endl;
        double ST = f1 * exp(f2 * Z); // Calcul de ST = f1 * exp(f2 * Z)
        double payoff = (ST > K) * (ST - K); // Payoff de l'option
        sum_payoffs += payoff;
    }

    vslDeleteStream(&stream);

    return f3 * (sum_payoffs / num_simulations);
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <num_simulations> <num_runs>" << std::endl;
        return 1;
    }

    ui64 num_simulations = std::stoull(argv[1]);
    ui64 num_runs = std::stoull(argv[2]);

    // Paramètres de l'option
    ui64 S0 = 100;
    ui64 K = 110;
    double T = 1.0;
    double r = 0.06;
    double sigma = 0.2;
    double q = 0.03;

    // Paramètres précalculés
    double factor1 = S0 * exp((r - q - 0.5 * sigma * sigma) * T);
    double factor2 = sigma * sqrt(T);
    double factor3 = exp(-r * T);

    std::vector<double> bms(num_runs);
    double t1 = dml_micros();

    // Lancer les simulations en parallèle
    //#pragma omp parallel for
    for (ui64 run = 0; run < num_runs; ++run) {
        bms[run] = black_scholes_monte_carlo_mkl(factor1, factor2, factor3, K, num_simulations);
    }

    double t2 = dml_micros();

    std::random_device rd;
    unsigned long long global_seed = rd();  // This will be the global seed

    std::cout << "Global initial seed: " << global_seed << "      argv[1]= " << argv[1] << "     argv[2]= " << argv[2] <<  std::endl;

    // Calculer les statistiques
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