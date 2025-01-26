#include <iostream>
#include <armpl.h> // ARM Performance Libraries for random number generation
#include <openrng.h>
#include <sys/time.h>
#include <iomanip>

#define ui64 uint64_t // Alias for unsigned 64-bit integer

// Function to get current time in microseconds
double dml_micros() {
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return (tv.tv_sec * 1e6) + tv.tv_usec;
}

// Function to check for errors in MKL operations
void check_error(int errcode, const char* message) {
    if (errcode != VSL_ERROR_OK) {
        std::cerr << "Error: " << message << " (" << errcode << ")" << std::endl;
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char** argv) {
    // Ensure the program is run with the correct number of arguments
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <num_simulations> <num_runs>" << std::endl;
        return 1;
    }

    // Read input arguments for the number of simulations and runs
    ui64 num_simulations = std::stoull(argv[1]);
    ui64 num_runs        = std::stoull(argv[2]);

    // Define constants for the simulation
    constexpr double S0    = 100.0; // Initial stock price
    constexpr double K     = 110.0; // Strike price
    constexpr double T     = 1.0;   // Time to maturity in years
    constexpr double r     = 0.06;  // Risk-free interest rate
    constexpr double sigma = 0.2;   // Volatility of the stock
    constexpr double q     = 0.03;  // Dividend yield

    // Derived parameters for the lognormal distribution
    constexpr double f1 = S0 * exp((r - q - 0.5 * sigma * sigma) * T); // Adjusted stock price factor
    constexpr double f2 = sigma * sqrt(T); // Standard deviation of the underlying normal distribution
    constexpr double f3 = exp(-r * T);     // Discount factor for present value
    constexpr double f4 = log(K / f1) / f2; // Threshold Z value for truncation

    int n = num_simulations; // Number of simulations

    // Parameters for the lognormal random number generator
    double alpha = 0.0;  // Mean of the underlying normal distribution
    double sig = f2;     // Standard deviation of the underlying normal distribution
    double b = K;        // Displacement (minimum value for the distribution)
    double beta = 1 / f4; // Scale factor for the distribution (adjusted inversely to f4)

    // Allocate memory for storing generated random numbers
    double* vec = (double*)malloc(sizeof(double) * n);

    double sum = 0.0; // Accumulator for the sum of generated numbers

    // Measure start time
    double t1 = dml_micros();

    // Initialize the random number stream using Philox random number generator
    VSLStreamStatePtr stream;
    int errcode = vslNewStream(&stream, VSL_BRNG_PHILOX4X32X10, 42);
    check_error(errcode, "vslNewStream failed");

    // Generate lognormal random numbers
    errcode = vdRngLognormal(VSL_RNG_METHOD_LOGNORMAL_BOXMULLER2, stream, n, vec, alpha, sig, b, beta);
    check_error(errcode, "vdRngLognormal failed");

    // Compute the sum of the generated numbers using parallel reduction
    #pragma omp parallel for reduction(+:sum)
    for (size_t i = 0; i < n; i++) {
        sum += vec[i];
    }

    // Output the mean of the generated random numbers
    std::cout << sum / n << std::endl;

    // Measure end time
    double t2 = dml_micros();

    // Clean up the random number stream
    errcode = vslDeleteStream(&stream);
    check_error(errcode, "vslDeleteStream failed");

    // Output the time taken for the simulation
    std::cout << std::fixed << std::setprecision(6)
              << "Time taken: " << (t2 - t1) / 1e6 << " seconds" << std::endl;

    // Free allocated memory
    free(vec);

    return 0;
}
