#include <iostream>
#include <chrono>
#include <vector>
#include <random>
#include <cmath>
#include <numc>

using namespace numc;
using namespace std::chrono;

void run_matrix_benchmark() {
    std::cout << "--- Matrix Multiplication Benchmark ---\n";
    auto start = high_resolution_clock::now();
    
    // Matrix size 1000x1000
    int N = 1000;
    tensor A = tensor<>::eye(N);
    tensor B = tensor<>::eye(N);
    
    auto mid = high_resolution_clock::now();
    
    tensor C = A * B;
    
    auto end = high_resolution_clock::now();
    
    std::cout << "Initialization (" << N << "x" << N << "): " 
              << duration_cast<milliseconds>(mid - start).count() << " ms\n";
    std::cout << "Multiplication (" << N << "x" << N << "): " 
              << duration_cast<milliseconds>(end - mid).count() << " ms\n\n";
}

void run_lu_benchmark() {
    std::cout << "--- LU Decomposition Benchmark ---\n";
    int N = 500;
    tensor A = tensor<>::eye(N);
    // Make it not perfectly identity to avoid trivial zeros
    for(int i=0; i<N; ++i) A(i, i) = 2.0;
    
    auto start = high_resolution_clock::now();
    tensor LU = linear_algebra::lu_decomp(A);
    auto end = high_resolution_clock::now();
    
    std::cout << "LU Decomposition (" << N << "x" << N << "): " 
              << duration_cast<milliseconds>(end - start).count() << " ms\n\n";
}

void run_fft_benchmark() {
    std::cout << "--- FFT Benchmark ---\n";
    int N = 1048576; // 2^20 (roughly 1 million points)
    numc::vector<double> signal(N);
    for(int i = 0; i < N; ++i) signal[i] = std::sin(2.0 * PI * i / 100.0);
    
    auto start = high_resolution_clock::now();
    
    auto spectrum = analysis::fft(signal);
    
    auto end = high_resolution_clock::now();
    
    std::cout << "FFT (" << N << " points): " 
              << duration_cast<milliseconds>(end - start).count() << " ms\n\n";
}

void run_calculus_benchmark() {
    std::cout << "--- Calculus & Integration Benchmark ---\n";
    auto x = funcv::x;
    func f = sin(x) * cos(x) + exp(-0.1 * x);
    
    // 1. Function Evaluation Benchmark
    auto start = high_resolution_clock::now();
    double sum = 0.0;
    int steps = 1000000; // 1M evaluations
    for(int i = 0; i < steps; ++i) {
        sum += f(i * 0.0001);
    }
    auto end = high_resolution_clock::now();
    std::cout << "Function Evaluation (1M eval): " 
              << duration_cast<milliseconds>(end - start).count() << " ms (sum check: " << sum << ")\n";
              
    // 2. Integration Benchmark
    start = high_resolution_clock::now();
    double area = analysis::simpson_1_3(f, 0.0, 100.0);
    end = high_resolution_clock::now();
    
    std::cout << "Numerical Integration (Simpson, [0, 100]): " 
              << duration_cast<milliseconds>(end - start).count() << " ms\n\n";
}

void run_rk4_benchmark() {
    std::cout << "--- Differential Equations (RK4) Benchmark ---\n";
    // 3-variable coupled ODE system
    std::function<numc::vector<double>(double, const numc::vector<double>&)> F = 
        [](double t, const numc::vector<double>& y) -> numc::vector<double> {
            return {-0.1 * y[0], -0.2 * y[1], -0.3 * y[2]};
        };
    
    auto start = high_resolution_clock::now();
    
    // 100,000 steps (t_end=100.0, h=0.001)
    numc::vector<double> y0 = {100.0, 50.0, 25.0};
    auto result = differential_equations::run_kut4(F, 0.0, y0, 100.0, 0.001);
    
    auto end = high_resolution_clock::now();
    
    std::cout << "Runge-Kutta 4th Order (100k steps, 3-eq system): " 
              << duration_cast<milliseconds>(end - start).count() << " ms\n\n";
}

void run_kmeans_benchmark() {
    std::cout << "--- Machine Learning (K-Means) Benchmark ---\n";
    int N = 10000;
    std::vector<vector<double>> points(N, vector<double>(2));
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(-10.0, 10.0);
    
    for(int i = 0; i < N; ++i) {
        points[i][0] = dist(gen);
        points[i][1] = dist(gen);
    }
    
    auto start = high_resolution_clock::now();
    
    auto result = machine_learning::kmeans(points, 5); // K=5
    
    auto end = high_resolution_clock::now();
    
    std::cout << "K-Means Clustering (" << N << " pts, K=5): " 
              << duration_cast<milliseconds>(end - start).count() << " ms\n\n";
}

void run_primes_benchmark() {
    std::cout << "--- Number Theory (Primes) Benchmark ---\n";
    int N = 10000000; // 10 million
    auto start = high_resolution_clock::now();
    
    auto primes = number_theory::sieve_of_eratosthenes(N);
    
    auto end = high_resolution_clock::now();
    
    std::cout << "Sieve of Eratosthenes (up to " << N << "): " 
              << duration_cast<milliseconds>(end - start).count() << " ms\n\n";
}

void run_graph_benchmark() {
    std::cout << "--- Graph Algorithms Benchmark ---\n";
    int N = 10000; // 10k nodes
    graphs::Graph<double> g(N);
    
    // Create a chain graph with some random connections
    for(int i = 0; i < N - 1; ++i) {
        g.add_edge(i, i + 1, 1.0);
        if (i + 10 < N) g.add_edge(i, i + 10, 5.0);
    }
    
    auto start = high_resolution_clock::now();
    
    auto result = g.dijkstra(0);
    
    auto end = high_resolution_clock::now();
    
    std::cout << "Dijkstra's Algorithm (10k nodes, ~20k edges): " 
              << duration_cast<milliseconds>(end - start).count() << " ms\n\n";
}

int main() {
    std::cout << "==========================================\n";
    std::cout << "           NumC Benchmarks Suite          \n";
    std::cout << "==========================================\n\n";
    
    run_matrix_benchmark();
    run_lu_benchmark();
    run_fft_benchmark();
    run_calculus_benchmark();
    run_rk4_benchmark();
    run_kmeans_benchmark();
    run_primes_benchmark();
    run_graph_benchmark();
    
    std::cout << "Benchmarks finished.\n";
    return 0;
}
