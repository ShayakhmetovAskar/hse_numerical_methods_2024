#include "BSM.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

BSM::BSM(double S0, double K, double r, double T, double sigma)
        : S0(S0), K(K), r(r), T(T), sigma(sigma) {}

double BSM::payoff(double S) {
    return std::max(S - K, 0.0);
}

void BSM::setupGrid(std::vector<double> &prices, std::vector<double> &values, int N) {
    double dS = 2.0 * S0 / N;
    for (int i = 0; i <= N; ++i) {
        prices[i] = i * dS;
        values[i] = payoff(prices[i]);
    }
}

void BSM::stepBackwards(std::vector<double> &values, const std::vector<double> &prices, int M, int N) {
    double dt = T / M;
    double dS = 2.0 * S0 / N;

    for (int j = M - 1; j >= 0; --j) {
        for (int i = 1; i < N; ++i) {
            double delta = (values[i + 1] - values[i - 1]) / (2.0 * dS);
            double gamma = (values[i + 1] - 2.0 * values[i] + values[i - 1]) / (dS * dS);
            double theta = -0.5 * sigma * sigma * prices[i] * prices[i] * gamma - r * prices[i] * delta + r * values[i];

            values[i] += theta * dt;
        }
        values[0] = 0;
        values[N] = prices[N] - K * exp(-r * (T - j * dt));
    }
}

double BSM::optionPrice(int M, int N) {
    std::vector<double> prices(N + 1);
    std::vector<double> values(N + 1);

    setupGrid(prices, values, N);
    stepBackwards(values, prices, M, N);

    int idx = static_cast<int>(S0 / (2.0 * S0 / N));
    return values[idx];
}

int main() {
    double S0 = 100.0;
    double K = 100.0;
    double r = 0.05;
    double T = 1.0;
    double sigma = 0.2;

    int M = 100;
    int N = 100;

    BSM bsm(S0, K, r, T, sigma);
    double price = bsm.optionPrice(M, N);

    std::cout << "Price: " << price << std::endl;

    return 0;
}