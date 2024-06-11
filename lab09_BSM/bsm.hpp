#ifndef BSM_HPP
#define BSM_HPP

#include <vector>

class BSM {
public:
    BSM(double S0, double K, double r, double T, double sigma);
    double optionPrice(int M, int N);

private:
    double S0;
    double K;
    double r;
    double T;
    double sigma;

    double payoff(double S);
    void setupGrid(std::vector<double>& prices, std::vector<double>& values, int N);
    void stepBackwards(std::vector<double>& values, const std::vector<double>& prices, int M, int N);
};

#endif // BSM_HPP

