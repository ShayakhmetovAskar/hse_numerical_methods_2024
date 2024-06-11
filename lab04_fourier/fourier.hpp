#ifndef ADAAI_FOURIER_HPP
#define ADAAI_FOURIER_HPP

double f(double x);
double f_sin(double x, void *p);
double f_cos(double x, void *p);
double calculateFourierCoefficientA(int n);
double calculateFourierCoefficientB(int n);

#endif // ADAAI_FOURIER_HPP
