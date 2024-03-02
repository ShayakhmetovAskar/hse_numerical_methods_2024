#include "fourier.hpp"
#include <iostream>
#include <cmath>
#include <gsl/gsl_integration.h>

// Структура для передачи параметров в функцию f
struct params {
    int n;
};

// Функция для вычисления функции f(x) = e^x * cos(n*x) или e^x * sin(n*x)
double f(double x, void* params) {
    return std::exp(x);
}

double f_sin(double x, void* p) {
    struct params* par = (struct params*)p;
    int n = par->n;

    return std::exp(x) * std::sin(n * x);
}

double f_cos(double x, void* p) {
    struct params* par = (struct params*)p;
    int n = par->n;

    return std::exp(x) * std::cos(n * x);
}

// Функция для вычисления коэффициента a_n ряда Фурье
double calculateFourierCoefficientA(int n) {
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
    gsl_function F;
    F.function = &f_cos; // if n = 0 then cos(nx)=1 so f_cos=f - that's what we need
    struct params p = {n};
    F.params = &p;
    double result, error;
    gsl_integration_qags(&F, -M_PI, M_PI, 0, 1e-7, n, w, &result, &error);
    gsl_integration_workspace_free(w);
    return result / M_PI;
}

// Функция для вычисления коэффициента b_n ряда Фурье
double calculateFourierCoefficientB(int n) {
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
    gsl_function F;
    F.function = &f_sin;
    struct params p = {n};
    F.params = &p;
    double result, error;
    gsl_integration_qags(&F, -M_PI, M_PI, 0, 1e-7, n, w, &result, &error);
    gsl_integration_workspace_free(w);
    return result / M_PI;
}

