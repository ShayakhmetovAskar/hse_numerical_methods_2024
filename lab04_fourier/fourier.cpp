#include "fourier.hpp"
#include <cmath>
#include <gsl/gsl_integration.h>
#include <iostream>

// Структура для передачи параметров в функции f_sin, f_cos
struct params {
  int n;
};

double f(double x) { return std::exp(x); }

double f_sin(double x, void *p) {
  struct params *par = (struct params *)p;
  int n = par->n;

  return f(x) * std::sin(n * x);
}

double f_cos(double x, void *p) {
  struct params *par = (struct params *)p;
  int n = par->n;

  return f(x) * std::cos(n * x);
}

// to calculate a_n
double calculateFourierCoefficientA(int n) {
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  F.function =
      &f_cos; // if n = 0 then cos(nx)=1 so f_cos=f - that's what we need
  struct params p = {n};
  F.params = &p;
  double result, error;
  gsl_integration_qags(&F, -M_PI, M_PI, 0, 1e-7, n, w, &result, &error);
  gsl_integration_workspace_free(w);
  return result / M_PI;
}

// to calculate b_n
double calculateFourierCoefficientB(int n) {
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  F.function = &f_sin;
  struct params p = {n};
  F.params = &p;
  double result, error;
  gsl_integration_qags(&F, -M_PI, M_PI, 0, 1e-7, n, w, &result, &error);
  gsl_integration_workspace_free(w);
  return result / M_PI;
}
