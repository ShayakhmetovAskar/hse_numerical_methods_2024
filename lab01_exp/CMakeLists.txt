cmake_minimum_required(VERSION 3.24)
project(hse_numerical_methods_2024/lab01_exp)

set(CMAKE_CXX_STANDARD 20)

add_executable(hse_numerical_methods_2024 testexp.cpp exp.cpp ../lab02_pade_exp/testpade.cpp ../lab03_Chebyshev/testChebyshev.cpp ../lab03_Chebyshev/Chebyshev_using_built_in_funcs.cpp)

find_package(GSL REQUIRED)
target_link_libraries(hse_numerical_methods_2024 GSL::gsl GSL::gslcblas)