cmake_minimum_required(VERSION 3.24)
project(hse_numerical_methods_2024/lab03_Chebyshev)

set(CMAKE_CXX_STANDARD 20)

add_executable(hse_numerical_methods_2024 Chebyshev.cpp testChebyshev.cpp)

find_package(GSL REQUIRED)
target_link_libraries(hse_numerical_methods_2024 GSL::gsl GSL::gslcblas)