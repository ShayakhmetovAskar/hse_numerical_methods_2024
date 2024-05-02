cmake_minimum_required(VERSION 3.24)
project(hse_numerical_methods_2024/lab06_parisgun)

set(CMAKE_CXX_STANDARD 20)

add_executable(hse_numerical_methods_2024 pg.cpp rkf45.cpp atm.cpp)

find_package(GSL REQUIRED)
target_link_libraries(hse_numerical_methods_2024 GSL::gsl GSL::gslcblas)