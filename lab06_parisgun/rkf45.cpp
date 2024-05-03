#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include "consts.hpp"
#include "atm.cpp"
#include "rkf45.hpp"
#include "pg.hpp"

namespace ADAAI {
    std::vector<double> Timestepper::step(double t, const std::vector<double> &y, double &dt) {
        return {};
    }


    std::vector<double> HarmonicOscillator::f(double t, const std::vector<double> &u) {
        double x = u[0];
        double y = u[1];
        double v_x = u[2];
        double v_y = u[3];
        double v = std::sqrt(v_x * v_x + v_y * v_y);

        double density = ADAAI::SAM::get_density(y);
        double A = ADAAI::SAM::get_A(y);
        double M = v / A;
        double c_d_M = ADAAI::ParisGun::interpolate(M);
        double Q = density * v * v * ADAAI::S * c_d_M / 2;
        double ax = Q * v_x / v / ADAAI::mass;
        double ay = Q * v_y / v / ADAAI::mass - ADAAI::g;

        double a_coriolis_x = 2 * ADAAI::omega * v_y;
        double a_coriolis_y = -2 * ADAAI::omega * v_x;

        // Centrifugal acceleration
        double a_centrifugal = v * v / (R + y);

        return {v_x, v_y, ax + a_coriolis_x - a_centrifugal * sin(alpha_rad),
                ay + a_coriolis_y - a_centrifugal * cos(alpha_rad), v_x / (R + y) / cos(alpha_rad)};
    }

    RKF45::RKF45(HarmonicOscillator* system_ptr) : system(system_ptr) {}
    std::vector<double> RKF45::step(double t, const std::vector<double> &y, double &dt) {

        // Stage calculations
        std::vector<std::vector<double> > k(6, std::vector<double>(y.size()));
        for (int i = 0; i < 6; ++i) {
            std::vector<double> y_temp(y.size());
            for (size_t j = 0; j < y.size(); ++j) {
                double sum = 0;
                for (int l = 0; l < i; ++l) {
                    sum += a[i][l] * k[l][j];
                }
                y_temp[j] = y[j] + dt * sum;
            }
            k[i] = system->f(t + c[i] * dt, y_temp);
        }

        // 4th and 5th order approximations
        std::vector<double> y_next(y.size());
        std::vector<double> y_next_star(y.size());
        for (size_t i = 0; i < y.size(); ++i) {
            double sum = 0;
            double sum_star = 0;
            for (int j = 0; j < 6; ++j) {
                sum += b[j] * k[j][i];
                sum_star += b_star[j] * k[j][i];
            }
            y_next[i] = y[i] + dt * sum;
            y_next_star[i] = y[i] + dt * sum_star;
        }

        // Error estimation
        double error = 0;
        for (size_t i = 0; i < y.size(); ++i) {
            error += (y_next[i] - y_next_star[i]) * (y_next[i] - y_next_star[i]);
        }
        error = sqrt(error / y.size());

        // Adaptive step size control
        double tol = 1e-6; // Tolerance
        double s = 0.9 * std::pow(tol / error, 1.0 / 5);
        if (s < 0.75) dt *= 0.75;
        else if (s > 1.5) dt *= 1.2;

        // Ensure dt doesn't get too small
        if (dt < std::numeric_limits<double>::epsilon()) {
            std::cerr << "Warning: Time step became too small, integration might be inaccurate." << std::endl;
        }

        return y_next;
    }


    void Observer::observe(double t, const std::vector<double> &y) {}


    void MaxHeightRangeObserver::observe(double t, const std::vector<double> &y) {
        if (y[1] > max_height) {
            max_height = y[1];
        }
        if (y[0] > max_range && y[1] <= 0) {
            max_range = y[0];
        }
    }

    double MaxHeightRangeObserver::getMaxHeight() {
        return max_height;
    }

    double MaxHeightRangeObserver::getMaxRange() {
        return max_range;
    }


    Integrator::Integrator(Timestepper *timestepper, HarmonicOscillator *system, Observer *observer) :
            timestepper(timestepper), system(system), observer(observer) {}

    void Integrator::integrate(double t_start, const std::vector<double> &y_start, double t_end, double dt) {
        double t = t_start;
        std::vector<double> y = y_start;
        while (t < t_end) {
            y = timestepper->step(t, y, dt);
            observer->observe(t, y);
            t += dt;
        }
    }

}

int main() {
    std::vector<double> y0 = {0, 0, ADAAI::v_0 * cos(ADAAI::alpha_rad), ADAAI::v_0 * sin(ADAAI::alpha_rad)};
    ADAAI::HarmonicOscillator system;
    ADAAI::RKF45 timestepper(&system);
    ADAAI::MaxHeightRangeObserver observer;
    ADAAI::Integrator integrator(&timestepper, &system, &observer);

    integrator.integrate(0, y0, 200, 0.05);

    std::cout << "Max height: " << observer.getMaxHeight() << " m" << std::endl;
    std::cout << "Max range: " << observer.getMaxRange() << " m" << std::endl;

    return 0;
}