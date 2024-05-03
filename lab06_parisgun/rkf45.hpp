#ifndef ADAAI_RKF45_HPP
#define ADAAI_RKF45_HPP
#include <iostream>
#include <vector>

namespace ADAAI {
    class Timestepper { // Class representing the timestepper
    public:
        virtual std::vector<double> step(double t, const std::vector<double> &y, double &dt) = 0;
    };

    class HarmonicOscillator { // Class representing the harmonic oscillator (projectile motion)
    public:
        std::vector<double> f(double t, const std::vector<double> &u);
    };

    class RKF45 : public Timestepper { // Class for Runge-Kutta-Fehlberg 4(5) method
    public:
        RKF45(HarmonicOscillator* system_ptr);

        std::vector<double> step(double t, const std::vector<double> &y, double &dt) override;
        // RKF45 coeffs (Butcher)
        static constexpr double a[6][6] = {
                {0,             0,              0,              0,             0,          0},
                {1.0 / 4,       0,              0,              0,             0,          0},
                {3.0 / 32,      9.0 / 32,       0,              0,             0,          0},
                {1932.0 / 2197, -7200.0 / 2197, 7296.0 / 2197,  0,             0,          0},
                {439.0 / 216,   -8,             3680.0 / 513,   -845.0 / 4104, 0,          0},
                {-8.0 / 27,     2,              -3544.0 / 2565, 1859.0 / 4104, -11.0 / 40, 0}
        };
        static constexpr double b[6] = {16.0 / 135, 0, 6656.0 / 12825, 28561.0 / 56430, -9.0 / 50, 2.0 / 55};
        static constexpr double b_star[6] = {25.0 / 216, 0, 1408.0 / 2565, 2197.0 / 4104, -1.0 / 5, 0};
        static constexpr double c[6] = {0, 1.0 / 4, 3.0 / 8, 12.0 / 13, 1, 1.0 / 2};

    private:
        HarmonicOscillator *system;
    };


    class Observer { // Observer class to track and record data during simulation
    public:
        virtual void observe(double t, const std::vector<double> &y) = 0;
    };

    class MaxHeightRangeObserver : public Observer { // Example observer to print maximum height and range
    public:
        void observe(double t, const std::vector<double> &y) override;

        double getMaxHeight();
        double getMaxRange();

    private:
        double max_height = 0;
        double max_range = 0;
    };

    class Integrator { // Integrator class to manage the simulation
    public:
        Integrator(Timestepper *timestepper, HarmonicOscillator *system, Observer *observer);

        void integrate(double t_start, const std::vector<double> &y_start, double t_end, double dt);

    private:
        Timestepper *timestepper;
        HarmonicOscillator *system;
        Observer *observer;
    };


}
#endif //ADAAI_RKF45_HPP
