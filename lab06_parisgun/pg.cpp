#include "atm.hpp"
#include "atm.cpp"
#include "consts.hpp"
#include "pg.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <cmath>

namespace ADAAI {
    double ParisGun::interpolate(double v_M){
        int idx = -1;
        if (v_M < speed[0]){
            return c_d[0];
        }
        else if (v_M > speed[12]){
            return c_d[12];
        }
        else {
            for (int i = 0; i < n - 1; ++i) {
                if (v_M >= speed[i] && v_M <= speed[i + 1]) {
                    idx = i;
                    break;
                }
            }
        }

        if (idx != -1) {
            double c_d_interp = c_d[idx] + (c_d[idx + 1] - c_d[idx]) * (v_M - speed[idx]) / (speed[idx + 1] - speed[idx]);
//            std::cout << "c_d_interp=" << c_d_interp << '\n';
            return c_d_interp;
        } else {
            std::cerr << "Speed out of range." << std::endl;
            return -1;
        }
    }

    int ParisGun::projectile_equations(double t, const double u[], double dudt[], void* params) {
        (void)(t);
        double Q;
        double x = u[0];
        double v_x = u[1];
        double y = u[2];
        double v_y = u[3];
        double v = std::sqrt(v_x * v_x + v_y * v_y);
        double altitude = y;

        double density = ADAAI::SAM::get_density(altitude);
        double A = ADAAI::SAM::get_A(altitude);
        double M = v / A;
//        std::cout << "M=" << M << '\n';
        double c_d_M = ParisGun::interpolate(M);
//        std::cout << "c_d_M=" << c_d_M << '\n';
        Q = density * v * v * S * c_d_M / 2;
//        std::cout << "Q=" << p->Q << '\n';

        dudt[0] = u[1]; // dx/dt = v_x
        dudt[1] = - Q * u[1] / v / mass; // dv_x/dt
        dudt[2] = u[3]; // dy/dt = v_y
        dudt[3] = - Q * u[3] / v / mass - g; // dv_y/dt
//        p->Q = Q;
        return GSL_SUCCESS;
    }

    void ParisGun::simulate_projectile_motion() {
        const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rkf45;
        gsl_odeiv2_system sys = {ParisGun::projectile_equations, nullptr, 4, nullptr};
        gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, T, 1e-6, 1e-6, 0.0);

        double t = 0.0;
        double dt = 5.0;

        double u[4];
        u[0] = 0.0; // x
        u[1] = v_0 * std::cos(alpha_rad); // v_x
        u[2] = 0.0; // y
        u[3] = v_0 * std::sin(alpha_rad); // v_y


        while (u[2] >= 0) {
            std::cout << "t=" << t << " x=" << u[0] << " y=" << u[2] << " vx=" << u[1] << " vy=" << u[3] << std::endl;
            int status = gsl_odeiv2_driver_apply(d, &t, t+dt, u);

            if (status != GSL_SUCCESS) {
                std::cerr << "Error: integration failed" << std::endl;
                break;
            }
        }
        gsl_odeiv2_driver_free(d);
    }
}

int main(){
    ADAAI::ParisGun::simulate_projectile_motion();
//    std::cout << ADAAI::v_0 * std::sin(ADAAI::alpha_rad);
    return 0;
}