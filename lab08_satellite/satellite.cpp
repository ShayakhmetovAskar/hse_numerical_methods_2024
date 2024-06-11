#include "consts.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

const double mu = 398600.4418; // km^3/s^2
const double J2 = 1.0827e-3;
const double Re = 6378.137; // km, радиус Земли

class EverhartIntegrator {
public:
  EverhartIntegrator(double step) : dt(step) { calculateCoefficients(); }

  void integrate(std::vector<double> &position, std::vector<double> &velocity,
                 int steps, const std::string &filename) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
    }

    outfile
        << "Step,PositionX,PositionY,PositionZ,VelocityX,VelocityY,VelocityZ\n";
    for (int step = 1; step < steps; ++step) {
      everhartInterpolation(position, velocity);
      if (step % 86400 == 0) { // Печать результатов каждый день
        outfile << step << "," << position[0] << "," << position[1] << ","
                << position[2] << "," << velocity[0] << "," << velocity[1]
                << "," << velocity[2] << "\n";
        std::cout << "Step: " << step << " Position: (" << position[0] << ", "
                  << position[1] << ", " << position[2] << ")"
                  << " Velocity: (" << velocity[0] << ", " << velocity[1]
                  << ", " << velocity[2] << ")\n";
      }
    }
    outfile.close();
    std::cout << "Data written to " << filename << std::endl;
  }

private:
  double dt;
  std::vector<double> alpha;
  std::vector<double> coeff_position;
  std::vector<double> coeff_velocity;

  void calculateCoefficients() {
    // Коэффициенты alpha для 6-го порядка
    alpha = {1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125};

    // Коэффициенты интерполяции для позиции и скорости
    coeff_position = {1.0 / 720.0, 1.0 / 120.0, 1.0 / 24.0,
                      1.0 / 6.0,   1.0 / 2.0,   1.0};
    coeff_velocity = {1.0 / 144.0, 1.0 / 24.0, 1.0 / 6.0,
                      1.0 / 2.0,   5.0 / 6.0,  1.0};
  }

  std::vector<double> computeAcceleration(const std::vector<double> &position) {
    double r2 = position[0] * position[0] + position[1] * position[1] +
                position[2] * position[2];
    double r = std::sqrt(r2);
    double r3 = r * r * r;
    double r5 = r2 * r3;

    double z2 = position[2] * position[2];
    double factor = ADAAI::J2 * ADAAI::mu * ADAAI::Re * ADAAI::Re / r5;

    std::vector<double> acceleration(3);
    //        acceleration[0] = -mu * position[0] * r * r / r3;
    //        acceleration[0] *= 1 - J2 * mu * Re * Re * 3 * z2 / 2 * ((-5 *
    //        position[0]) / pow(r, 7.0 / 2));  // dU/dx acceleration[0] -= mu *
    //        Re * Re / 2 * ((-3 * position[0]) / pow(r, 5.0 / 2));  // dU/dx
    //        acceleration[1] = -mu * position[1] * r * r / r3;
    //        acceleration[1] *= 1 - J2 * mu * Re * Re * 3 * z2 / 2 * ((-5 *
    //        position[1]) / pow(r, 7.0 / 2)); acceleration[1] -= mu * Re * Re /
    //        2 * ((-3 * position[1]) / pow(r, 5.0 / 2)); acceleration[2] = mu *
    //        position[2] * r * r / r3;  // dU/dz acceleration[2] *= 1 - J2 * mu
    //        * Re * Re * 3 / 2 *
    //                           ((2 * position[0] * position[0] * position[2] +
    //                           2 * position[1] * position[1] * position[2] -
    //                             3 * position[2] * position[2] * position[2])
    //                             / pow(r, 7.0 / 2));
    //        acceleration[2] -= mu * Re * Re / 2 * ((-3 * position[2]) /
    //        pow(r, 5.0 / 2));
    //
    acceleration[0] =
        -ADAAI::mu * position[0] / r3 * (1 - factor * (5 * z2 / r2 - 1));
    acceleration[1] =
        -ADAAI::mu * position[1] / r3 * (1 - factor * (5 * z2 / r2 - 1));
    acceleration[2] =
        -ADAAI::mu * position[2] / r3 * (1 - factor * (5 * z2 / r2 - 3));

    return acceleration;
  }

  void everhartInterpolation(std::vector<double> &position,
                             std::vector<double> &velocity) {
    const int m = alpha.size();
    std::vector<std::vector<double>> positions(m, std::vector<double>(3));
    std::vector<std::vector<double>> velocities(m, std::vector<double>(3));
    std::vector<std::vector<double>> accelerations(m, std::vector<double>(3));

    positions[0] = position;
    velocities[0] = velocity;
    accelerations[0] = computeAcceleration(position);

    for (int i = 1; i < m; ++i) {
      for (int j = 0; j < 3; ++j) {
        positions[i][j] =
            positions[0][j] + alpha[i] * dt * velocities[i - 1][j] +
            0.5 * alpha[i] * alpha[i] * dt * dt * accelerations[i - 1][j];
      }
      accelerations[i] = computeAcceleration(positions[i]);
      for (int j = 0; j < 3; ++j) {
        velocities[i][j] =
            velocities[0][j] + alpha[i] * dt * accelerations[i - 1][j];
      }
    }

    for (int i = 0; i < 3; ++i) {
      position[i] += dt * (coeff_position[0] * velocities[0][i] +
                           coeff_position[1] * velocities[1][i] +
                           coeff_position[2] * velocities[2][i] +
                           coeff_position[3] * velocities[3][i] +
                           coeff_position[4] * velocities[4][i] +
                           coeff_position[5] * velocities[5][i]);
      velocity[i] += dt * (coeff_velocity[0] * accelerations[0][i] +
                           coeff_velocity[1] * accelerations[1][i] +
                           coeff_velocity[2] * accelerations[2][i] +
                           coeff_velocity[3] * accelerations[3][i] +
                           coeff_velocity[4] * accelerations[4][i] +
                           coeff_velocity[5] * accelerations[5][i]);
    }
  }
};

int main() {
  // Начальные условия для круговой орбиты на высоте 7500 км
  std::vector<double> position = {0.0, 0.0, 7500.0}; // km
  std::vector<double> velocity = {
      7.12, 0.0, 0.0}; // km/s, начальная скорость в направлении оси X

  double dt = 1.0;             // Шаг времени (секунды)
  int steps = 365 * 24 * 3600; // Количество шагов за один год

  EverhartIntegrator integrator(dt);
  integrator.integrate(position, velocity, steps, "satellite_orbit.csv");

  return 0;
}
