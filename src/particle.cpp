#include "particle.hpp"
#include <Eigen/Core>
#include <random>
Particle::Particle(const Eigen::Vector3d &in_position,
                   const Eigen::Vector3d &in_velocity,
                   const Eigen::Vector3d &in_acceleration, double in_mass)
    : position(in_position), velocity(in_velocity),
      acceleration(in_acceleration), mass(in_mass) {}

Eigen::Vector3d Particle::getPosition() const { return position; }

Eigen::Vector3d Particle::getVelocity() const { return velocity; }

Eigen::Vector3d Particle::getAcceleration() const { return acceleration; }

double Particle::getMass() const { return mass; }

void Particle::update(double dt) {
  position += dt * velocity;
  velocity += dt * acceleration;
}

void Particle::setAcceleration(const Eigen::Vector3d &in_acceleration) {
  acceleration = in_acceleration;
}

Eigen::Vector3d calcAcceleration(const Particle &p1, const Particle &p2,
                                 double epsilon) {
  Eigen::Vector3d r = p2.getPosition() - p1.getPosition();
  double distanceSquared = r.squaredNorm();
  double softeningFactorSquared = epsilon * epsilon;
  // Check when the distance is too small between too particle
  if (distanceSquared < 1e-20) {
    return Eigen::Vector3d(0, 0, 0);
  }
  double denominator = pow(distanceSquared + softeningFactorSquared, 1.5);
  double a = p2.getMass() / denominator;
  return a * r;
}

void Particle::SumAccelerations(const std::vector<Particle> &particles,
                                double epsilon) {
  acceleration.setZero();
  for (const auto &p : particles) {
    // Sum the accelerations without itself
    if (&p != this) {
      acceleration += calcAcceleration(*this, p, epsilon);
    }
  }
}

std::vector<Particle> SolarSystemGenerator::generateInitialConditions() {
  std::vector<double> masses = {
      1.0,           1.0 / 6023600, 1.0 / 408524, 1.0 / 332946.038,
      1.0 / 3098710, 1.0 / 1047.55, 1.0 / 3499,   1.0 / 22962,
      1.0 / 19352};
  std::vector<double> distances = {0.0, 0.4, 0.7,  1.0, 1.5,
                                   5.2, 9.5, 19.2, 30.1};
  std::vector<Particle> particles;
  particles.reserve(masses.size());
  std::random_device rd;
  // Developer can change 80, it is a seed control the random generation
  std::mt19937 gen(80);
  std::uniform_real_distribution<> dist(0, 2 * M_PI);
  Eigen::Vector3d position(0.0, 0.0, 0.0);
  Eigen::Vector3d velocity(0.0, 0.0, 0.0);
  Eigen::Vector3d acceleration(0.0, 0.0, 0.0);
  // Put the central particle into the first place of the list
  particles.emplace_back(position, velocity, acceleration, masses[0]);
  for (size_t i = 1; i < masses.size(); ++i) {
    double theta = dist(gen);
    double r = distances[i];
    // Randomly generate angles and put them into list
    Eigen::Vector3d position(r * sin(theta), r * cos(theta), 0.0);
    Eigen::Vector3d velocity(-sqrt(1.0 / r) * cos(theta),
                             sqrt(1.0 / r) * sin(theta), 0.0);
    Eigen::Vector3d acceleration(0.0, 0.0, 0.0);
    particles.emplace_back(position, velocity, acceleration, masses[i]);
  }
  return particles;
}

std::vector<Particle> RandomGenerator::generateInitialConditions() {
  std::random_device rd;
  // Developer can change 80, it is a seed control the random generation
  std::mt19937 gen(80);
  std::uniform_real_distribution<> mass_dist(1.0 / 6000000, 1.0 / 1000);
  std::uniform_real_distribution<> distance_dist(0.4, 30);
  std::uniform_real_distribution<> angle_dist(0, 2 * M_PI);
  std::vector<Particle> particles;
  particles.reserve(num_particles_);
  Eigen::Vector3d position(0.0, 0.0, 0.0);
  Eigen::Vector3d velocity(0.0, 0.0, 0.0);
  Eigen::Vector3d acceleration(0.0, 0.0, 0.0);
  // Put the central particle into the first place of the list
  particles.emplace_back(position, velocity, acceleration, 1.0);
  for (size_t i = 1; i < num_particles_; ++i) {
    double mass = mass_dist(gen);
    double distance = distance_dist(gen);
    double angle = angle_dist(gen);
    // Randomly generate angles and put them into list
    Eigen::Vector3d position(distance * sin(angle), distance * cos(angle), 0.0);
    Eigen::Vector3d velocity(-sqrt(1.0 / distance) * cos(angle),
                             sqrt(1.0 / distance) * sin(angle), 0.0);
    Eigen::Vector3d acceleration(0.0, 0.0, 0.0);
    particles.emplace_back(position, velocity, acceleration, mass);
  }
  return particles;
}