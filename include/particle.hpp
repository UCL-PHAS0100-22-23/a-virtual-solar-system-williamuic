#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <Eigen/Core>
// Particle class simulates a particle in space
class Particle {
public:
  // Constructor that takes in the initial position, velocity, acceleration, and
  // mass of the particle
  Particle(const Eigen::Vector3d &in_position,
           const Eigen::Vector3d &in_velocity,
           const Eigen::Vector3d &in_acceleration, double in_mass);

  // Getter functions for the particle's position, velocity, acceleration, and
  // mass
  Eigen::Vector3d getPosition() const;
  Eigen::Vector3d getVelocity() const;
  Eigen::Vector3d getAcceleration() const;
  double getMass() const;

  // Setter function for the particle's acceleration
  void setAcceleration(const Eigen::Vector3d &in_acceleration);

  // Update the particle's position and velocity based on its current state and
  // elapsed time dt
  void update(double dt);

  // Sum the accelerations of all particles in a given vector
  void SumAccelerations(const std::vector<Particle> &particles, double epsilon);

private:
  Eigen::Vector3d position;     // The particle's position
  Eigen::Vector3d velocity;     // The particle's velocity
  Eigen::Vector3d acceleration; // The particle's acceleration
  double mass;                  // The particle's mass
};

// Calculate the acceleration between two particles using their positions and
// masses
Eigen::Vector3d calcAcceleration(const Particle &p1, const Particle &p2,
                                 double epsilon = 0.0);

// The InitialConditionGenerator abstract class defines a method for generating
// initial conditions for a simulation
class InitialConditionGenerator {
public:
  virtual std::vector<Particle> generateInitialConditions() = 0;
};

// The SolarSystemGenerator class generates initial conditions for a solar
// system simulation
class SolarSystemGenerator : public InitialConditionGenerator {
public:
  std::vector<Particle> generateInitialConditions() override;
};

// The RandomGenerator class generates initial conditions for a simulation with
// a specified number of particles
class RandomGenerator : public InitialConditionGenerator {
public:
  // Constructor that takes in the number of particles to generate
  RandomGenerator(int num_particles) : num_particles_(num_particles) {}

  // Generate initial conditions for the simulation
  std::vector<Particle> generateInitialConditions() override;

private:
  int num_particles_; // The number of particles to generate
};
#endif