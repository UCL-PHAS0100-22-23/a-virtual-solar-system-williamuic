#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <Eigen/Core>

class Particle {
public:
    Particle(const Eigen::Vector3d& in_position, const Eigen::Vector3d& in_velocity, const Eigen::Vector3d& in_acceleration, double in_mass);

    Eigen::Vector3d getPosition() const;
    Eigen::Vector3d getVelocity() const;
    Eigen::Vector3d getAcceleration() const;
    double getMass() const;

    void setAcceleration(const Eigen::Vector3d& in_acceleration);
    void update(double dt);
    void SumAccelerations(const std::vector<Particle>& particles, double epsilon);
 

private:
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    Eigen::Vector3d acceleration;
    double mass;
};
Eigen::Vector3d calcAcceleration(const Particle& p1, const Particle& p2, double epsilon=0.0);
class InitialConditionGenerator {
public:
    virtual std::vector<Particle> generateInitialConditions() = 0;
};

class SolarSystemGenerator : public InitialConditionGenerator {
public:
    std::vector<Particle> generateInitialConditions() override;
};

class RandomGenerator : public InitialConditionGenerator {
public:
    RandomGenerator(int num_particles) : num_particles_(num_particles) {}
    std::vector<Particle> generateInitialConditions() override;

private:
    int num_particles_;
};
#endif