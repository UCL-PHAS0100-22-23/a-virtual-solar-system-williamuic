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
#endif