#include "particle.hpp"
#include <Eigen/Core>

Particle::Particle(const Eigen::Vector3d& in_position, const Eigen::Vector3d& in_velocity, const Eigen::Vector3d& in_acceleration, double in_mass)
    : position(in_position), velocity(in_velocity), acceleration(in_acceleration), mass(in_mass) {}

Eigen::Vector3d Particle::getPosition() const {
    return position;
}

Eigen::Vector3d Particle::getVelocity() const {
    return velocity;
}

Eigen::Vector3d Particle::getAcceleration() const {
    return acceleration;
}

double Particle::getMass() const {
    return mass;
}

void Particle::update(double dt) {
    position += dt * velocity;
    velocity += dt * acceleration;
}


void Particle::setAcceleration(const Eigen::Vector3d& in_acceleration) {
    acceleration = in_acceleration;
}