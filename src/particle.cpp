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

Eigen::Vector3d calcAcceleration(const Particle& p1, const Particle& p2, double epsilon) {
    Eigen::Vector3d r = p2.getPosition() - p1.getPosition();
    double distanceSquared = r.squaredNorm();
    double softeningFactorSquared = epsilon * epsilon;

    if (distanceSquared < 1e-20) {
        return Eigen::Vector3d(0, 0, 0);
    }

    double denominator = pow(distanceSquared + softeningFactorSquared, 1.5);
    double a = p2.getMass() / denominator;
    return a * r;
}

void Particle::SumAccelerations(const std::vector<Particle>& particles, double epsilon) {
    acceleration.setZero();
    for (const auto& p : particles) {
        if (&p != this) {
            acceleration += calcAcceleration(*this, p, epsilon);
        }
    }
}