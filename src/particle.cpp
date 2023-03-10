#include "particle.hpp"
#include <Eigen/Core>

Particle::Particle(double in_mass) :
    mass{in_mass}
    {}

double Particle::getMass() const {
    return mass;
}