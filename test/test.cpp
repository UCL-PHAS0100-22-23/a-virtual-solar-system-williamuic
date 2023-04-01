#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "particle.hpp"
#include <Eigen/Core>

using Catch::Matchers::WithinRel;

TEST_CASE("Particle moves with no acceleration", "[particle]") {
    Particle p(Eigen::Vector3d(1.0, 2.0, 3.0), Eigen::Vector3d(4.0, 5.0, 6.0), Eigen::Vector3d(0.0, 0.0, 0.0), 1.0);
    double dt = 0.1;
    p.update(dt);
    REQUIRE(p.getPosition().isApprox(Eigen::Vector3d(1.4, 2.5, 3.6)));
    REQUIRE(p.getVelocity().isApprox(Eigen::Vector3d(4.0, 5.0, 6.0)));
}

TEST_CASE("Particle moves with constant acceleration", "[particle]") {
    Particle p(Eigen::Vector3d(1.0, 2.0, 3.0), Eigen::Vector3d(4.0, 5.0, 6.0), Eigen::Vector3d(0.0, 0.0, -9.8), 1.0);
    double dt = 0.1;
    p.update(dt);
    REQUIRE(p.getPosition().isApprox(Eigen::Vector3d(1.4, 2.5, 3.6)));
    REQUIRE(p.getVelocity().isApprox(Eigen::Vector3d(4.0, 5.0, 5.02)));
}

TEST_CASE("Particle orbits around origin with artificial acceleration", "[particle]") {
    Particle p(Eigen::Vector3d(1.0, 0.0, 0.0), Eigen::Vector3d(0.0, 1.0, 0.0), Eigen::Vector3d(-1.0, 0.0, 0.0), 1.0);
    double dt = 0.001;
    int t_final = 2.0 * M_PI / dt;
    for (int t = 0.0; t < t_final; t++) {
        Eigen::Vector3d x = p.getPosition();
        Eigen::Vector3d a = -x.normalized();
        p.setAcceleration(a);
        p.update(dt);
    }
    REQUIRE(p.getPosition().isApprox(Eigen::Vector3d(1.0, 0.0, 0.0),1e-2));
    REQUIRE(p.getVelocity().isApprox(Eigen::Vector3d(0.0, 1.0, 0.0),1e-2));
}

TEST_CASE("Particle has correct mass", "[particle]") {
    Particle p(Eigen::Vector3d(1.0, 2.0, 3.0), Eigen::Vector3d(4.0, 5.0, 6.0), Eigen::Vector3d(0.0, 0.0, -9.8), 2.0);
    REQUIRE(p.getMass() == 2.0);
}