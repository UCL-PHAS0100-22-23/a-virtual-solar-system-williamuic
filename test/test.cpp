#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "particle.hpp"
#include <Eigen/Core>
#include <iostream>
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

// Gravitational Force
TEST_CASE("Gravitational force between two particles", "[particle]") {
    // Test case 1: Two particles with same mass and distance of 1 unit
    Particle p1(Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 1.0);
    Particle p2(Eigen::Vector3d(1.0, 0.0, 0.0), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 1.0);
    Eigen::Vector3d expected_acceleration(1.0, 0.0, 0.0);
    double epsilon = 0.0;
    Eigen::Vector3d calculated_acceleration = calcAcceleration(p1, p2, epsilon);
    REQUIRE(calculated_acceleration.isApprox(expected_acceleration, 1e-5));

    // Test case 2: Two particles with same mass and distance of 2 units
    Particle p3(Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 1.0);
    Particle p4(Eigen::Vector3d(2.0, 0.0, 0.0), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 1.0);
    expected_acceleration = Eigen::Vector3d(0.25, 0.0, 0.0);
    calculated_acceleration = calcAcceleration(p3, p4, epsilon);
    REQUIRE(calculated_acceleration.isApprox(expected_acceleration, 1e-5));
}


TEST_CASE("SumAccelerations for all particles", "[particle]") {
    // Test case 1: Single particle in empty space
    Particle p1(Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 1.0);
    std::vector<Particle> particles = {p1};
    double epsilon = 0.0;
    p1.SumAccelerations(particles, epsilon);

    REQUIRE(p1.getAcceleration().isZero());

    // Test case 2: Two particles with same mass and distance of 1 unit
    Particle p2(Eigen::Vector3d(1.0, 0.0, 0.0), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 1.0);
    particles = {p1, p2};
    p1.SumAccelerations(particles, epsilon);
    Eigen::Vector3d expected_acceleration(1.0, 0.0, 0.0);

    REQUIRE(p1.getAcceleration().isApprox(expected_acceleration, 1e-5));

    // Test case 3: Three particles in a row with equal mass and spacing
    Particle p3(Eigen::Vector3d(-1.0, 0.0, 0.0), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 1.0);
    particles = {p1, p2, p3};
    p1.SumAccelerations(particles, epsilon);
    expected_acceleration = Eigen::Vector3d(0.0, 0.0, 0.0);
    REQUIRE(p2.getAcceleration().isApprox(expected_acceleration, 1e-5));
}

TEST_CASE("Gravitational force edge cases", "[particle]") {
    // Test case 1: Two particles with same position
    Particle p1(Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 1.0);
    Particle p2(Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 1.0);
    Eigen::Vector3d expected_acceleration(0.0, 0.0, 0.0);
    double epsilon = 0.0;
    Eigen::Vector3d calculated_acceleration = calcAcceleration(p1, p2, epsilon);
    REQUIRE(calculated_acceleration.isApprox(expected_acceleration, 1e-5));

    // Test case 2: Two particles with very small distance and non-zero softening factor
    Particle p3(Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 1.0);
    Particle p4(Eigen::Vector3d(1e-10, 0.0, 0.0), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 1.0);
    expected_acceleration = Eigen::Vector3d(1e5, 0.0, 0.0);
    epsilon = 1e-5;
    calculated_acceleration = calcAcceleration(p3, p4, epsilon);
    std::cout << "Calculated acceleration: " << calculated_acceleration << std::endl; // Print out calculated acceleration

    REQUIRE(calculated_acceleration.isApprox(expected_acceleration, 1e-5));
}

TEST_CASE("SumAccelerations edge cases", "[particle]") {
    // Test case 1: Two particles with opposite mass and same position
    Particle p1(Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 1.0);
    Particle p2(Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), -1.0);
    std::vector<Particle> particles = {p1, p2};
    double epsilon = 0.0;
    p1.SumAccelerations(particles, epsilon);
    Eigen::Vector3d expected_acceleration(0.0, 0.0, 0.0);
    std::cout << "Particle acceleration: " << p1.getAcceleration() << std::endl; // Print out particle acceleration

    REQUIRE(p1.getAcceleration().isApprox(expected_acceleration, 1e-5));

    // Test case 2: Three particles in a line with very small distance and non-zero softening factor
    Particle p3(Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 1.0);
    Particle p4(Eigen::Vector3d(1e-10, 0.0, 0.0), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 1.0);
    Particle p5(Eigen::Vector3d(2e-10, 0.0, 0.0), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 1.0);
    particles = {p3, p4, p5};
    epsilon = 1e-5;
    p3.SumAccelerations(particles, epsilon);
    std::cout << "Particle acceleration: " << p3.getAcceleration() << std::endl; // Print out particle acceleration
    expected_acceleration = Eigen::Vector3d(3e5, 0.0, 0.0);
    Eigen::Vector3d expected_acceleration_p4 = Eigen::Vector3d(-1e30, 0.0, 0.0);
    Eigen::Vector3d expected_acceleration_p5 = Eigen::Vector3d(5e29, 0.0, 0.0);
    REQUIRE(p3.getAcceleration().isApprox(expected_acceleration, 1e-5));
}