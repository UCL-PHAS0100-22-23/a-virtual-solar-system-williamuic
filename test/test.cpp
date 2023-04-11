#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "particle.hpp"
#include <Eigen/Core>
#include <iostream>
#include <random>
using Catch::Matchers::WithinRel;
std::vector<Particle> create_initial_particles() {
    std::vector<double> masses = {1.0, 1.0/6023600, 1.0/408524, 1.0/332946.038, 1.0/3098710, 1.0/1047.55, 1.0/3499, 1.0/22962, 1.0/19352};
    std::vector<double> distances = {0.0, 0.4, 0.7, 1.0, 1.5, 5.2, 9.5, 19.2, 30.1};
    std::vector<Particle> particles;
    particles.reserve(masses.size());
    //Developer could choose using rd or use a fix number so that every time it will generate same random number
    std::random_device rd;
    std::mt19937 gen(80);
    std::uniform_real_distribution<> dist(0, 2 * M_PI);
    // Create particle for the Sun
    Eigen::Vector3d position(0.0, 0.0, 0.0);
    Eigen::Vector3d velocity(0.0, 0.0, 0.0);
    Eigen::Vector3d acceleration(0.0, 0.0, 0.0);
    particles.emplace_back(position, velocity, acceleration, masses[0]);
    // Create particles for each planet
    for (size_t i = 1; i < masses.size(); ++i) {
        double theta = dist(gen);
        double r = distances[i];

        Eigen::Vector3d position(r * sin(theta), r * cos(theta), 0.0);
        Eigen::Vector3d velocity(-sqrt(1.0 / r) * cos(theta), sqrt(1.0 / r) * sin(theta), 0.0);
        Eigen::Vector3d acceleration(0.0, 0.0, 0.0);
        particles.emplace_back(position, velocity, acceleration, masses[i]);
    }
    return particles;
}
double calcKineticEnergy(const Particle& p) {
    double v_squared = p.getVelocity().squaredNorm();
    return 0.5 * p.getMass() * v_squared;
}

double calcPotentialEnergy(const Particle& p1, const Particle& p2) {
    Eigen::Vector3d r = p2.getPosition() - p1.getPosition();
    double distance = r.norm();

    // Avoid division by zero
    if (distance < 1e-10 || &p1 == &p2) {
        return 0.0;
    }

    double energy = -0.5 * p1.getMass() * p2.getMass() / distance;
    return energy;
}

double calcTotalEnergy(const std::vector<Particle>& particles) {
    double kinetic_energy = 0.0;
    double potential_energy = 0.0;

    // Calculate kinetic and potential energy for each particle
    for (size_t i = 0; i < particles.size(); ++i) {
        const Particle& p1 = particles[i];
        kinetic_energy += calcKineticEnergy(p1);

        for (size_t j = i + 1; j < particles.size(); ++j) {
            const Particle& p2 = particles[j];
            potential_energy += calcPotentialEnergy(p1, p2);
        }
    }
    return kinetic_energy + potential_energy;
}
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

//1.3 Simulator
TEST_CASE("create_initial_particles with correct number of particles", "[particle]") {
    std::vector<Particle> particles = create_initial_particles();
    REQUIRE(particles.size() == 9);
}

TEST_CASE("create_initial_particles with correct masses", "[particle]") {
    std::vector<Particle> particles = create_initial_particles();
    std::vector<double> expected_masses = {1.0, 1.0/6023600, 1.0/408524, 1.0/332946.038, 1.0/3098710, 1.0/1047.55, 1.0/3499, 1.0/22962, 1.0/19352};
    for (size_t i = 0; i < expected_masses.size(); ++i) {
        REQUIRE(particles[i].getMass() == expected_masses[i]);
    }
}

TEST_CASE("single body orbits around the Sun", "[solarSystem]") {
    // Set up a single particle orbiting around the Sun
    Particle p1(Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 1.0);
    Particle p2(Eigen::Vector3d(1.0, 0.0, 0.0), Eigen::Vector3d(0.0, 1.0, 0.0), Eigen::Vector3d::Zero(), 1e-6);
    std::vector<Particle> particles = {p1, p2};
    int num_steps = 1;
    int dt = 0.0001;
    // Run the simulation
    for (int step = 0; step < num_steps; ++step) {
        for (Particle& p : particles) {
            p.SumAccelerations(particles, 1e-5);
        }
        for (Particle& p : particles) {
            p.update(dt);
        }
    }
    Eigen::Vector3d expected_pos(1, 0, 0);
    REQUIRE(particles[1].getPosition().isApprox(expected_pos, 1e-6));
}

TEST_CASE("calculate kinetic energy", "[energy]") {
    Particle p(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 0, 0), 1);
    double expected_energy = 0.5;
    double energy = calcKineticEnergy(p);
    REQUIRE(energy == expected_energy);
}

TEST_CASE("calculate potential energy between two particles", "[energy]") {
    Particle p1(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), 1);
    Particle p2(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), 0.1);
    double expected_energy = -0.05;
    double energy = calcPotentialEnergy(p1, p2);
    REQUIRE(energy == expected_energy);
}

TEST_CASE("calculate total energy of system", "[energy]") {
    Particle p1(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), 1);
    Particle p2(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), 0.1);
    std::vector<Particle> particles = {p1, p2};
    double expected_energy = -0.05;
    double energy = calcTotalEnergy(particles);
    REQUIRE(energy == expected_energy);
}




TEST_CASE("Energy calculation is correct (parallel vs serial)", "[energy]") {
    const int num_particles = 100;
    const double dt = 0.001;
    const double total_time = 10.0;
    const int num_steps = static_cast<int>(total_time / dt);
    const double epsilon = 1e-5;

    // Generate particles
    std::vector<Particle> particles = create_initial_particles();


    // Calculate initial energies using serial version
    double initial_kinetic_energy_serial = 0.0;
    double initial_potential_energy_serial = 0.0;
    double initial_total_energy_serial = 0.0;
    for (const Particle& p : particles) {
        initial_kinetic_energy_serial += calcKineticEnergy(p);
    }
    for (size_t i = 0; i < particles.size(); ++i) {
        const Particle& p1 = particles[i];
        for (size_t j = i + 1; j < particles.size(); ++j) {
            const Particle& p2 = particles[j];
            initial_potential_energy_serial += calcPotentialEnergy(p1, p2);
        }
    }
    initial_total_energy_serial = initial_kinetic_energy_serial + initial_potential_energy_serial;

    // Calculate final energies using serial version
    std::vector<Particle> particles_serial(particles);
    for (int step = 0; step < num_steps; ++step) {
        for (size_t i = 0; i < particles_serial.size(); ++i) {
            particles_serial[i].SumAccelerations(particles_serial, epsilon);
        }
        for (size_t i = 0; i < particles_serial.size(); ++i) {
            particles_serial[i].update(dt);
        }
    }
    double final_kinetic_energy_serial = 0.0;
    double final_potential_energy_serial = 0.0;
    double final_total_energy_serial = 0.0;
    for (const Particle& p : particles_serial) {
        final_kinetic_energy_serial += calcKineticEnergy(p);
    }
    for (size_t i = 0; i < particles_serial.size(); ++i) {
        const Particle& p1 = particles_serial[i];
        for (size_t j = i + 1; j < particles_serial.size(); ++j) {
            const Particle& p2 = particles_serial[j];
            final_potential_energy_serial += calcPotentialEnergy(p1, p2);
        }
    }
    final_total_energy_serial = final_kinetic_energy_serial + final_potential_energy_serial;
    // Calculate final energies using parallel version
    std::vector<Particle> particles_parallel(particles);
    for (int step = 0; step < num_steps; ++step) {
        #pragma omp parallel for
        for (size_t i = 0; i < particles_parallel.size(); ++i) {
            particles_parallel[i].SumAccelerations(particles_parallel, epsilon);
        }
        #pragma omp parallel for
        for (size_t i = 0; i < particles_parallel.size(); ++i) {
            particles_parallel[i].update(dt);
        }
    }
    double final_kinetic_energy_parallel = 0.0;
    double final_potential_energy_parallel = 0.0;
    double final_total_energy_parallel = 0.0;
    #pragma omp parallel for reduction(+:final_kinetic_energy_parallel, final_potential_energy_parallel)
    for (size_t i = 0; i < particles_parallel.size(); ++i) {
        final_kinetic_energy_parallel += calcKineticEnergy(particles_parallel[i]);
        for (size_t j = i + 1; j < particles_parallel.size(); ++j) {
            final_potential_energy_parallel += calcPotentialEnergy(particles_parallel[i], particles_parallel[j]);
        }
    }
    final_total_energy_parallel = final_kinetic_energy_parallel + final_potential_energy_parallel;

    // Compare energies
    REQUIRE(std::abs(final_kinetic_energy_parallel - final_kinetic_energy_serial) < 1e-10);
    REQUIRE(std::abs(final_potential_energy_parallel - final_potential_energy_serial) < 1e-10);
    REQUIRE(std::abs(final_total_energy_parallel - final_total_energy_serial) < 1e-10);
}
