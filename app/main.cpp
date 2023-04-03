#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include "particle.hpp"
#include <omp.h>
#include <chrono>
#include <memory>

void print_help() {
    std::cout << "Usage: build/solarSystemSimulator [options]" << std::endl;
    std::cout << "  -h, --help            Help message" << std::endl;
    std::cout << "  -t, --timestep        Time step" << std::endl;
    std::cout << "  -n, --numsteps        Number of time steps" << std::endl;
}

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

int main(int argc, char *argv[]) {
    double dt = 0.0001;
    double total_time =  2  * M_PI;
    int num_steps = total_time/dt;
    if (argc == 1) {
        print_help();
        return 0;
    }

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            print_help();
            return 0;
        } else if (arg == "-t" || arg == "--timestep") {
            if (i + 1 < argc) {
                dt = std::stod(argv[++i]);
                num_steps = total_time/dt;
            } else {
                std::cerr << "Error: " << arg << " requires an argument" << std::endl;
                return 1;
            }
        } else if (arg == "-n" || arg == "--numsteps") {
            if (i + 1 < argc) {
                num_steps = std::stoi(argv[++i]);

            } else {
                std::cerr << "Error: " << arg << " requires an argument" << std::endl;
                return 1;
            }
        } else {
            std::cerr << "Error: Invalid argument '" << arg << "'" << std::endl;
            print_help();
            return 1;
        }
    }

    // Initialize particles
    std::vector<Particle> particles = create_initial_particles();

    // Print initial positions
    std::cout << "Initial positions:" << std::endl;
    for (const Particle& p : particles) {
        std::cout << p.getPosition().transpose() << std::endl;
    }
    for (int step = 0; step < num_steps; ++step) {
        // Update accelerations
        for (Particle& p : particles) {
            p.SumAccelerations(particles, 1e-5);
        }
        
        // Update positions and velocities
        for (Particle& p : particles) {
            p.update(dt);
        }
        
    }
    std::cout << "Final positions:" << std::endl;
    for (const Particle& p : particles) {
        std::cout << p.getPosition().transpose() << std::endl;
    }
}