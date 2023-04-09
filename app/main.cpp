#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include "particle.hpp"
#include <omp.h>
#include <chrono>
#include <memory>
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
void print_help() {
    std::cout << "Usage: build/solarSystemSimulator [options]" << std::endl;
    std::cout << "  -h, --help            Help message" << std::endl;
    std::cout << "  -t, --timestep        Time step" << std::endl;
    std::cout << "  -n, --numsteps        Number of time steps" << std::endl;
    std::cout << "  -s, --solar           solarSystemGenerator" << std::endl;
    std::cout << "  -r, --random          randomGenerator with number of particles" << std::endl;
    std::cout << "  -e,                   epsilon" << std::endl;
}



int main(int argc, char *argv[]) {
    double dt = 0.0001;
    double total_time =  2  * M_PI; //When running 100 years should use 2*100*M_PI
    int num_steps = total_time/dt;
    std::unique_ptr<InitialConditionGenerator> generator;
    double epsilon = 1e-5;
    int num_particles = 0;
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
        } else if (arg == "-s" || arg == "--solar") {
            generator = std::make_unique<SolarSystemGenerator>();
        } else if (arg == "-r" || arg == "--random") {
            num_particles = std::stoi(argv[++i]);
            generator = std::make_unique<RandomGenerator>(num_particles);
        } else if (arg == "-e") {
            if (i + 1 < argc) {
                epsilon = std::stod(argv[++i]);
            } else {
                std::cerr << "Error: " << arg << " requires an argument" << std::endl;
                return 1;
            }
        } else if (arg == "--threads") {
            if (i + 1 < argc) {
                int num_threads = std::stoi(argv[++i]);
                omp_set_num_threads(num_threads);
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
    if (!generator) {
        std::cerr << "Error: You should choose a generator" << std::endl;
        return 1;
    }
    // Initialize particles
    std::vector<Particle> particles = generator->generateInitialConditions();

    // Print initial positions
    //std::cout << "Initial positions:" << std::endl;
    /*for (const Particle& p : particles) {
        std::cout << p.getPosition().transpose() << std::endl;
    }*/
    double initial_kinetic_energy = 0.0;
    double initial_potential_energy = 0.0;
    double initial_total_energy = 0.0;

    // Calculate initial energies of the system
    for (const Particle& p : particles) {
        initial_kinetic_energy += calcKineticEnergy(p);
    }
    for (size_t i = 0; i < particles.size(); ++i) {
        const Particle& p1 = particles[i];
        for (size_t j = i + 1; j < particles.size(); ++j) {
            const Particle& p2 = particles[j];
            initial_potential_energy += calcPotentialEnergy(p1, p2);
        }
    }
    initial_total_energy = calcTotalEnergy(particles);

    
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int step = 0; step < num_steps; ++step) {
        std::vector<Particle> particles_copy(particles);
        #pragma omp parallel for 
        for (size_t i = 0; i < particles.size(); ++i) {
            particles[i].SumAccelerations(particles_copy, epsilon);
        }
        #pragma omp parallel for 
        for (size_t i = 0; i < particles.size(); ++i) {
            particles[i].update(dt);
        }
    
}
    auto end_time = std::chrono::high_resolution_clock::now();
    auto total_time1 = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    double avg_time = static_cast<double>(total_time1) / num_steps;

    //std::cout << "Final positions:" << std::endl;
    /* for (const Particle& p : particles) {
        std::cout << p.getPosition().transpose() << std::endl;
    } */
    double final_kinetic_energy = 0.0;
    double final_potential_energy = 0.0;
    double final_total_energy = 0.0;

    // Calculate final energies of the system
    for (const Particle& p : particles) {
        final_kinetic_energy += calcKineticEnergy(p);
    }
    for (size_t i = 0; i < particles.size(); ++i) {
        const Particle& p1 = particles[i];
        for (size_t j = i + 1; j < particles.size(); ++j) {
            const Particle& p2 = particles[j];
            final_potential_energy += calcPotentialEnergy(p1, p2);
        }
    }
    final_total_energy = calcTotalEnergy(particles);
    std::cout << "Initial kinetic energy = " << initial_kinetic_energy << std::endl;
    std::cout << "Initial potential energy = " << initial_potential_energy << std::endl;
    std::cout << "Initial total energy = " << initial_total_energy << std::endl;
    std::cout << "Final kinetic energy = " << final_kinetic_energy << std::endl;
    std::cout << "Final potential energy = " << final_potential_energy << std::endl;
    std::cout << "Final total energy = " << final_total_energy << std::endl;
    double energy_difference = std::abs(final_total_energy - initial_total_energy);
    std::cout << "Energy difference = " << energy_difference << std::endl;
    std::cout << "Timestep = " << dt << ", Total time = " << total_time1 << " ms, Average time = " << avg_time << " ms" << std::endl;

}