#ifndef DIRECT_SUM_H
#define DIRECT_SUM_H

#include "tree_config.h"
#include "treewalk.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace {
    // Helper function to compute squared distance and direction vector
    struct VectorInfo {
        double dx, dy, dz;
        double r2, r, rinv;
        double rhat[3];
    };
    
    VectorInfo compute_vector_info(const particle_t& p1, const particle_t& p0) {
        VectorInfo info;
        info.dx = static_cast<double>(p1.pos[0] - p0.pos[0]);
        info.dy = static_cast<double>(p1.pos[1] - p0.pos[1]);
        info.dz = static_cast<double>(p1.pos[2] - p0.pos[2]);
        
        info.r2 = info.dx * info.dx + info.dy * info.dy + info.dz * info.dz;
        info.r = std::sqrt(info.r2);
        info.rinv = 1.0 / info.r;
        
        info.rhat[0] = info.dx * info.rinv;
        info.rhat[1] = info.dy * info.rinv;
        info.rhat[2] = info.dz * info.rinv;
        
        return info;
    }
    
    double compute_rhat_magnitude(const double rhat[3]) {
        return std::sqrt(rhat[0] * rhat[0] + rhat[1] * rhat[1] + rhat[2] * rhat[2]);
    }
    
    void print_debug_info(const particle_t& p0, double r_min, double rhat_violation, 
                         double total_mass, const double com[3]) {
        std::cout << "input pos: " << p0.pos[0] << " " << p0.pos[1] << " " << p0.pos[2] << std::endl;
        std::cout << "r_min: " << r_min << std::endl;
        std::cout << "rhat_max_violation: " << rhat_violation << std::endl;
        std::cout << "total_mass: " << total_mass << std::endl;
        std::cout << "com: " << com[0] << " " << com[1] << " " << com[2] << std::endl;
    }
}

particle_t direct_sum_force(std::vector<particle_t> particles, particle_t p0, double eps, bool debug=false) {
    // Reset acceleration
    p0.acc[0] = p0.acc[1] = p0.acc[2] = 0.0;
    
    // Initialize tracking variables
    double total_mass = 0.0;
    double com[3] = {0.0, 0.0, 0.0};
    double r2_min = std::numeric_limits<double>::max();
    double rhat_max_violation = 0.0;
    
    // std::cout << "particles.size(): " << particles.size() << std::endl;
    
    // Compute forces from all particles
    for (const auto& particle : particles) {
        const VectorInfo vec = compute_vector_info(particle, p0);
        const double mass = static_cast<double>(particle.mass);
        
        // Update minimum distance
        r2_min = std::min(r2_min, vec.r2);
        
        // Compute and apply force
        const double W2 = softW2_gadget4(vec.r, vec.rinv, eps);
        double force = mass * W2;
        if(vec.r > 0.0) // guard against pathological case
        {
        p0.acc[0] -= force * vec.rhat[0];
        p0.acc[1] -= force * vec.rhat[1];
        p0.acc[2] -= force * vec.rhat[2];
        }
        
        // Track maximum unit vector violation for debugging
        const double rhat_mag = compute_rhat_magnitude(vec.rhat);
        rhat_max_violation = std::max(rhat_max_violation, std::abs(rhat_mag - 1.0));
        
        // Accumulate center of mass data
        total_mass += mass;
        const double pos_x = static_cast<double>(particle.pos[0]);
        const double pos_y = static_cast<double>(particle.pos[1]);
        const double pos_z = static_cast<double>(particle.pos[2]);
        
        com[0] += pos_x * mass;
        com[1] += pos_y * mass;
        com[2] += pos_z * mass;
    }
    
    // Finalize center of mass
    com[0] /= total_mass;
    com[1] /= total_mass;
    com[2] /= total_mass;
    
    // Output debug information
    if(debug) {
        print_debug_info(p0, std::sqrt(r2_min), rhat_max_violation, total_mass, com);
    }
    
    return p0;
}

#endif