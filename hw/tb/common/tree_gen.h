#ifndef TREE_GEN_H
#define TREE_GEN_H

#include "tree_config.h"
#include "tree_types.h"
#include "peano_hilbert.h"
#include "treecon.h"
#include "test_util.h"
#include <random>
#include <algorithm>
#include <iostream>
#include <vector>

struct TreeAndParticles {
    std::vector<nodeleaf> tree;
    std::vector<particle_t> particles;
};

TreeAndParticles generate_random_tree(int num_particles, int max_depth, bool verbose){
    const int NUM_PARTICLES = 1000;
    
    // Create random number generator with fixed seed for reproducibility
    std::mt19937 gen(42);
    std::uniform_real_distribution<float> pos_dis(0.0f, 1.0f);
    
    // Create particles array and fill with random positions
    PeanoHilbert ph(max_depth);
    std::vector<particle_t> particles(NUM_PARTICLES);
    for (int i = 0; i < NUM_PARTICLES; i++) {
        particles[i].pos[0] = pos_dis(gen);
        particles[i].pos[1] = pos_dis(gen);
        particles[i].pos[2] = pos_dis(gen);
        particles[i].mass = 1.0f / NUM_PARTICLES;
        particles[i].idx = i;
        
        // Compute and store PH key
        particles[i].key = ph.generate_key(particles[i].pos[0], 
                                                   particles[i].pos[1], 
                                                   particles[i].pos[2]);
    }
    
    // Sort particles by PH key
    std::sort(particles.begin(), particles.end(), 
              [](const particle_t& a, const particle_t& b) {
                  return a.key < b.key;
              });

    // Print first 10 particles
    if (verbose) {
        std::cout << "First 10 particles:" << std::endl;
    for (int i = 0; i < 10; i++) {
            print_particle(particles[i], i);
        }
        std::cout << std::endl;
    }
    
    // Construct tree with kernel
    std::vector<ap_uint<512>> tree_raw(NUM_PARTICLES * 2);
    create_bhtree_kernel(particles.data(), tree_raw.data(), NUM_PARTICLES);

    int num_nodes = 0;
    std::vector<nodeleaf> tree;
    for (int i = 0; i < tree_raw.size(); i++) {
        nodeleaf node = convert_output_node(tree_raw[i]);
        tree.push_back(node);
        num_nodes++;
        if (node.is_last) break;
    }

    return {tree, particles};
}


#endif