#ifndef TREE_IO_H
#define TREE_IO_H

#include "tree_config.h"
#include "tree_types.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

// File I/O functions for particles and tree data
void save_particles_and_tree(const std::vector<particle_t>& particles, 
                             const std::vector<nodeleaf>& tree,
                             const std::string& filename_prefix) {
    // Create full path with data directory using git repo base
    std::string base_path = "/home/abeane/Projects/BHTreeFPGA/hw/tb/data/" + filename_prefix;
    
    // Save particles
    std::ofstream pfile(base_path + "_particles.bin", std::ios::binary);
    size_t num_particles = particles.size();
    pfile.write(reinterpret_cast<const char*>(&num_particles), sizeof(num_particles));
    pfile.write(reinterpret_cast<const char*>(particles.data()), num_particles * sizeof(particle_t));
    pfile.close();
    
    // Save tree
    std::ofstream tfile(base_path + "_tree.bin", std::ios::binary);
    size_t num_nodes = tree.size();
    tfile.write(reinterpret_cast<const char*>(&num_nodes), sizeof(num_nodes));
    tfile.write(reinterpret_cast<const char*>(tree.data()), num_nodes * sizeof(nodeleaf));
    tfile.close();
    
    std::cout << "Saved " << num_particles << " particles and " << num_nodes << " tree nodes to " << base_path << "_*.bin" << std::endl;
}

std::pair<std::vector<particle_t>, std::vector<nodeleaf>> load_particles_and_tree(const std::string& filename_prefix) {
    std::vector<particle_t> particles;
    std::vector<nodeleaf> tree;
    
    // Create full path with data directory using git repo base
    std::string base_path = "/home/abeane/Projects/BHTreeFPGA/hw/tb/data/" + filename_prefix;
    
    // Load particles
    std::ifstream pfile(base_path + "_particles.bin", std::ios::binary);
    if (pfile) {
        size_t num_particles;
        pfile.read(reinterpret_cast<char*>(&num_particles), sizeof(num_particles));
        particles.resize(num_particles);
        pfile.read(reinterpret_cast<char*>(particles.data()), num_particles * sizeof(particle_t));
        pfile.close();
    }
    
    // Load tree
    std::ifstream tfile(base_path + "_tree.bin", std::ios::binary);
    if (tfile) {
        size_t num_nodes;
        tfile.read(reinterpret_cast<char*>(&num_nodes), sizeof(num_nodes));
        tree.resize(num_nodes);
        tfile.read(reinterpret_cast<char*>(tree.data()), num_nodes * sizeof(nodeleaf));
        tfile.close();
    }
    
    std::cout << "Loaded " << particles.size() << " particles and " << tree.size() << " tree nodes from " << base_path << "_*.bin" << std::endl;
    return {particles, tree};
}

#endif
