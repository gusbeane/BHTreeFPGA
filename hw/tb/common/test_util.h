#ifndef TEST_UTIL_H
#define TEST_UTIL_H

#include "tree_config.h"
#include "tree_types.h"

// Helper function to print a particle for debugging
void print_particle(const particle_t& p, int index) {
    std::cout << "Particle " << index << ": "
              << "pos=(" << double(p.pos[0]) << "," << double(p.pos[1]) << "," << double(p.pos[2]) << ") "
              << "mass=" << double(p.mass) << " "
              << "key=0x" << std::hex << int(p.key) << std::dec << " "
              << "idx=" << int(p.idx) << std::endl;
}

// Helper function to print a node for debugging  
void print_node(const nodeleaf& node, int index) {
    std::cout << "Node " << index << ": "
              << "level=" << int(node.level) << " "
              << "key=0x" << std::hex << int(node.key) << std::dec << " "
              << "npart=" << int(node.num_particles) << " "
              << "start_idx=" << int(node.start_idx) << " "
              << "pos=(" << double(node.pos[0]) << "," << double(node.pos[1]) << "," << double(node.pos[2]) << ") "
              << "mass=" << double(node.mass) << " "
              << "leaf=" << (node.is_leaf ? "true" : "false") << " "
              << "last=" << (node.is_last ? "true" : "false") << " "
              << "next_sibling=" << int(node.next_sibling) << std::endl;
}

// Helper function to convert ap_uint<512> back to nodeleaf for analysis
nodeleaf convert_output_node(const ap_uint<512>& node_bits) {
    return *reinterpret_cast<const nodeleaf*>(&node_bits);
}

// Test particle structure
// Note: Using accessor methods for xint/yint/zint is efficient for small test code and CPU-side logic.
// For large arrays of this struct, especially in performance-critical or hardware contexts, 
// prefer to precompute and store xint/yint/zint in advance to avoid repeated computation and 
// potential cache/memory overhead from method calls and mutable state.

struct TestParticle {
    double x, y, z, mass;
    unsigned int ph_key;

    // Accessors to compute integer coordinates on demand
    unsigned int xint() const {
        constexpr unsigned int max_uint = std::numeric_limits<unsigned int>::max();
        return static_cast<unsigned int>(x * max_uint);
    }
    unsigned int yint() const {
        constexpr unsigned int max_uint = std::numeric_limits<unsigned int>::max();
        return static_cast<unsigned int>(y * max_uint);
    }
    unsigned int zint() const {
        constexpr unsigned int max_uint = std::numeric_limits<unsigned int>::max();
        return static_cast<unsigned int>(z * max_uint);
    }
};

// Helper function to create a test particle
particle_t create_test_particle(double x, double y, double z, double mass, 
                                uint32_t ph_key, int index) {
    particle_t p;
    p.pos[0] = pos_t(x);
    p.pos[1] = pos_t(y); 
    p.pos[2] = pos_t(z);
    p.mass = mass_t(mass);
    p.idx = count_t(index);
    p.key = phkey_t(ph_key & 0x3FFFFFFF); // Mask to 30 bits for phkey_t
    return p;
}

#endif