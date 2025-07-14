#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>

// Include HLS headers
#include "../include/bhtree_config_hls.h"
#include "../include/bhtree_types_hls.h"
#include "hls_stream.h"

// Forward declaration of HLS kernel
void create_bhtree_kernel(hls::stream<particle_t> &particle_stream,
                          ap_uint<512> *tree,
                          count_t num_particles);

// Test particle structure
struct TestParticle {
    double x, y, z, mass;
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

// Helper function to convert ap_uint<512> back to nodeleaf for analysis
nodeleaf convert_output_node(const ap_uint<512>& node_bits) {
    return *reinterpret_cast<const nodeleaf*>(&node_bits);
}

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
              << "last=" << (node.is_last ? "true" : "false") << std::endl;
}

// Simple function to compute Peano-Hilbert key (simplified 3D version)
uint32_t simple_peano_hilbert_key(double x, double y, double z, int depth = 10) {
    // Convert to integer coordinates (0 to 2^depth - 1)
    uint32_t max_coord = (1U << depth) - 1;
    uint32_t ix = (uint32_t)(x * max_coord);
    uint32_t iy = (uint32_t)(y * max_coord);
    uint32_t iz = (uint32_t)(z * max_coord);
    
    // Very simplified PH key - just interleave bits for now
    // This is not a true Peano-Hilbert curve but will work for basic testing
    uint32_t key = 0;
    for (int i = 0; i < depth; i++) {
        uint32_t bit_x = (ix >> i) & 1;
        uint32_t bit_y = (iy >> i) & 1;
        uint32_t bit_z = (iz >> i) & 1;
        key |= (bit_x << (3*i)) | (bit_y << (3*i + 1)) | (bit_z << (3*i + 2));
    }
    return key;
}

// Test with manually created particles 
bool test_simple_manual_particles() {
    std::cout << "\n=== Test: Simple Manual Particles ===" << std::endl;
    
    // Create a very simple test case with known positions
    const int NUM_PARTICLES = 4;
    
    // Manually create particles at different octants
    std::vector<TestParticle> test_particles = {
        {0.1, 0.1, 0.1, 0.25},  // octant 0
        {0.9, 0.1, 0.1, 0.25},  // octant 1  
        {0.1, 0.9, 0.1, 0.25},  // octant 2
        {0.9, 0.9, 0.9, 0.25}   // octant 7
    };
    
    // Compute PH keys and create sorted list
    std::vector<std::pair<uint32_t, int>> key_index_pairs;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        uint32_t ph_key = simple_peano_hilbert_key(test_particles[i].x, 
                                                   test_particles[i].y, 
                                                   test_particles[i].z);
        key_index_pairs.push_back({ph_key, i});
    }
    
    // Sort by PH key
    std::sort(key_index_pairs.begin(), key_index_pairs.end());
    
    std::cout << "Manual particles sorted by PH keys:" << std::endl;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int orig_idx = key_index_pairs[i].second;
        uint32_t key = key_index_pairs[i].first;
        std::cout << "  Particle " << i << " (orig " << orig_idx << "): pos=(" 
                  << test_particles[orig_idx].x << "," 
                  << test_particles[orig_idx].y << "," 
                  << test_particles[orig_idx].z 
                  << ") key=0x" << std::hex << key << std::dec << std::endl;
    }
    
    // Convert to HLS format and run kernel
    hls::stream<particle_t> particle_stream("manual_test");
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int orig_idx = key_index_pairs[i].second;
        uint32_t key = key_index_pairs[i].first;
        particle_t p = create_test_particle(test_particles[orig_idx].x,
                                           test_particles[orig_idx].y,
                                           test_particles[orig_idx].z,
                                           test_particles[orig_idx].mass,
                                           key, i);
        particle_stream.write(p);
        print_particle(p, i);
    }
    
    std::vector<ap_uint<512>> tree_output(50); // Small buffer for simple test
    
    std::cout << "\nRunning HLS kernel on manual particles..." << std::endl;
    create_bhtree_kernel(particle_stream, tree_output.data(), NUM_PARTICLES);
    
    // Analyze output
    int num_nodes = 0;
    std::cout << "\nOutput nodes:" << std::endl;
    for (int i = 0; i < 50; i++) {
        nodeleaf node = convert_output_node(tree_output[i]);
        
        // Stop if we encounter a node with zero mass and zero particles (likely uninitialized)
        if (node.mass == 0.0 && node.num_particles == 0) {
            break;
        }
        
        num_nodes++;
        print_node(node, i);
        if (node.is_last) break;
    }
    
    std::cout << "Simple test produced " << num_nodes << " nodes" << std::endl;
    
    // Basic sanity checks
    bool test_passed = true;
    
    if (num_nodes == 0) {
        std::cout << "ERROR: No output nodes generated!" << std::endl;
        test_passed = false;
    }
    
    // Check that masses are reasonable
    double total_mass = 0.0;
    for (int i = 0; i < num_nodes; i++) {
        nodeleaf node = convert_output_node(tree_output[i]);
        total_mass += double(node.mass);
    }
    
    std::cout << "Total mass in tree: " << total_mass << " (expected: 1.0)" << std::endl;
    if (std::abs(total_mass - 1.0) > 0.01) {
      std::cout << "ERROR: Total mass mismatch!" << std::endl;
      test_passed = false;
    }

    // Now we loop over the nodes and check that their center of mass matches
    // the position of their respective particle
    for (int i = 0; i < num_nodes; i++) {
      nodeleaf node = convert_output_node(tree_output[i]);
      if (abs(double(node.pos[0]) - test_particles[i].x) > 0.01 ||
          abs(double(node.pos[1]) - test_particles[i].y) > 0.01 ||
          abs(double(node.pos[2]) - test_particles[i].z) > 0.01) {
        std::cout << "ERROR: Node " << i
                  << " center of mass does not match particle position!"
                  << std::endl;
        test_passed = false;
      }
      else {
        std::cout << "Node " << i << " center of mass matches particle position!" << std::endl;
      }
    }

    return test_passed;
}

bool test_at_max_depth() {
    std::cout << "\n=== Test: At Max Depth ===" << std::endl;

    const int NUM_PARTICLES = 10;

    // Manually create particles at different octants
    std::vector<TestParticle> test_particles(NUM_PARTICLES);
    double step = (1.0 / pow(2, MAX_DEPTH));
    for (int i = 0; i < NUM_PARTICLES; i++) {
        test_particles[i] = {0.1 + i*step/NUM_PARTICLES, 0.1, 0.1, 1.0/double(NUM_PARTICLES)};
    }
    
    // Compute PH keys and create sorted list
    std::vector<std::pair<uint32_t, int>> key_index_pairs;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        uint32_t ph_key = simple_peano_hilbert_key(test_particles[i].x, 
                                                   test_particles[i].y, 
                                                   test_particles[i].z);
        key_index_pairs.push_back({ph_key, i});
    }
    
    // Sort by PH key
    std::sort(key_index_pairs.begin(), key_index_pairs.end());
    
    std::cout << "Manual particles sorted by PH keys:" << std::endl;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int orig_idx = key_index_pairs[i].second;
        uint32_t key = key_index_pairs[i].first;
        std::cout << "  Particle " << i << " (orig " << orig_idx << "): pos=(" 
                  << test_particles[orig_idx].x << "," 
                  << test_particles[orig_idx].y << "," 
                  << test_particles[orig_idx].z 
                  << ") key=0x" << std::hex << key << std::dec << std::endl;
    }
    
    // Convert to HLS format and run kernel
    hls::stream<particle_t> particle_stream("manual_test");
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int orig_idx = key_index_pairs[i].second;
        uint32_t key = key_index_pairs[i].first;
        particle_t p = create_test_particle(test_particles[orig_idx].x,
                                           test_particles[orig_idx].y,
                                           test_particles[orig_idx].z,
                                           test_particles[orig_idx].mass,
                                           key, i);
        particle_stream.write(p);
        print_particle(p, i);
    }
    std::vector<ap_uint<512>> tree_output(50); // Small buffer for simple test
    
    std::cout << "\nRunning HLS kernel on manual particles..." << std::endl;
    create_bhtree_kernel(particle_stream, tree_output.data(), NUM_PARTICLES);

    // Print output output
    int num_nodes = 0;
    std::cout << "\nOutput nodes:" << std::endl;
    for (int i = 0; i < 50; i++) {
        nodeleaf node = convert_output_node(tree_output[i]);
        print_node(node, i);
        num_nodes++;
        if (node.is_last) break;
    }

    // Now some sanity checks
    bool test_passed = true;

    // Number of nodes should be 11
    test_passed &= (num_nodes == 11);

    // Sum of leaf node num_particles should be equal to NUM_PARTICLES
    int num_particles = 0;
    for (int i = 0; i < num_nodes; i++) {
        nodeleaf node = convert_output_node(tree_output[i]);
        if (node.is_leaf) num_particles += node.num_particles;
    }
    test_passed &= (num_particles == NUM_PARTICLES);

    // Check that node 2 through 10 have same pos and that their pos match total center of mass
    // Compute center of mass
    TestParticle com = {0.0, 0.0, 0.0, 0.0};
    TestParticle com_node0 = {0.0, 0.0, 0.0, 0.0};
    TestParticle com_node1 = {0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < NUM_PARTICLES; i++) {
        com.x += test_particles[i].x * test_particles[i].mass;
        com.y += test_particles[i].y * test_particles[i].mass;
        com.z += test_particles[i].z * test_particles[i].mass;
        com.mass += test_particles[i].mass;

        if(i < 8) {
            com_node0.x += test_particles[i].x * test_particles[i].mass;
            com_node0.y += test_particles[i].y * test_particles[i].mass;
            com_node0.z += test_particles[i].z * test_particles[i].mass;
            com_node0.mass += test_particles[i].mass;
        }
        else {
            com_node1.x += test_particles[i].x * test_particles[i].mass;
            com_node1.y += test_particles[i].y * test_particles[i].mass;
            com_node1.z += test_particles[i].z * test_particles[i].mass;
            com_node1.mass += test_particles[i].mass;
        }
    }
    com.x /= com.mass;
    com.y /= com.mass;
    com.z /= com.mass;
    com_node0.x /= com_node0.mass;
    com_node0.y /= com_node0.mass;
    com_node0.z /= com_node0.mass;
    com_node1.x /= com_node1.mass;
    com_node1.y /= com_node1.mass;
    com_node1.z /= com_node1.mass;

    // Check that node 2 through 10 have same pos and that their pos match total center of mass to within 0.1%
    const double COM_TOL = 1e-6;
    for (int i = 2; i < 11; i++) {
        nodeleaf node = convert_output_node(tree_output[i]);
        test_passed &= (abs(double(node.pos[0]) - com.x) < COM_TOL);
        test_passed &= (abs(double(node.pos[1]) - com.y) < COM_TOL);
        test_passed &= (abs(double(node.pos[2]) - com.z) < COM_TOL);
    }

    // Check that node 0 has the same center of mass as the first 8 particles and node 1 has the same center of mass as the last 2 particles
    nodeleaf node_0 = convert_output_node(tree_output[0]);
    test_passed &= (abs(double(node_0.pos[0]) - com_node0.x) < COM_TOL);
    test_passed &= (abs(double(node_0.pos[1]) - com_node0.y) < COM_TOL);
    test_passed &= (abs(double(node_0.pos[2]) - com_node0.z) < COM_TOL);
    nodeleaf node_1 = convert_output_node(tree_output[1]);
    test_passed &= (abs(double(node_1.pos[0]) - com_node1.x) < COM_TOL);
    test_passed &= (abs(double(node_1.pos[1]) - com_node1.y) < COM_TOL);
    test_passed &= (abs(double(node_1.pos[2]) - com_node1.z) < COM_TOL);

    return test_passed;
}

int main() {
    std::cout << "=== HLS BHTree Manual Testbench ===" << std::endl;
    std::cout << "MAX_DEPTH = " << MAX_DEPTH << std::endl;
    std::cout << "NLEAF = " << NLEAF << std::endl;
    std::cout << "sizeof(nodeleaf) = " << sizeof(nodeleaf) << " bytes" << std::endl;
    
    bool test_passed = test_simple_manual_particles();

    test_passed &= test_at_max_depth();
    
    // Summary
    std::cout << "\n=== Test Summary ===" << std::endl;
    if (test_passed) {
        std::cout << "✅ Manual particle test PASSED!" << std::endl;
        return 0;
    } else {
        std::cout << "❌ Manual particle test FAILED!" << std::endl;
        return 1;
    }
}
