#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <array>
#include <random>

// Include HLS headers
#include "treecon_config_hls.h"
#include "treecon_types_hls.h"
#include "peano_hilbert.h"

const int TREE_DEPTH = 2048;
const int PARTICLE_DEPTH = 1024;

// Forward declaration of HLS kernel
void create_bhtree_kernel(const particle_t *particles,
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
              << "last=" << (node.is_last ? "true" : "false") << " "
              << "next_sibling=" << int(node.next_sibling) << std::endl;
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

bool test_peano_hilbert_key() {
// Create Peano-Hilbert generator with depth 10
PeanoHilbert ph(10);
    
// Single position key generation
uint64_t key = ph.generate_key(100, 200, 300);
std::cout << "Key for position (100, 200, 300): 0x" << std::hex << key << std::dec << std::endl;

// Multiple positions
std::vector<std::array<uint32_t, 3>> positions = {
    {0, 0, 0},
    {1000, 2000, 3000},
    {0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
    {
        static_cast<uint32_t>(0.4 * static_cast<double>(UINT32_MAX)),
        static_cast<uint32_t>(0.8 * static_cast<double>(UINT32_MAX)),
        static_cast<uint32_t>(0.3 * static_cast<double>(UINT32_MAX))
    },
    {1710436918, 3435321836, 4283210188}
};

std::vector<uint64_t> keys = ph.generate_keys(positions);

// Ensure keys match golden values
std::vector<uint64_t> golden_keys = {
    0x0,
    0x0,
    0x29249249, 
    0xdf03f03, 
    0x310f1b87,
};

bool all_match = true;

for (size_t i = 0; i < keys.size(); ++i) {
    if (keys[i] != golden_keys[i]) {
        std::cout << "PH key mismatch at position " << i << ": " << keys[i] << " != " << golden_keys[i] << std::endl;
        all_match = false;
    }
}

if(all_match) {
    std::cout << "All keys match golden values" << std::endl;
} else {
    std::cout << "Keys do not match golden values" << std::endl;
}

    return all_match;
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
    particle_t particles[PARTICLE_DEPTH];
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int orig_idx = key_index_pairs[i].second;
        uint32_t key = key_index_pairs[i].first;
        particle_t p = create_test_particle(test_particles[orig_idx].x,
                                           test_particles[orig_idx].y,
                                           test_particles[orig_idx].z,
                                           test_particles[orig_idx].mass,
                                           key, i);
        particles[i] = p;
        print_particle(p, i);
    }
    
    std::vector<ap_uint<512>> tree_output(TREE_DEPTH); // Small buffer for simple test
    
    std::cout << "\nRunning HLS kernel on manual particles..." << std::endl;
    create_bhtree_kernel(particles, tree_output.data(), NUM_PARTICLES);
    
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
    std::vector<TestParticle> test_particles(PARTICLE_DEPTH);
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
    particle_t particles[PARTICLE_DEPTH];
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int orig_idx = key_index_pairs[i].second;
        uint32_t key = key_index_pairs[i].first;
        particle_t p = create_test_particle(test_particles[orig_idx].x,
                                           test_particles[orig_idx].y,
                                           test_particles[orig_idx].z,
                                           test_particles[orig_idx].mass,
                                           key, i);
        particles[i] = p;
        print_particle(p, i);
    }
    std::vector<ap_uint<512>> tree_output(TREE_DEPTH); // Small buffer for simple test
    
    std::cout << "\nRunning HLS kernel on manual particles..." << std::endl;
    create_bhtree_kernel(particles, tree_output.data(), NUM_PARTICLES);

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
    bool test_tmp;

    // Number of nodes should be 19
    test_tmp = (num_nodes == 19);
    if(!test_tmp) {
        std::cout << "ERROR: Number of nodes should be 19, got " << num_nodes << std::endl;
    }
    test_passed &= test_tmp;

    // Sum of leaf node num_particles should be equal to NUM_PARTICLES
    int num_particles = 0;
    for (int i = 0; i < num_nodes; i++) {
        nodeleaf node = convert_output_node(tree_output[i]);
        if (node.is_leaf) num_particles += node.num_particles;
    }
    test_tmp = (num_particles == NUM_PARTICLES);
    if(!test_tmp) {
        std::cout << "ERROR: Number of particles in leaf nodes should be " << NUM_PARTICLES << ", got " << num_particles << std::endl;
    }
    test_passed &= test_tmp;

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
    test_tmp = true;
    for (int i = 10; i < 19; i++) {
        nodeleaf node = convert_output_node(tree_output[i]);
        test_tmp &= (abs(double(node.pos[0]) - com.x) < COM_TOL);
        test_tmp &= (abs(double(node.pos[1]) - com.y) < COM_TOL);
        test_tmp &= (abs(double(node.pos[2]) - com.z) < COM_TOL);
    }

    if(!test_tmp) {
        std::cout << "ERROR: Node 10 through 18 have different positions!" << std::endl;
    }
    test_passed &= test_tmp;

    return test_passed;
}

bool test_random_particles() {
    std::cout << "\n=== Test: Random Particles ===" << std::endl;
    bool test_passed = true;

    const int NUM_PARTICLES = 1000;
    
    // Create random number generator with fixed seed for reproducibility
    std::mt19937 gen(42);
    std::uniform_real_distribution<float> pos_dis(0.0f, 1.0f);
    
    // Create particles array and fill with random positions
    std::vector<particle_t> particles(NUM_PARTICLES);
    for (int i = 0; i < NUM_PARTICLES; i++) {
        particles[i].pos[0] = pos_dis(gen);
        particles[i].pos[1] = pos_dis(gen);
        particles[i].pos[2] = pos_dis(gen);
        particles[i].mass = 1.0f / NUM_PARTICLES;
        particles[i].idx = i;
        
        // Compute and store PH key
        particles[i].key = simple_peano_hilbert_key(particles[i].pos[0], 
                                                   particles[i].pos[1], 
                                                   particles[i].pos[2]);
    }
    
    // Sort particles by PH key
    std::sort(particles.begin(), particles.end(), 
              [](const particle_t& a, const particle_t& b) {
                  return a.key < b.key;
              });

    // Print first 10 particles
    std::cout << "First 10 particles:" << std::endl;
    for (int i = 0; i < 10; i++) {
        print_particle(particles[i], i);
    }
    std::cout << std::endl;
    
    // Construct tree with kernel
    std::vector<ap_uint<512>> tree_output(TREE_DEPTH);
    create_bhtree_kernel(particles.data(), tree_output.data(), NUM_PARTICLES);
    
    int num_nodes = 0;
    for (int i = 0; i < tree_output.size(); i++) {
        nodeleaf node = convert_output_node(tree_output[i]);
        num_nodes++;
        if (node.is_last) break;
    }

    std::cout << "Number of particles: " << NUM_PARTICLES << std::endl;
    std::cout << "Number of nodes: " << num_nodes << std::endl;

    // Reverse tree output
    std::reverse(tree_output.begin(), tree_output.begin() + num_nodes);

    // Print first 10 nodes
    std::cout << "First 10 nodes:" << std::endl;
    for (int i = 0; i < 10; i++) {
        nodeleaf node = convert_output_node(tree_output[i]);
        print_node(node, i);
    }
    std::cout << std::endl;

    // return false;

    // Check that the next sibling pointers are correct
    std::cout << "\nChecking next sibling pointers..." << std::endl;
    int next_sibling = 0;
    bool sibling_test_passed = true;
    for (int i = 0; i < num_nodes; i++) {
        nodeleaf node = convert_output_node(tree_output[i]);
        if(node.level > 1) continue;

        // check that the next sibling is correct
        sibling_test_passed &= (next_sibling == i);
        if(next_sibling == i) {
            std::cout << "Node " << i << " next sibling is correct, level=" << node.level << std::endl;
        } else {
            std::cout << "Node " << i << " next sibling is incorrect, expected=" << next_sibling << ", level=" << node.level << std::endl;
        }

        // print_node(node, i);

        if(node.next_sibling != -1u) {
            next_sibling = i + node.next_sibling;
        }
        else {
            next_sibling = -1u;
        }
    }

    std::cout << "Last next sibling = " << next_sibling << std::endl;
    sibling_test_passed &= (next_sibling == -1u);

    test_passed &= sibling_test_passed;

    // check sibling of a max depth node
    for(int i = 0; i < num_nodes; i++) {
        nodeleaf node = convert_output_node(tree_output[i]);
        if(node.level == MAX_DEPTH) {
            print_node(node, i);
            break;
        }
    }

    return test_passed;
}

int main() {
    std::cout << "=== HLS BHTree Manual Testbench ===" << std::endl;
    std::cout << "MAX_DEPTH = " << MAX_DEPTH << std::endl;
    std::cout << "NLEAF = " << NLEAF << std::endl;
    std::cout << "sizeof(nodeleaf) = " << sizeof(nodeleaf) << " bytes" << std::endl;
    
    bool test0, test1, test2, test3, test_passed;
    test0 = test_peano_hilbert_key();
    test1 = test_simple_manual_particles();
    test2 = test_at_max_depth();
    test3 = test_random_particles();
    test_passed = test0 && test1 && test2 && test3;

    if(!test0) {
        std::cout << "❌ Test peano hilbert key FAILED!" << std::endl;
    }
    if(!test1) {
        std::cout << "❌ Test simple manual particles FAILED!" << std::endl;
    }
    if(!test2) {
        std::cout << "❌ Test at max depth FAILED!" << std::endl;
    }
    if(!test3) {
        std::cout << "❌ Test random particles FAILED!" << std::endl;
    }

    
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
