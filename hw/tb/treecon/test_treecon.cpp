#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <array>
#include <random>

// Include HLS headers
#include "peano_hilbert.h"
#include "treecon.h"
#include "test_peano_hilbert.h"
#include "test_util.h"
#include "tree_gen.h"

const int TREE_DEPTH = 2048;
const int PARTICLE_DEPTH = 1024;

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
    PeanoHilbert ph(10);
    std::vector<std::pair<uint32_t, int>> key_index_pairs;
    for (int i = 0; i < NUM_PARTICLES; i++) {
      uint32_t ph_key = ph.generate_key(
          test_particles[i].xint(), test_particles[i].yint(), test_particles[i].zint());
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
    int idx_map[] = {0, 2, 1, 3};
    int part_idx = idx_map[i];
      nodeleaf node = convert_output_node(tree_output[i]);
      if (abs(double(node.pos[0]) - test_particles[part_idx].x) > 0.01 ||
          abs(double(node.pos[1]) - test_particles[part_idx].y) > 0.01 ||
          abs(double(node.pos[2]) - test_particles[part_idx].z) > 0.01) {
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
    PeanoHilbert ph(10);
    std::vector<std::pair<uint32_t, int>> key_index_pairs;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        uint32_t ph_key = ph.generate_key(test_particles[i].xint(), 
                                                   test_particles[i].yint(), 
                                                   test_particles[i].zint());
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

bool test_random_particles(bool verbose) {
    std::cout << "\n=== Test: Random Particles ===" << std::endl;
    bool test_passed = true;

    const int NUM_PARTICLES = 1000;
    
    auto result = generate_random_tree(NUM_PARTICLES, MAX_DEPTH, verbose);
    std::vector<nodeleaf> tree = result.tree;
    std::vector<particle_t> particles = result.particles;

    // check expected number of nodes
    int expected_nodes = 1494;
    bool nodes_test_passed = (tree.size() == expected_nodes);
    if (verbose && !nodes_test_passed) {
      std::cout << "❌ Number of nodes test FAILED" << std::endl;
      std::cout << "Number of nodes: " << tree.size() << std::endl;
      std::cout << "Expected number of nodes: " << expected_nodes << std::endl;
    }
    test_passed &= nodes_test_passed;

    // Check that the next sibling pointers are correct
    if (verbose) {
        std::cout << "\nChecking next sibling pointers..." << std::endl;
    }
    int next_sibling = 0;
    bool sibling_test_passed = true;
    for (int i = 0; i < tree.size(); i++) {
        nodeleaf node = tree[i];
        if(node.level > 1) continue;

        // check that the next sibling is correct
        sibling_test_passed &= (next_sibling == i);
        if (verbose) {
          if (next_sibling == i) {
            std::cout << "Node " << i
                      << " next sibling is correct, level=" << node.level
                      << std::endl;
          } else {
            std::cout << "Node " << i
                      << " next sibling is incorrect, expected=" << next_sibling
                      << ", level=" << node.level << std::endl;
            print_node(node, i);
          }
        }

        if(node.next_sibling != -1u) {
            next_sibling = i + node.next_sibling;
        }
        else {
            next_sibling = -1u;
        }
    }
    if (verbose) {
      std::cout << "Last next sibling = " << next_sibling << std::endl;
    }

    sibling_test_passed &= (next_sibling == -1u);
    test_passed &= sibling_test_passed;

    // check sibling of a max depth node and mass of leaves
    mass_t total_leaf_mass = 0.0;
    bool max_depth_node_test_passed = true;
    for(int i = 0; i < tree.size(); i++) {
        nodeleaf node = tree[i];
        if(node.is_leaf) {
            total_leaf_mass += node.mass;
            if(node.next_sibling != 1 && node.next_sibling != -1u) {
                max_depth_node_test_passed = false;
                if (verbose) {
                    std::cout << "Node " << i
                              << " next sibling is incorrect, expected=1 or -1u, got " << node.next_sibling << std::endl;
                    print_node(node, i);
                }
            }
        }
    }
    bool leaf_mass_test_passed = (std::abs(double(total_leaf_mass) - 1.0) < 1e-4);

    if (verbose) {
      if(leaf_mass_test_passed) {
        std::cout << "Total leaf mass ✅ PASSED = " << double(total_leaf_mass) << std::endl;
      }
      else {
        std::cout << "Total leaf mass ❌ FAILED = " << double(total_leaf_mass) << " (expected: 1.0)" << std::endl;
      }
    }

    test_passed &= leaf_mass_test_passed;

    test_passed &= max_depth_node_test_passed;

    std::cout << "Random particles test " << (test_passed ? "✅ PASSED" : "❌ FAILED") << std::endl;

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
    test3 = test_random_particles(false);
    test_passed = test0 && test1 && test2 && test3;

    if(!test3) {
        test_random_particles(true);
    }

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
