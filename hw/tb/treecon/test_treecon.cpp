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
#include "tree_io.h"

const int TREE_DEPTH = 16384;
const int PARTICLE_DEPTH = 8192;

// Test with manually created particles
bool test_simple_manual_particles(bool verbose) {
  std::cout << "\n=== Test: Simple Manual Particles ===" << std::endl;

  auto result = generate_simple_tree(MAX_DEPTH, verbose);
  std::vector<nodeleaf> tree = result.tree;
  std::vector<particle_t> particles = result.particles;

  std::cout << "Simple test produced " << tree.size() << " nodes" << std::endl;

  // Basic sanity checks
  bool test_passed = true;

  if (tree.size() == 0) {
    std::cout << "ERROR: No output nodes generated!" << std::endl;
    test_passed = false;
  }

  // Check that masses are reasonable
  double total_mass = 0.0;
  for (int i = 0; i < tree.size(); i++) {
    total_mass += double(tree[i].mass);
  }

  if (verbose) {
    std::cout << "Total mass in tree: " << total_mass << " (expected: 1.0)" << std::endl;
  }
    if (std::abs(total_mass - 1.0) > 0.01) {
      test_passed = false;
    }

    // Now we loop over the nodes and check that their center of mass matches
    // the position of their respective particle
    for (int i = 0; i < tree.size(); i++) {
      int part_idx = tree.size() - i - 1;
    nodeleaf node = tree[i];
    if (abs(double(node.pos[0] - particles[part_idx].pos[0])) > 0.01 ||
        abs(double(node.pos[1] - particles[part_idx].pos[1])) > 0.01 ||
        abs(double(node.pos[2] - particles[part_idx].pos[2])) > 0.01) {
      if (verbose) {
        std::cout << "ERROR: Node " << i
                  << " center of mass does not match particle position!"
                  << std::endl;
        print_node(node, i);
      }
      test_passed = false;
    } else {
      if (verbose) {
        std::cout << "Node " << i
                  << " center of mass matches particle position!" << std::endl;
      }
    }
    }

    return test_passed;
}

bool test_at_max_depth(bool verbose) {
    std::cout << "\n=== Test: At Max Depth ===" << std::endl;
    
    auto result = generate_max_depth_tree(MAX_DEPTH, verbose);
    std::vector<nodeleaf> tree = result.tree;
    std::vector<particle_t> particles = result.particles;
    


    if (verbose) {
        std::cout << "\nOutput nodes:" << std::endl;
        for (int i = 0; i < std::min(50, (int)tree.size()); i++) {
            print_node(tree[i], i);
        }
    }

    // Now some sanity checks
    bool test_passed = true;
    bool test_tmp = true;

    // Number of nodes should be 19
    test_tmp = (tree.size() == 19);
    if(!test_tmp) {
        std::cout << "ERROR: Number of nodes should be 19, got " << tree.size() << std::endl;
    }
    test_passed &= test_tmp;

    // Sum of leaf node num_particles should be equal to NUM_PARTICLES
    int num_particles_in_leaves = 0;
    for (int i = 0; i < tree.size(); i++) {
        nodeleaf node = tree[i];
        if (node.is_leaf) num_particles_in_leaves += node.num_particles;
    }
    test_tmp = (num_particles_in_leaves == particles.size());
    if(!test_tmp) {
        std::cout << "ERROR: Number of particles in leaf nodes should be " << particles.size() << ", got " << num_particles_in_leaves << std::endl;
    }
    test_passed &= test_tmp;

    // Check that node 10 through 18 have same pos and that their pos match total center of mass
    // Compute center of mass from particles
    double com_x = 0.0, com_y = 0.0, com_z = 0.0, total_mass = 0.0;
    double com_node0_x = 0.0, com_node0_y = 0.0, com_node0_z = 0.0, mass_node0 = 0.0;
    double com_node1_x = 0.0, com_node1_y = 0.0, com_node1_z = 0.0, mass_node1 = 0.0;
    
    for (int i = 0; i < particles.size(); i++) {
        double x = double(particles[i].pos[0]);
        double y = double(particles[i].pos[1]);
        double z = double(particles[i].pos[2]);
        double mass = double(particles[i].mass);
        
        com_x += x * mass;
        com_y += y * mass;
        com_z += z * mass;
        total_mass += mass;

    }
    
    com_x /= total_mass;
    com_y /= total_mass;
    com_z /= total_mass;

    // Check that node 10 through 18 have same pos and that their pos match total center of mass to within tolerance
    const double COM_TOL = 1e-6;
    test_tmp = true;
    for (int i = 0; i < MAX_DEPTH-1; i++) {
        nodeleaf node = tree[i];
        test_tmp &= (std::abs(double(node.pos[0]) - com_x) < COM_TOL);
        test_tmp &= (std::abs(double(node.pos[1]) - com_y) < COM_TOL);
        test_tmp &= (std::abs(double(node.pos[2]) - com_z) < COM_TOL);
    }

    if(!test_tmp) {
        std::cout << "ERROR: Node 0 through " << MAX_DEPTH-1 << " have different positions!" << std::endl;
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

    // Print level 1 nodes
    // if (verbose) {
        std::cout << "Level 1 nodes:" << std::endl;
        for (int i = 0; i < tree.size(); i++) {
            nodeleaf node = tree[i];
            if (node.level == 1) {
                print_node(node, i);
            }
        }
    // }

    save_particles_and_tree(particles, tree, "random_test_data");

    return test_passed;
}

bool test_quasirandom_particles(bool verbose) {
    std::cout << "\n=== Test: Quasirandom Particles ===" << std::endl;
    bool test_passed = true;
    const int NUM_PARTICLES = 8000;
    double pos0[3] = {0.5, 0.5, 0.5};
    auto result = generate_quasirandom_tree(NUM_PARTICLES, MAX_DEPTH, verbose, pos0);
    
    std::vector<nodeleaf> tree = result.tree;
    std::vector<particle_t> particles = result.particles;

    save_particles_and_tree(particles, tree, "quasirandom_test_data");

    // check expected number of nodes
    int expected_nodes = 11487;
    bool nodes_test_passed = (tree.size() == expected_nodes);
    if (verbose && !nodes_test_passed) {
      std::cout << "❌ Number of nodes test for quasirandom particles FAILED" << std::endl;
      std::cout << "Number of nodes: " << tree.size() << std::endl;
      std::cout << "Expected number of nodes: " << expected_nodes << std::endl;
    }
    test_passed &= nodes_test_passed;

    std::cout << "Quasirandom particles test " << (test_passed ? "✅ PASSED" : "❌ FAILED") << std::endl;

    return test_passed;
}

int main() {
    std::cout << "=== HLS BHTree Manual Testbench ===" << std::endl;
    std::cout << "MAX_DEPTH = " << MAX_DEPTH << std::endl;
    std::cout << "NLEAF = " << NLEAF << std::endl;
    std::cout << "sizeof(nodeleaf) = " << sizeof(nodeleaf) << " bytes" << std::endl;

    bool test0, test1, test2, test3, test4, test_passed;

    test0 = test_peano_hilbert_key();
    test1 = test_simple_manual_particles(false);
    test2 = test_at_max_depth(false);
    test3 = test_random_particles(false);
    test4 = test_quasirandom_particles(false);
    test_passed = test0 && test1 && test2 && test3 && test4;


    std::cout << std::endl;

    if(!test0) {
        std::cout << "❌ Test peano hilbert key FAILED!" << std::endl;
    }
    if(!test1) {
        std::cout << "❌ Test simple manual particles FAILED!" << std::endl;
        test_simple_manual_particles(true);
    }
    if(!test2) {
        std::cout << "❌ Test at max depth FAILED!" << std::endl;
        test2 = test_at_max_depth(true);
    }
    if(!test3) {
        std::cout << "❌ Test random particles FAILED!" << std::endl;
        test_random_particles(true);
    }

    if(!test4) {
        std::cout << "❌ Test quasirandom particles FAILED!" << std::endl;
        test_quasirandom_particles(true);
    }

    // Summary
    std::cout << "\n=== Test Summary ===" << std::endl;
    if (test_passed) {
        std::cout << "✅ Tree construction tests PASSED!" << std::endl << std::endl;
        return 0;
    } else {
        std::cout << "❌ Tree construction tests FAILED!" << std::endl << std::endl;
        return 1;
    }
}
