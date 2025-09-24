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
#include <cmath>

struct TreeAndParticles {
    std::vector<nodeleaf> tree;
    std::vector<particle_t> particles;
};

std::vector<nodeleaf> treecon_wrapper(std::vector<particle_t> particles,
                                      int NUM_PARTICLES, bool verbose) {
  // Construct tree with kernel
  std::vector<ap_uint<512>> tree_raw(NUM_PARTICLES * 3);
  create_bhtree_kernel(particles.data(), tree_raw.data(), NUM_PARTICLES);

  int num_nodes = 0;
  std::vector<nodeleaf> tree;
  for (int i = 0; i < tree_raw.size(); i++) {
    nodeleaf node = convert_output_node(tree_raw[i]);
    tree.push_back(node);
    num_nodes++;
    if (node.is_last)
      break;
  }

  // Reverse tree output
  std::reverse(tree.begin(), tree.begin() + tree.size());

  // Print first 10 nodes
  if (verbose) {
    std::cout << "First 10 nodes:" << std::endl;
    for (int i = 0; i < 10; i++) {
      print_node(tree[i], i);
    }
    std::cout << std::endl;
  }

  return tree;
}

std::vector<particle_t> phsort(std::vector<particle_t> particles, int max_depth,
                               bool verbose) {
  PeanoHilbert ph(max_depth);
  for (int i = 0; i < particles.size(); i++) {
    // Compute and store PH key
    int w = particles[i].pos[0].width;
    uint64_t x_int = uint64_t(particles[i].pos[0].range(w-1, 0));
    uint64_t y_int = uint64_t(particles[i].pos[1].range(w-1, 0));
    uint64_t z_int = uint64_t(particles[i].pos[2].range(w-1, 0));
    particles[i].key = phkey_t(ph.generate_key(x_int, y_int, z_int));
  }

  // Sort particles by PH key
  std::sort(
      particles.begin(), particles.end(),
      [](const particle_t &a, const particle_t &b) { return a.key < b.key; });

  // Print first 10 particles
  if (verbose) {
    std::cout << "First 10 particles:" << std::endl;
    for (int i = 0; i < 10; i++) {
      print_particle(particles[i], i);
    }
    std::cout << std::endl;
  }

  return particles;
}

TreeAndParticles generate_simple_tree(int max_depth, bool verbose) {
  // Create a very simple test case with exactly 4 particles at known positions
  const int NUM_PARTICLES = 4;

  // Create particles array with predefined positions at different octants
  std::vector<particle_t> particles(NUM_PARTICLES);
  
  // Create simple test case with particles in different octants
  particles[0].pos[0] = pos_t(0.1f); particles[0].pos[1] = pos_t(0.1f); particles[0].pos[2] = pos_t(0.1f); // octant 0
  particles[1].pos[0] = pos_t(0.9f); particles[1].pos[1] = pos_t(0.1f); particles[1].pos[2] = pos_t(0.1f); // octant 1
  particles[2].pos[0] = pos_t(0.1f); particles[2].pos[1] = pos_t(0.9f); particles[2].pos[2] = pos_t(0.1f); // octant 2
  particles[3].pos[0] = pos_t(0.9f); particles[3].pos[1] = pos_t(0.9f); particles[3].pos[2] = pos_t(0.9f); // octant 7
  
  // Set mass and index for all particles
  for (int i = 0; i < NUM_PARTICLES; i++) {
    particles[i].mass = mass_t(0.25f); // 1.0 / 4 particles
    particles[i].idx = count_t(i);
  }

  if (verbose) {
    std::cout << "Created 4 particles for simple tree test" << std::endl;
  }

  // Sort particles by PH key
  particles = phsort(particles, max_depth, verbose);

  // Print first and last ph key
  if (verbose) {
    std::cout << "First ph key: " << particles[0].key << std::endl;
    std::cout << "Last ph key: " << particles[particles.size() - 1].key << std::endl;
  }

  std::vector<nodeleaf> tree = treecon_wrapper(particles, NUM_PARTICLES, verbose);

  return {tree, particles};
}

TreeAndParticles generate_max_depth_tree(int max_depth, bool verbose) {
  // Create particles designed to test maximum depth behavior
  const int NUM_PARTICLES = 10;

  // Create particles array with positions that will trigger max depth
  std::vector<particle_t> particles(NUM_PARTICLES);
  
  // Create particles with slightly different x coordinates but same y,z
  // This creates a scenario that pushes the tree to maximum depth
  double step = (1.0 / pow(2, max_depth));
  for (int i = 0; i < NUM_PARTICLES; i++) {
    particles[i].pos[0] = pos_t(0.1f + float(i * step / NUM_PARTICLES));
    particles[i].pos[1] = pos_t(0.1f);
    particles[i].pos[2] = pos_t(0.1f);
    particles[i].mass = mass_t(1.0f / NUM_PARTICLES);
    particles[i].idx = count_t(i);
  }

  if (verbose) {
    std::cout << "Created 10 particles for max depth tree test" << std::endl;
  }

  // Sort particles by PH key
  particles = phsort(particles, max_depth, verbose);

  // Print first and last ph key
  if (verbose) {
    std::cout << "First ph key: " << particles[0].key << std::endl;
    std::cout << "Last ph key: " << particles[particles.size() - 1].key << std::endl;
  }

  std::vector<nodeleaf> tree = treecon_wrapper(particles, NUM_PARTICLES, verbose);

  return {tree, particles};
}

TreeAndParticles generate_random_tree(int num_particles, int max_depth,
                                      bool verbose) {
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
  }

  // sort particles by PH key
  particles = phsort(particles, max_depth, verbose);

  // Set idx
  for (int i = 0; i < NUM_PARTICLES; i++) {
    particles[i].idx = i;
  }

  // print first and last ph key
  std::cout << "First ph key: " << particles[0].key << std::endl;
  std::cout << "Last ph key: " << particles[particles.size() - 1].key << std::endl;

  std::vector<nodeleaf> tree = treecon_wrapper(particles, NUM_PARTICLES, verbose);

  return {tree, particles};
}

TreeAndParticles generate_quasirandom_tree(int num_particles, int max_depth,
  bool verbose, double pos0[3]) {
const int NUM_PARTICLES = 8000;

  double g = 1.22074408460575947536;
  double a1 = 1.0/g;
  double a2 = 1.0/(g*g);
  double a3 = 1.0/(g*g*g);
std::vector<particle_t> particles(NUM_PARTICLES);

  double pos_dbl[3] = {pos0[0], pos0[1], pos0[2]};

  for (int i = 0; i < NUM_PARTICLES; i++) {
    particles[i].pos[0] = pos_t(pos_dbl[0]);
    particles[i].pos[1] = pos_t(pos_dbl[1]);
    particles[i].pos[2] = pos_t(pos_dbl[2]);

    particles[i].mass = 1.0f / NUM_PARTICLES;

    pos_dbl[0] = std::fmod(pos_dbl[0]+a1, 1.0);
    pos_dbl[1] = std::fmod(pos_dbl[1]+a2, 1.0);
    pos_dbl[2] = std::fmod(pos_dbl[2]+a3, 1.0);
  }

  // sort particles by PH key
  particles = phsort(particles, max_depth, verbose);

  // Set idx
  for (int i = 0; i < NUM_PARTICLES; i++) {
    particles[i].idx = i;
  }

  // print first and last ph key
  std::cout << "First ph key: " << particles[0].key << std::endl;
  std::cout << "Last ph key: " << particles[particles.size() - 1].key << std::endl;

  std::vector<nodeleaf> tree = treecon_wrapper(particles, NUM_PARTICLES, verbose);

  return {tree, particles};
}

#endif