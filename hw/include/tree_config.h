#ifndef TREECON_CONFIG_HLS_H
#define TREECON_CONFIG_HLS_H

#include "tree_types.h"
#include "hls_stream.h"

const unsigned int MAX_DEPTH = 10;
const unsigned int MAX_PARTICLES = (1U << 31);
const unsigned int MAX_NODES = (2 * MAX_PARTICLES);
const unsigned int NLEAF = 1;

const unsigned int NODE_WRITE_BURST_SIZE = 16;

struct particle_t {
    pos_t pos[3];
    acc_t acc[3];
    mass_t mass;
    count_t idx;
    phkey_t key;
    count_t next_tree_idx;
    double h;
    bool valid;
} __attribute__((aligned(64))); // padded to 64 bytes

struct nodeleaf {
    phkey_t key; // 30 bits
    level_t level; // 32 bits
    
    pos_t pos[3]; // center of mass if node, position of first particle if leaf // 96 bits

    mass_t mass; // 32 bits

    count_t start_idx; // 32 bits
    count_t num_particles; // 32 bits
    count_t next_sibling; // 32 bits

    ap_uint<1> is_leaf; // 1 bit
    ap_uint<1> is_last; // 1 bit

    count_t target_cell; // 32 bits
    bool valid; // 1 bit

    // min and max softenings
    float hmin;
    float hmax;

} __attribute__((aligned(64))); // padded to 64 bytes

#endif