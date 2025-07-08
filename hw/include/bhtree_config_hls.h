#ifndef BHTREE_CONFIG_HLS_H
#define BHTREE_CONFIG_HLS_H

#include "bhtree_types_hls.h"
#include "hls_stream.h"

const unsigned int MAX_DEPTH = 10;
const unsigned int MAX_PARTICLES = (1U << 31);
const unsigned int MAX_NODES = (2 * MAX_PARTICLES);
const unsigned int NLEAF = 1;

const unsigned int NODE_WRITE_BURST_SIZE = 16;

struct particle_t {
    pos_t pos[3];
    mass_t mass;
    count_t idx;
    phkey_t key;
};

struct nodeleaf {
    phkey_t key; // 30 bits
    level_t level; // 32 bits
    
    pos_t pos[3]; // center of mass if node, position of first particle if leaf // 96 bits

    mass_t mass; // 32 bits

    count_t start_idx; // 32 bits
    count_t num_particles; // 32 bits

    bool is_leaf; // 1 bit

    // total is 255 bits, pad to 256 bits
    ap_uint<1> padding;
};

struct nodeleaf_stack {
    nodeleaf nodes[MAX_DEPTH];
    level_t num_nodes;
};

#endif