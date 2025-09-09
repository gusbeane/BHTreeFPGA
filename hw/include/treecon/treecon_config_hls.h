#ifndef TREECON_CONFIG_HLS_H
#define TREECON_CONFIG_HLS_H

#include "treecon_types_hls.h"
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

    ap_uint<1> is_leaf; // 1 bit
    ap_uint<1> is_last; // 1 bit

    // total is 256 bits, no need to pad to 256 bits
    // ap_uint<1> padding;
} __attribute__((aligned(64)));

// union nodeleaf {
//     struct {
//         phkey_t key;
//         level_t level;
//         pos_t pos[3];
//         mass_t mass;
//         count_t start_idx;
//         count_t num_particles;
//         ap_uint<1> is_leaf;
//         ap_uint<1> is_last;
//     };
//     ap_uint<8> raw_data[64]; // Forces 512-bit (64-byte) size
// };

#endif