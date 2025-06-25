#ifndef BHTREE_CONFIG_HLS_H
#define BHTREE_CONFIG_HLS_H

#include "bhtree_types_hls.h"

#define MAX_DEPTH 10
#define MAX_PARTICLES (1U << 31)
#define MAX_NODES (2 * MAX_PARTICLES)
#define NLEAF 1

struct particle_t {
    pos_t pos[3];
    mass_t mass;
    count_t idx;
    phkey_t key;
};

struct nodeleaf {
    phkey_t key;
    level_t level;
    
    pos_t pos[3]; // center of mass if node, position of first particle if leaf

    mass_t mass;

    count_t start_idx;
    count_t num_particles;

    bool is_leaf;
};

#endif