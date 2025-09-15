#ifndef TREECON_H
#define TREECON_H

#include "tree_config.h"
#include "tree_types.h"

void create_bhtree_kernel(const particle_t *particles, ap_uint<512> *tree,
    count_t num_particles);

#endif