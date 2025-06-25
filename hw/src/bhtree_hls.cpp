#include "bhtree_types_hls.h"
#include "bhtree_config_hls.h"

nodeleaf add_particle_to_node(const nodeleaf& input_node, 
                              const pos_t pos_x, const pos_t pos_y, const pos_t pos_z,
                              const mass_t mass, const count_t idx) {
    // #pragma HLS INLINE
    #pragma HLS PIPELINE II=1
    
    nodeleaf output_node = input_node;
    
    output_node.pos[0] += pos_x * mass;
    output_node.pos[1] += pos_y * mass;
    output_node.pos[2] += pos_z * mass;
    output_node.mass += mass;
    output_node.num_particles++;

    if(idx != -1) {
        // max particle count is 2^31, so -1 will never be a valid particle index
        output_node.start_idx = idx;
    }

    return output_node;
}