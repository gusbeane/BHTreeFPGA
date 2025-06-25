#include "bhtree_types_hls.h"
#include "bhtree_config_hls.h"

nodeleaf add_particle_to_node(const nodeleaf& input_node, 
                              const pos_t pos_x, const pos_t pos_y, const pos_t pos_z,
                              const mass_t mass, const count_t idx) {
    #pragma HLS INLINE
    // #pragma HLS PIPELINE II=1
    
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

struct tree_output {
    nodeleaf output_nodes[MAX_DEPTH];     // Nodes to write to DRAM tree
    nodeleaf new_stack[MAX_DEPTH];        // Updated stack for next particle
    count_t num_output_nodes;             // Number of valid nodes in output_nodes
};

tree_output add_particle_to_tree(
    const nodeleaf input_stack[MAX_DEPTH],  // Current node stack
    const pos_t pos_x, const pos_t pos_y, const pos_t pos_z,  // Particle position
    const mass_t mass,                      // Particle mass
    const phkey_t particle_key,             // Particle's Peano-Hilbert key
    const count_t particle_idx              // Particle index
) {
    #pragma HLS INTERFACE ap_none port=input_stack
    #pragma HLS ARRAY_PARTITION variable=input_stack complete dim=1
    #pragma HLS PIPELINE II=1
    
    tree_output result;
    #pragma HLS ARRAY_PARTITION variable=result.output_nodes complete dim=1
    #pragma HLS ARRAY_PARTITION variable=result.new_stack complete dim=1
    
    // Initialize output
    result.num_output_nodes = 0;
    
    // Copy input stack to working stack
    nodeleaf temp_stack[MAX_DEPTH];
    #pragma HLS ARRAY_PARTITION variable=temp_stack complete dim=1
    
    COPY_STACK: for(int i = 0; i < MAX_DEPTH; i++) {
        #pragma HLS UNROLL
        temp_stack[i] = input_stack[i];
    }
    
    bool flush_mode = false;
    bool quiet_mode = false;
    
    PROCESS_LEVELS: for(int level = 1; level <= MAX_DEPTH; level++) {
        #pragma HLS UNROLL
        
        // Calculate node key for this level
        phkey_t node_key = particle_key >> (3 * (MAX_DEPTH - level));
        
        if(node_key == temp_stack[level-1].key && !flush_mode) {
            // Key matches and not in flush mode - just add particle
            temp_stack[level-1] = add_particle_to_node(
                temp_stack[level-1], pos_x, pos_y, pos_z, mass, particle_idx
            );
        } else {
            // Key mismatch or in flush mode - need to flush
            flush_mode = true;
            
            if(!quiet_mode) {
                // Add current node to output
                result.output_nodes[result.num_output_nodes] = temp_stack[level-1];
                result.num_output_nodes++;
            }

            // Check if we should enter quiet mode
            if(temp_stack[level-1].num_particles <= NLEAF) {
                quiet_mode = true;
            }
            
            // Create new node for this level
            nodeleaf new_node;
            new_node.key = node_key;
            new_node.level = level;
            new_node.pos[0] = 0;
            new_node.pos[1] = 0;
            new_node.pos[2] = 0;
            new_node.mass = 0;
            new_node.num_particles = 0;
            new_node.start_idx = 0;
            new_node.is_leaf = true;
            
            // Add particle to new node
            temp_stack[level-1] = add_particle_to_node(
                new_node, pos_x, pos_y, pos_z, mass, particle_idx
            );
        }
    }
    
    // Copy working stack to output stack
    COPY_OUTPUT_STACK: for(int i = 0; i < MAX_DEPTH; i++) {
        #pragma HLS UNROLL
        result.new_stack[i] = temp_stack[i];
    }
    
    return result;
}