#include "bhtree_types_hls.h"
#include "bhtree_config_hls.h"

nodeleaf add_particle_to_node(const nodeleaf& input_node, const particle_t& particle, bool set_idx) {
    #pragma HLS INLINE
    // #pragma HLS PIPELINE II=1
    #pragma HLS INTERFACE ap_none port=particle
    
    nodeleaf output_node = input_node;
    
    output_node.pos[0] += particle.pos[0] * particle.mass;
    output_node.pos[1] += particle.pos[1] * particle.mass;
    output_node.pos[2] += particle.pos[2] * particle.mass;
    output_node.mass += particle.mass;
    output_node.num_particles++;

    if(set_idx) {
        // max particle count is 2^31, so -1 will never be a valid particle index
        output_node.start_idx = particle.idx;
    }

    return output_node;
}

nodeleaf_stack add_particle_to_tree(
    const nodeleaf_stack& input_stack,  // Current node stack
    const particle_t& particle,
    hls::stream<nodeleaf_stack>& node_stream
) {
    // #pragma HLS INTERFACE ap_none port=particle
    // #pragma HLS INTERFACE ap_none port=input_stack
    #pragma HLS ARRAY_PARTITION variable=input_stack.nodes complete dim=1
    #pragma HLS PIPELINE II=1
    
    nodeleaf_stack result;
    #pragma HLS ARRAY_PARTITION variable=result.nodes complete dim=1
    
    // Initialize output
    result.num_nodes = 0;
    
    // Copy input stack to working stack
    nodeleaf_stack temp_stack;
    #pragma HLS ARRAY_PARTITION variable=temp_stack complete dim=1
    
    COPY_STACK: for(int i = 0; i < MAX_DEPTH; i++) {
        #pragma HLS UNROLL
        temp_stack.nodes[i] = input_stack.nodes[i];
    }
    
    bool flush_mode = false;
    bool quiet_mode = false;
    
    PROCESS_LEVELS: for(int level = 1; level <= MAX_DEPTH; level++) {
        #pragma HLS UNROLL
        
        // Calculate node key for this level
        phkey_t node_key = particle.key >> (3 * (MAX_DEPTH - level));
        
        if(node_key == temp_stack.nodes[level-1].key && !flush_mode) {
            // Key matches and not in flush mode - just add particle
            temp_stack.nodes[level-1] = add_particle_to_node(
                temp_stack.nodes[level-1], particle, false
            );
        } else {
            // Key mismatch or in flush mode - need to flush
            flush_mode = true;
            
            if(!quiet_mode) {
                // Add current node to output
                result.nodes[result.num_nodes] = temp_stack.nodes[level-1];
                result.num_nodes++;
            }

            // Check if we should enter quiet mode
            if(temp_stack.nodes[level-1].num_particles <= NLEAF) {
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
            temp_stack.nodes[level-1] = add_particle_to_node(
                new_node, particle, true
            );
        }
    }
    
    // Copy working stack to output stack
    COPY_OUTPUT_STACK: for(int i = 0; i < MAX_DEPTH; i++) {
        #pragma HLS UNROLL
        result.nodes[i] = temp_stack.nodes[i];
    }
    
    node_stream.write(result);

    return result;
}

