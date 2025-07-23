#include "bhtree_config_hls.h"
#include "bhtree_types_hls.h"
#include <cstring>

#ifndef __SYNTHESIS__
#include <iostream>
#endif

nodeleaf add_particle_to_node(const nodeleaf &input_node,
                              const particle_t &particle, bool set_idx) {
#pragma HLS INLINE
#pragma HLS INTERFACE ap_none port = particle

  nodeleaf output_node = input_node;

  output_node.pos[0] += particle.pos[0] * particle.mass;
  output_node.pos[1] += particle.pos[1] * particle.mass;
  output_node.pos[2] += particle.pos[2] * particle.mass;
  output_node.mass += particle.mass;
  output_node.num_particles++;

  if (set_idx) {
    output_node.start_idx = particle.idx;
  }

  return output_node;
}

nodeleaf generate_empty_node(phkey_t node_key, level_t level) {
#pragma HLS INLINE

  nodeleaf node;
  node.key = node_key;
  node.level = level;
  node.pos[0] = 0;
  node.pos[1] = 0;
  node.pos[2] = 0;
  node.mass = 0;
  node.num_particles = 0;
  node.start_idx = 0;
  node.is_leaf = true;
  node.is_last = false;

  return node;
}

void particle_processor(const particle_t *particles,
                        hls::stream<nodeleaf> &node_stream,
                        count_t num_particles) {

  particle_t p = particles[0];

  nodeleaf active_stack[MAX_DEPTH];
  nodeleaf result_stack[MAX_DEPTH];
  int result_stack_num_nodes = 0;
#pragma HLS ARRAY_PARTITION variable = active_stack complete dim = 1
#pragma HLS ARRAY_PARTITION variable = result_stack complete dim = 1

INIT_STACK:
  for (int i = 0; i < MAX_DEPTH; i++) {
#pragma HLS UNROLL
    nodeleaf new_node = generate_empty_node(p.key >> (3 * (MAX_DEPTH - (i + 1))), i + 1);
    active_stack[i] = add_particle_to_node(new_node, p, true);
  }

PROCESS_PARTICLES:
  for (int i = 1; i < num_particles; i++) {
    p = particles[i];

    bool flush_mode = false;
    bool quiet_mode = false;
    result_stack_num_nodes = 0;

  PROCESS_LEVELS:
    for (int level = 1; level <= MAX_DEPTH; level++) {
#pragma HLS UNROLL

      phkey_t node_key = p.key >> (3 * (MAX_DEPTH - level));

      if (node_key == active_stack[level - 1].key && !flush_mode) {
        active_stack[level - 1] =
            add_particle_to_node(active_stack[level - 1], p, false);
      } else {
        flush_mode = true;

        if (!quiet_mode) {
          result_stack[result_stack_num_nodes] = active_stack[level - 1];
          result_stack_num_nodes++;
        }

        if (active_stack[level - 1].num_particles <= NLEAF) {
          quiet_mode = true;
        }

        nodeleaf new_node = generate_empty_node(node_key, level);
        active_stack[level - 1] = add_particle_to_node(new_node, p, true);
      }
    }

    // Write result nodes to stream one at a time
    for (int j = 0; j < result_stack_num_nodes; j++) {
      node_stream.write(result_stack[result_stack_num_nodes - j - 1]);
    }
  }

  // Write final stack
  active_stack[0].is_last = true;
  
  bool quiet_mode = false;
  result_stack_num_nodes = 0;
  for (int level = 1; level <= MAX_DEPTH; level++) {
#pragma HLS UNROLL
    if (!quiet_mode) {
      result_stack[result_stack_num_nodes] = active_stack[level - 1];
      result_stack_num_nodes++;
    }
    if (active_stack[level - 1].num_particles <= NLEAF) {
      quiet_mode = true;
    }
  }

  for (int j = 0; j < result_stack_num_nodes; j++) {
    node_stream.write(result_stack[result_stack_num_nodes - j - 1]);
  }
}

void node_writer(hls::stream<nodeleaf> &node_stream,
                 ap_uint<512> *tree) {
  
  long long int idx = 0;

  NODE_WRITE_LOOP: while(1) {
#pragma HLS PIPELINE II=1
    
    nodeleaf node = node_stream.read();

    // Now we close the node
    node.pos[0] /= node.mass;
    node.pos[1] /= node.mass;
    node.pos[2] /= node.mass;

    // Mark the node as a leaf if it has less than NLEAF particles
    node.is_leaf = (node.num_particles <= NLEAF) || node.level == MAX_DEPTH;
#ifndef __SYNTHESIS__
    std::cout << "Node " << idx << " is_leaf: " << node.is_leaf << std::endl;
#endif
    // Reinterpret nodeleaf as ap_uint<512> for efficient AXI write
    ap_uint<512> node_bits;
    node_bits = *reinterpret_cast<ap_uint<512>*>(&node);
    
    tree[idx] = node_bits;
    idx++;
    
    if (node.is_last) break;
  }
}

void create_bhtree_kernel(const particle_t *particles, ap_uint<512> *tree,
                          count_t num_particles) {
#pragma HLS INTERFACE m_axi port = particles bundle = gmem0 depth = 1024
#pragma HLS INTERFACE m_axi port = tree offset = slave bundle = gmem1 depth = 1024 max_widen_bitwidth = 512 max_write_burst_length = 64 num_write_outstanding = 16 latency = 64
#pragma HLS INTERFACE s_axilite port=particles bundle=control
#pragma HLS INTERFACE s_axilite port=num_particles bundle=control
#pragma HLS INTERFACE s_axilite port=tree bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control

#pragma HLS DATAFLOW

  hls::stream<nodeleaf> node_stream("node_stream");
#pragma HLS STREAM variable=node_stream depth=64

  particle_processor(particles, node_stream, num_particles);
  node_writer(node_stream, tree);
}