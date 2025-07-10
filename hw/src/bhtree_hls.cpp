#include "bhtree_config_hls.h"
#include "bhtree_types_hls.h"
#include <cstring>

const int NODE_BUFFER_DEPTH = 128;

// The small control struct to be passed over the synchronization stream.
struct token_t {
    ap_uint<1> buffer_idx; // The buffer that is full (0 or 1)
    ap_uint<16> num_nodes;  // How many nodes are in it
    bool is_last;          // Flag to signal termination
};

nodeleaf add_particle_to_node(const nodeleaf &input_node,
                              const particle_t &particle, bool set_idx) {
#pragma HLS INLINE
// #pragma HLS PIPELINE II=1
#pragma HLS INTERFACE ap_none port = particle

  nodeleaf output_node = input_node;

  output_node.pos[0] += particle.pos[0] * particle.mass;
  output_node.pos[1] += particle.pos[1] * particle.mass;
  output_node.pos[2] += particle.pos[2] * particle.mass;
  output_node.mass += particle.mass;
  output_node.num_particles++;

  if (set_idx) {
    // max particle count is 2^31, so -1 will never be a valid particle index
    output_node.start_idx = particle.idx;
  }

  return output_node;
}

void particle_processor(hls::stream<particle_t> &particle_stream,
                        hls::stream<token_t> &token_stream,
                        nodeleaf node_buffer_0[NODE_BUFFER_DEPTH],
                        nodeleaf node_buffer_1[NODE_BUFFER_DEPTH],
                        count_t num_particles) {
#pragma HLS INTERFACE ap_memory port=node_buffer_0
#pragma HLS INTERFACE ap_memory port=node_buffer_1

  particle_t p = particle_stream.read();

  int buf_idx = 0;
  int current_buf = 0;

  nodeleaf active_stack[MAX_DEPTH];
  int active_stack_num_nodes = 0;
  nodeleaf result_stack[MAX_DEPTH];
  int result_stack_num_nodes = 0;
#pragma HLS ARRAY_PARTITION variable = active_stack complete dim = 1
#pragma HLS ARRAY_PARTITION variable = result_stack complete dim = 1

INIT_STACK:
  for (int i = 0; i < MAX_DEPTH; i++) {
#pragma HLS UNROLL
    active_stack[i].key = p.key >> (3 * (MAX_DEPTH - (i + 1)));
    active_stack[i].level = i + 1;
    active_stack[i].pos[0] = p.pos[0];
    active_stack[i].pos[1] = p.pos[1];
    active_stack[i].pos[2] = p.pos[2];
    active_stack[i].mass = p.mass;
    active_stack[i].num_particles = 1;
    active_stack[i].start_idx = 0;
    active_stack[i].is_leaf = true;
    active_stack[i].is_last = false;
  }
  active_stack_num_nodes = MAX_DEPTH;

PROCESS_PARTICLES:
  for (int i = 1; i < num_particles; i++) {
#pragma HLS PIPELINE II = 1
    p = particle_stream.read();

    bool flush_mode = false;
    bool quiet_mode = false;

  PROCESS_LEVELS:
    for (int level = 1; level <= MAX_DEPTH; level++) {
#pragma HLS UNROLL

      // Calculate node key for this level
      phkey_t node_key = p.key >> (3 * (MAX_DEPTH - level));

      if (node_key == active_stack[level - 1].key && !flush_mode) {
        // Key matches and not in flush mode - just add particle
        active_stack[level - 1] =
            add_particle_to_node(active_stack[level - 1], p, false);
      } else {
        // Key mismatch or in flush mode - need to flush
        flush_mode = true;

        if (!quiet_mode) {
          // Add current node to output
          result_stack[result_stack_num_nodes] =
              active_stack[level - 1];
          result_stack_num_nodes++;
        }

        // Check if we should enter quiet mode
        if (active_stack[level - 1].num_particles <= NLEAF) {
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
        new_node.is_last = false;

        // Add particle to new node
        active_stack[level - 1] = add_particle_to_node(new_node, p, true);
      }
    } // end PROCESS_LEVELS

    // now we place the result stack into the node stream in reverse order
    for (int j = 0; j < MAX_DEPTH; j++) {
#pragma HLS UNROLL
      if (current_buf == 0) {
          node_buffer_0[buf_idx + j] = result_stack[j];
      } else {
          node_buffer_1[buf_idx + j] = result_stack[j];
      }
    }

    // we do a fixed size copy to make this hw friendly, but increment the index
    // only by the number of valid nodes.
    buf_idx += result_stack_num_nodes;

    if (NODE_BUFFER_DEPTH - buf_idx < MAX_DEPTH) {
      // next particle could in theory fill the buffer, so we need to flush it.
      token_stream.write(token_t{(unsigned int)current_buf, (unsigned short)buf_idx, false});
      current_buf = (current_buf + 1) % 2;
      buf_idx = 0;
    }
    result_stack_num_nodes = 0;
  } // end PROCESS_PARTICLES

  // write out the last stack. we just need to count the number of nodes in the
  // stack
  active_stack[0].is_last = true;

  bool quiet_mode = false;
  result_stack_num_nodes = 0;
  for (int level = 1; level <= MAX_DEPTH; level++) {
#pragma HLS UNROLL
    if (!quiet_mode) {
      result_stack[result_stack_num_nodes] =
          active_stack[level - 1];
      result_stack_num_nodes++;
    }
    if (active_stack[level - 1].num_particles <= NLEAF) {
      quiet_mode = true;
    }
  }

  // now we place the result stack into the node stream in reverse order
  for (int j = 0; j < MAX_DEPTH; j++) {
#pragma HLS UNROLL
    if (current_buf == 0) {
        node_buffer_0[buf_idx + j] = result_stack[j];
    } else {
        node_buffer_1[buf_idx + j] = result_stack[j];
    }
  }

  // flush bc last particle
  token_stream.write(token_t{(unsigned int)current_buf, (unsigned short)buf_idx, true});
}

// Write nodes to two independent AXI ports (tree0 / tree1) so that we can
// sustain two writes per cycle when the interface width is eventually bumped
// to 512-bit.  Nodes are written in an interleaved fashion: even index ➜ tree0,
// odd index ➜ tree1.  Each tree pointer is connected to a distinct AXI master
// interface bundle in the kernel top function.
void node_writer(hls::stream<token_t> &token_stream,
                 nodeleaf node_buffer_0[NODE_BUFFER_DEPTH],
                 nodeleaf node_buffer_1[NODE_BUFFER_DEPTH],
                 nodeleaf *tree0,
                 nodeleaf *tree1) {
  long long int idx = 0;
  long long int idx0 = 0; // write index for tree0
  long long int idx1 = 0; // write index for tree1

  NODE_WRITE_LOOP: while(1) {
    token_t t = token_stream.read();

    // for(int i = 0; i < t.num_nodes; i++) {
    //   switch(t.buffer_idx) {
    //     case 0:
    //       tree[idx++] = node_buffer_0[i];
    //       break;
    //     case 1:
    //       tree[idx++] = node_buffer_1[i];
    //       break;
    //   }
    // }

        // Use fixed loop bounds for better HLS synthesis
    // We stream out two nodes per iteration so that the scheduler can map
    // each write to a different AXI interface (tree0 / tree1) and therefore
    // maintain II=1.
    WRITE_NODES: for(int i = 0; i < NODE_BUFFER_DEPTH; i += 2) {
#pragma HLS PIPELINE II=1
#pragma HLS LOOP_TRIPCOUNT min=1 max=64 avg=32
      if (idx+i < t.num_nodes) {
        // write even element to tree0
        nodeleaf node_even = (t.buffer_idx == 0) ? node_buffer_0[i]
                                                : node_buffer_1[i];
        tree0[idx+i] = node_even;
      }

      if (idx+i + 1 < t.num_nodes) {
        // write odd element to tree1
        nodeleaf node_odd = (t.buffer_idx == 0) ? node_buffer_0[i + 1]
                                               : node_buffer_1[i + 1];
        tree1[idx+i+1] = node_odd;
      }
    }

    idx += t.num_nodes; // keep global count unchanged for now (if used later)

    if(t.is_last) break;
  } 
}

void create_bhtree_kernel(hls::stream<particle_t> &particle_stream,
                          nodeleaf *tree0,
                          nodeleaf *tree1,
                          count_t num_particles) {
#pragma HLS DATAFLOW

#pragma HLS INTERFACE m_axi port=tree0 offset=slave bundle=gmem0 depth=MAX_NODES max_widen_bitwidth=1024 max_write_burst_length=256
#pragma HLS INTERFACE m_axi port=tree1 offset=slave bundle=gmem1 depth=MAX_NODES max_widen_bitwidth=1024 max_write_burst_length=256

//   static nodeleaf node_buffer[2][NODE_BUFFER_DEPTH];
// #pragma HLS ARRAY_PARTITION variable=node_buffer complete dim=1
// #pragma HLS ARRAY_RESHAPE variable=node_buffer factor=MAX_DEPTH dim=2
// #pragma HLS BIND_STORAGE variable=node_buffer type=RAM_2P impl=BRAM

  nodeleaf node_buffer_0[NODE_BUFFER_DEPTH];
  nodeleaf node_buffer_1[NODE_BUFFER_DEPTH];
#pragma HLS ARRAY_PARTITION variable=node_buffer_0 cyclic factor=10 dim=1
#pragma HLS ARRAY_PARTITION variable=node_buffer_1 cyclic factor=10 dim=1
#pragma HLS BIND_STORAGE variable=node_buffer_0 type=RAM_2P impl=BRAM
#pragma HLS BIND_STORAGE variable=node_buffer_1 type=RAM_2P impl=BRAM

  hls::stream<token_t> token_stream("token_stream");
#pragma HLS STREAM variable=token_stream depth=8

  particle_processor(particle_stream, token_stream, node_buffer_0, node_buffer_1, num_particles);

  node_writer(token_stream, node_buffer_0, node_buffer_1, tree0, tree1);
}