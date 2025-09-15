#include "hls_math.h"
#include "treecon_config_hls.h"
#include "treecon_types_hls.h"
#include <cmath>
#include <cstring>

#ifndef __SYNTHESIS__
#include <iostream>
#endif

const int ARRAY_SIZE = 32;

struct nodepart_packet {
  nodeleaf node;
  particle_t particle;
};

void force_kernel(hls::stream<nodepart_packet> &nodepart_toforce, hls::stream<particle_t> &part_toread, pos_t theta) {
  // read in the node and particle
  while(1) {
    #pragma HLS PIPELINE II=1
    nodepart_packet nodepart = nodepart_toforce.read();
    nodeleaf node = nodepart.node;
    particle_t particle = nodepart.particle;

    part_toread.write(particle);
  }
}

void read_kernel(hls::stream<particle_t> &part_toread, hls::stream<nodepart_packet> &nodepart_toforce, nodeleaf *tree) {
  while(1) {
    #pragma HLS PIPELINE II=1
    particle_t particle = part_toread.read();
    nodeleaf node = tree[particle.next_tree_idx];
    nodepart_packet nodepart = {node, particle};
    nodepart_toforce.write(nodepart);
  }
}

void treewalk_kernel(hls::stream<particle_t> &particle_in,
                     hls::stream<particle_t> &particle_out, nodeleaf *tree,
                     pos_t theta, double G) {

#pragma HLS DATAFLOW

  hls::stream<nodepart_packet> nodepart_toforce("nodepart_toforce");
  hls::stream<particle_t> part_toread("part_toread");

  force_kernel(nodepart_toforce, part_toread, theta);
  read_kernel(part_toread, nodepart_toforce, tree);
}