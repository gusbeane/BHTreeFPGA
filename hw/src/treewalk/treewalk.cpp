#include "hls_math.h"
#include "tree_config.h"
#include "tree_types.h"
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

particle_t particle_acc(particle_t particle, nodeleaf node, pos_t dx, pos_t dy,
                        pos_t dz, pos_t r2) {
#pragma HLS INLINE
  // compute the force
  double r2_double = (double)r2;
  double rinv_double = hls::rsqrt(r2_double);
  double rinv3_double = rinv_double * rinv_double * rinv_double;

  double force = (double)(node.mass * particle.mass);
  force *= rinv3_double;

  particle.acc[0] += force * (double)(dx);
  particle.acc[1] += force * (double)(dy);
  particle.acc[2] += force * (double)(dz);

  return particle;
}

void force_kernel(hls::stream<nodepart_packet> &nodepart_toforce,
                  hls::stream<particle_t> &part_toread, hls::stream<particle_t> &part_towrite, pos_t inv_thetasq) {
  // read in the node and particle
  static const pos_t node_sizes_sq[MAX_DEPTH] = {
      1 / (1 << 2),  1 / (1 << 4),  1 / (1 << 6),  1 / (1 << 8),
      1 / (1 << 10), 1 / (1 << 12), 1 / (1 << 14), 1 / (1 << 16),
      1 / (1 << 18), 1 / (1 << 20)};

  while(1) {
    #pragma HLS PIPELINE II=1
    nodepart_packet nodepart = nodepart_toforce.read();
    nodeleaf node = nodepart.node;
    particle_t particle = nodepart.particle;

    // do distance calculation
    pos_t dx = particle.pos[0] - node.pos[0];
    pos_t dy = particle.pos[1] - node.pos[1];
    pos_t dz = particle.pos[2] - node.pos[2];
    pos_t r2 = dx * dx + dy * dy + dz * dz;

    // do opening criterion
    pos_t nodesize_thetasq = node_sizes_sq[node.level] * inv_thetasq;
    if (r2 < nodesize_thetasq || node.is_leaf) {
      // we need to open this node
      particle.next_tree_idx += 1;
    } else {
      // we accept the opening criterion, or we are at a leaf
      particle = particle_acc(particle, node, dx, dy, dz, r2);
      particle.next_tree_idx += node.next_sibling;
    }

    if (node.next_sibling == -1u) {
      part_towrite.write(particle);
    } else {
      part_toread.write(particle);
    }
    continue;
  }
}

void write_kernel(hls::stream<particle_t> &part_towrite, particle_t *particle_out) {
  long idx = 0;
  while(1) {
    #pragma HLS PIPELINE II=1
    particle_t particle = part_towrite.read();
    particle_out[idx] = particle;
    idx++;
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
                     particle_t *particle_out, nodeleaf *tree,
                     pos_t theta) {

#pragma HLS DATAFLOW

  hls::stream<nodepart_packet> nodepart_toforce("nodepart_toforce");
  hls::stream<particle_t> part_toread("part_toread");
  hls::stream<particle_t> part_towrite("part_towrite");

  force_kernel(nodepart_toforce, part_toread, part_towrite, theta);
  read_kernel(part_toread, nodepart_toforce, tree);
  write_kernel(part_towrite, particle_out);
}