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
  ap_uint<1> is_particle;
  ap_uint<512> nodeleaf_or_particle;
};

template <int cell_index>
void cell(pos_t theta, double G, hls::stream<nodepart_packet> &data_in,
          hls::stream<nodepart_packet> &data_out) {

  //   const pos_t node_sizes[MAX_DEPTH] = {
  //   0.5,      0.25,      0.125,      0.0625,      0.03125,
  //   0.015625, 0.0078125, 0.00390625, 0.001953125, 0.0009765625};

  static nodeleaf local_node;
  static pos_t local_node_size_theta_sq;
  static int has_node = 0;

  nodepart_packet pkt = data_in.read();

  if (pkt.is_particle && has_node) {
    particle_t particle =
        *reinterpret_cast<particle_t *>(&pkt.nodeleaf_or_particle);

    if (particle.next_cell > cell_index) {
      data_out.write(pkt);
    } else if (particle.next_cell == cell_index) {
      // now we do MAC test
      pos_t dx = particle.pos[0] - local_node.pos[0];
      pos_t dy = particle.pos[1] - local_node.pos[1];
      pos_t dz = particle.pos[2] - local_node.pos[2];
      pos_t rsq = dx * dx + dy * dy + dz * dz;
      if (rsq > local_node_size_theta_sq || local_node.is_leaf) {
        // we accept the multipole criterion or we are at a leaf, so we can
        // accumulate
        pos_t r = hls::sqrtf(rsq);
        pos_t rcubed = r * rsq;
        double rcubed_float = (double)(rcubed);
        double dx_float = (double)(dx);
        double dy_float = (double)(dy);
        double dz_float = (double)(dz);
        double local_node_mass = (double)local_node.mass;

        particle.acc[0] += G * local_node_mass * (dx_float) / (rcubed_float);
        particle.acc[1] += G * local_node_mass * (dy_float) / (rcubed_float);
        particle.acc[2] += G * local_node_mass * (dz_float) / (rcubed_float);

        particle.next_cell += local_node.next_sibling;

        nodepart_packet pkt_out;
        pkt_out.is_particle = 1;
        pkt_out.nodeleaf_or_particle =
            *reinterpret_cast<ap_uint<512> *>(&particle);
        data_out.write(pkt_out);
      } else {
        // we do not accept the multipole criterion, so we just pass the
        // particle along, incrementing its next cell
        particle.next_cell++;
        nodepart_packet pkt_out;
        pkt_out.is_particle = 1;
        pkt_out.nodeleaf_or_particle =
            *reinterpret_cast<ap_uint<512> *>(&particle);
        data_out.write(pkt_out);
      }

    } else {
      // this should not happen, flush the particle out with an invalid flag
      particle.valid = 0;
      particle.next_cell = -1;
      nodepart_packet pkt_out;
      pkt_out.is_particle = 1;
      pkt_out.nodeleaf_or_particle =
          *reinterpret_cast<ap_uint<512> *>(&particle);
      data_out.write(pkt_out);
    }
  } else {
    nodeleaf node = *reinterpret_cast<nodeleaf *>(&pkt.nodeleaf_or_particle);
    if (node.target_cell == cell_index) {
      local_node = node;
      //   local_node_size_theta_sq = node_sizes[node.level - 1] *
      //  node_sizes[node.level - 1] / (theta * theta);
      local_node_size_theta_sq = pos_t(1.0) / (theta * theta);
      has_node = 1;
    } else if (node.target_cell > cell_index) {
      data_out.write(pkt);
    } else if (node.target_cell < cell_index) {
      // this should not happen, flush the node out with an invalid flag
      node.valid = 0;
      node.target_cell = -1;
      data_out.write(pkt);
    }
  }
}

void systolic_array(hls::stream<nodepart_packet> &data_in,
                    hls::stream<nodepart_packet> &data_out, pos_t theta,
                    double G) {

#pragma HLS DATAFLOW

  hls::stream<nodepart_packet> data_streams[ARRAY_SIZE - 1];

  cell<0>(theta, G, data_in, data_streams[0]);
  cell<1>(theta, G, data_streams[0], data_streams[1]);
  cell<2>(theta, G, data_streams[1], data_streams[2]);
  cell<3>(theta, G, data_streams[2], data_streams[3]);
  cell<4>(theta, G, data_streams[3], data_streams[4]);
  cell<5>(theta, G, data_streams[4], data_streams[5]);
  cell<6>(theta, G, data_streams[5], data_streams[6]);
  cell<7>(theta, G, data_streams[6], data_streams[7]);
  cell<8>(theta, G, data_streams[7], data_streams[8]);
  cell<9>(theta, G, data_streams[8], data_streams[9]);
  cell<10>(theta, G, data_streams[9], data_streams[10]);
  cell<11>(theta, G, data_streams[10], data_streams[11]);
  cell<12>(theta, G, data_streams[11], data_streams[12]);
  cell<13>(theta, G, data_streams[12], data_streams[13]);
  cell<14>(theta, G, data_streams[13], data_streams[14]);
  cell<15>(theta, G, data_streams[14], data_streams[15]);
  cell<16>(theta, G, data_streams[15], data_streams[16]);
  cell<17>(theta, G, data_streams[16], data_streams[17]);
  cell<18>(theta, G, data_streams[17], data_streams[18]);
  cell<19>(theta, G, data_streams[18], data_streams[19]);
  cell<20>(theta, G, data_streams[19], data_streams[20]);
  cell<21>(theta, G, data_streams[20], data_streams[21]);
  cell<22>(theta, G, data_streams[21], data_streams[22]);
  cell<23>(theta, G, data_streams[22], data_streams[23]);
  cell<24>(theta, G, data_streams[23], data_streams[24]);
  cell<25>(theta, G, data_streams[24], data_streams[25]);
  cell<26>(theta, G, data_streams[25], data_streams[26]);
  cell<27>(theta, G, data_streams[26], data_streams[27]);
  cell<28>(theta, G, data_streams[27], data_streams[28]);
  cell<29>(theta, G, data_streams[28], data_streams[29]);
  cell<30>(theta, G, data_streams[29], data_streams[30]);
  cell<31>(theta, G, data_streams[30], data_out);
}