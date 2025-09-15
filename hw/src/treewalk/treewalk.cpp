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

void force_kernel(hls::stream<nodepart_packet> &data_in,
                    hls::stream<nodepart_packet> &data_out, pos_t theta,
                    double G) {

#pragma HLS DATAFLOW

  
}