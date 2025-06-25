#ifndef BHTREE_TYPES_HLS_H
#define BHTREE_TYPES_HLS_H

// Type definitions for Vitis HLS

#include "ap_int.h"
#include "ap_fixed.h"

typedef ap_uint<32> posint_t[3];
typedef ap_fixed<32,0> pos_t;  // 32-bit fixed point, 0 fractional bits
typedef ap_uint<30> phkey_t;
typedef ap_uint<4> level_t; // max depth is 10
typedef ap_fixed<32,0> mass_t; // mass type

typedef ap_uint<32> count_t;

#endif // BHTREE_TYPES_HLS_H