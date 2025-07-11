#ifndef BHTREE_TYPES_HLS_H
#define BHTREE_TYPES_HLS_H

// Type definitions for Vitis HLS

#include "ap_int.h"
#include "ap_fixed.h"

typedef ap_uint<32> posint_t[3];
typedef ap_ufixed<32,0> pos_t;  // 32-bit fixed point, 0 fractional bits
typedef ap_uint<30> phkey_t;
typedef int level_t; // max depth is 10, this allows for up to 64 levels
typedef ap_ufixed<32,0> mass_t; // mass type

typedef int count_t;

#endif // BHTREE_TYPES_HLS_H