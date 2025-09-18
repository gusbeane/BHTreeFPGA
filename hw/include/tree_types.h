#ifndef TREECON_TYPES_HLS_H
#define TREECON_TYPES_HLS_H

// Type definitions for Vitis HLS

#include "ap_int.h"
#include "ap_fixed.h"

typedef ap_ufixed<32,0> pos_t;  // 32-bit fixed point, 0 fractional bits
typedef double acc_t; // is double necessary?
typedef ap_uint<30> phkey_t;
typedef int level_t; // max depth is 10, this allows for up to 64 levels
typedef ap_ufixed<32,0, AP_RND, AP_SAT> mass_t; // mass type

typedef int count_t;

#endif // TREECON_TYPES_HLS_H