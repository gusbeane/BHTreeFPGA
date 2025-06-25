#ifndef BHTREE_TYPES_HLS_H
#define BHTREE_TYPES_HLS_H

// Type definitions for Vitis HLS

#include "ap_int.h"

typedef ap_uint<32> posint_t;
typedef ap_fixed<32,0> pos_t;  // 32-bit fixed point, 0 fractional bits
typedef ap_uint<30> phkey_t;

#endif // BHTREE_TYPES_HLS_H