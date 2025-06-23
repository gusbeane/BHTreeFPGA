#pragma once

#include <vector>
#include <cstdint>
#include <stdexcept>
#include "bhtree_config.h"
#include "bhtree_types.h"

class PeanoHilbert {
private:
    int dimensions;
    int max_depth;
    
    // Rotation table for 3D space (48 x 8 matrix)
    static const uint8_t rottable3[48][8];
    
    // Subpixel lookup table for 3D space (48 x 8 matrix)
    static const uint8_t subpix3[48][8];

public:
    PeanoHilbert(int dimensions = 3, int max_depth = MAX_DEPTH);
    
    // Generate Peano-Hilbert key for a single position (32-bit)
    uint32_t generate_key(uint32_t x_in, uint32_t y_in, uint32_t z_in, int depth = -1);
    
    // Generate Peano-Hilbert key for a single position (64-bit)
    uint64_t generate_key_64(uint32_t x_in, uint32_t y_in, uint32_t z_in, int depth = -1);
    
    // Generate keys for multiple positions (32-bit)
    std::vector<uint32_t> generate_keys(const PositionVector& positions, int depth = -1);
    
    // Generate keys for multiple positions (64-bit)
    std::vector<uint64_t> generate_keys_64(const PositionVector& positions, int depth = -1);
}; 