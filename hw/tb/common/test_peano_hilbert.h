#ifndef TEST_PEANO_HILBERT_H
#define TEST_PEANO_HILBERT_H

#include <iostream>
#include <vector>
#include <array>

// Include HLS headers
#include "peano_hilbert.h"

bool test_peano_hilbert_key() {
    // Create Peano-Hilbert generator with depth 10
    PeanoHilbert ph(10);
        
    // Single position key generation
    uint64_t key = ph.generate_key(100, 200, 300);
    std::cout << "Key for position (100, 200, 300): 0x" << std::hex << key << std::dec << std::endl;
    
    // Multiple positions
    std::vector<std::array<uint32_t, 3>> positions = {
        {0, 0, 0},
        {1000, 2000, 3000},
        {0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF},
        {
            static_cast<uint32_t>(0.4 * static_cast<double>(UINT32_MAX)),
            static_cast<uint32_t>(0.8 * static_cast<double>(UINT32_MAX)),
            static_cast<uint32_t>(0.3 * static_cast<double>(UINT32_MAX))
        },
        {1710436918, 3435321836, 4283210188}
    };
    
    std::vector<uint64_t> keys = ph.generate_keys(positions);
    
    // Ensure keys match golden values
    std::vector<uint64_t> golden_keys = {
        0x0,
        0x0,
        0x29249249, 
        0xdf03f03, 
        0x310f1b87,
    };
    
    bool all_match = true;
    
    for (size_t i = 0; i < keys.size(); ++i) {
        if (keys[i] != golden_keys[i]) {
            std::cout << "PH key mismatch at position " << i << ": " << keys[i] << " != " << golden_keys[i] << std::endl;
            all_match = false;
        }
    }
    
    if(all_match) {
        std::cout << "All keys match golden values" << std::endl;
    } else {
        std::cout << "Keys do not match golden values" << std::endl;
    }
    
    return all_match;
}

#endif // TEST_PEANO_HILBERT_H