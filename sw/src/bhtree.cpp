#include <iostream>
#include "../include/peano_hilbert.h"
#include "../include/bhtree_types.h"

int main(){
    // Test the PeanoHilbert class
    PeanoHilbert ph;
    
    // Test single key generation
    uint32_t key1 = ph.generate_key(100, 200, 300);
    std::cout << "Generated key for (100, 200, 300): " << key1 << std::endl;
    
    // Test multiple key generation
    PositionVector positions;
    Position3D pos1 = {UINT32_MAX/4, UINT32_MAX/4, UINT32_MAX/4};
    Position3D pos2 = {3*UINT32_MAX/4, UINT32_MAX/4, UINT32_MAX/4};
    Position3D pos3 = {UINT32_MAX/4, 3*UINT32_MAX/4, UINT32_MAX/4};
    Position3D pos4 = {3*UINT32_MAX/4, 3*UINT32_MAX/4, UINT32_MAX/4};
    Position3D pos5 = {UINT32_MAX/4, UINT32_MAX/4, 3*UINT32_MAX/4};
    Position3D pos6 = {3*UINT32_MAX/4, UINT32_MAX/4, 3*UINT32_MAX/4};
    Position3D pos7 = {UINT32_MAX/4, 3*UINT32_MAX/4, 3*UINT32_MAX/4};
    Position3D pos8 = {3*UINT32_MAX/4, 3*UINT32_MAX/4, 3*UINT32_MAX/4};
    
    positions.push_back(pos1);
    positions.push_back(pos2);
    positions.push_back(pos3);
    positions.push_back(pos4);
    positions.push_back(pos5);
    positions.push_back(pos6);
    positions.push_back(pos7);
    positions.push_back(pos8);
    
    auto keys = ph.generate_keys(positions);
    std::cout << "Generated keys for multiple positions:" << std::endl;
    for (size_t i = 0; i < keys.size(); i++) {
        std::cout << "Position (" << positions[i][0] << ", " << positions[i][1] << ", " << positions[i][2] 
                  << ") -> Key: " << keys[i] << std::endl;
    }
    
    // Test with different depth
    uint32_t key_depth5 = ph.generate_key(100, 200, 300, 5);
    std::cout << "Generated key for (100, 200, 300) with depth 5: " << key_depth5 << std::endl;
    
    // Test 64-bit version for higher depth
    uint64_t key64 = ph.generate_key_64(100, 200, 300, 15);
    std::cout << "Generated 64-bit key for (100, 200, 300) with depth 15: " << key64 << std::endl;
    
    // Print out the size of the NodeOrLeaf struct
    std::cout << "Size of NodeOrLeaf: " << sizeof(NodeOrLeaf) << std::endl;
    
    return 0;
} 