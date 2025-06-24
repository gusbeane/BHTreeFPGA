#pragma once

#include <array>
#include <vector>
#include <cstdint>

// Type definitions for BHTree
typedef std::array<uint32_t, 3> Position3D;
typedef std::vector<Position3D> PositionVector;
typedef std::array<double, 3> RealPosition3D;

// Node/Leaf structure for BHTree
struct NodeOrLeaf {
    // node information. in principle level and key can be in the same data type, but we keep them separate for now
    uint16_t level;
    unsigned int key; // uses 3 * MAX_DEPTH bits

    // particle information
    unsigned int start_idx; // start index in the particle array
    unsigned int Nparticles; // number of particles in the node
    
    // center of mass. if the node is a leaf and has only 1 particle, these store the particle's data
    // otherwise, for leaves you need to go back to the particle array
    double com_x;
    double com_y;
    double com_z;
    double mass;

    // whether or not the node is a leaf.
    // if the node is a leaf, then we simply read off the particle data from start_idx to start_idx + Nparticles
    bool is_leaf;
};

typedef std::vector<NodeOrLeaf> Tree;
