#pragma once

#include <cstdint>

// Global constants for BHTree
const int MAX_DEPTH = 10;
const uint32_t MAX_PARTICLES = (1U << 31);
const uint32_t MAX_NODES = (2 * MAX_PARTICLES); 
const uint32_t NLEAF = 1;