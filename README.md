# BHTreeFPGA

> Barnes-Hut Gravity Tree on an FPGA

## Overview

**BHTreeFPGA** implements the [Barnes-Hut algorithm](https://www.nature.com/articles/324446a0) for gravitational N-body simulations on FPGA hardware.

The Barnes-Hut method organizes particles into a hierarchical tree structure (an octree in 3D, or quadtree in 2D). Each node in the tree represents a spatial region and stores properties such as its center of mass. The domain is recursively subdivided until each leaf node contains fewer than `NLEAF` particles (typically 1–8).

Traditional CPU implementations build the tree by inserting particles one at a time, leading to inefficient memory access (pointer chasing). GPU algorithms like [Bonsai](https://github.com/treecode/Bonsai) improve on this with breadth-first, level-wise tree construction. However, they can be inefficient in highly clustered distributions or with very deep timebin hierarchies.

This project introduces a hardware-friendly streaming approach:

- **Space-filling curve ordering:** Particles are sorted by a Peano-Hilbert (or similar) key on the host CPU before being streamed to the FPGA. This guarantees that particles belonging to the same node are contiguous in memory.
- **On-FPGA node management:** The FPGA maintains a working list of tree nodes (up to `MAX_DEPTH`). As each particle arrives, it updates the relevant nodes. Completed nodes are communicated back to the host. This processing can be efficiently pipelined.
- **Tree-walking:** Work is in progress to implement tree traversal (force calculation) on the FPGA.

**Todo:**
- The design can successfully implement and meet the current timing target of `100 MHz`. Next steps: generate a bitstream and write sample host-side code.
- Nodes can be generated with a pipeline initiation interval (II) of 1. However, it is not yet clear how to efficiently transfer these nodes off the board. They are produced by the particle processing pipeline at a variable output rate (sometimes 0, sometimes up to `MAX_DEPTH`, with an average of ~1.4).
- Implement tree walking.

## Repository Structure
```
/
├── proto/      # Archival Python version
├── sw/         # Unmaintained software (CPU) version
├── hw/         # FPGA hardware implementation
│   ├── config/     # HLS configuration files
│   ├── src/        # Source code
│   ├── include/    # Headers
│   ├── scripts/    # HLS build scripts
│   └── tb/         # Testbenches
├── LICENSE
└── README.md   # Project overview (this file)
```
## License

This project is licensed under the [MIT License](LICENSE).