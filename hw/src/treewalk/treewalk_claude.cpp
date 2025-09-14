#include "treecon_types_hls.h"
#include "treecon_config_hls.h"
#include <cstring>
#include <cmath>

#ifndef __SYNTHESIS__
#include <iostream>
#endif

// Constants for the systolic array
const int ARRAY_SIZE = 64;  // Number of processing cells in the array
const double THETA = 0.5;   // Opening angle criterion
const double G = 1.0;       // Gravitational constant

// Unified packet type that can carry either particles or nodes
struct unified_packet {
    // Packet metadata
    ap_uint<2> packet_type;   // 0=invalid, 1=particle, 2=node_update
    count_t target_cell;      // Which cell should process this packet
    
    // Union-like data payload (both fields present, but only one used based on packet_type)
    particle_t particle;
    nodeleaf node;
    
    // Constructor for invalid packet
    unified_packet() {
        packet_type = 0;  // Invalid
        target_cell = 0;
    }
    
    // Static factory methods for type safety
    static unified_packet make_particle(particle_t p, count_t target) {
        unified_packet pkt;
        pkt.packet_type = 1;  // Particle
        pkt.target_cell = target;
        pkt.particle = p;
        return pkt;
    }
    
    static unified_packet make_node(nodeleaf n, count_t target) {
        unified_packet pkt;
        pkt.packet_type = 2;  // Node update
        pkt.target_cell = target;
        pkt.node = n;
        return pkt;
    }
    
    // Helper methods
    bool is_valid() const { return packet_type != 0; }
    bool is_particle() const { return packet_type == 1; }
    bool is_node() const { return packet_type == 2; }
};

// Single processing cell in the systolic array
void processing_cell(
    count_t cell_index,
    hls::stream<unified_packet>& data_in,
    hls::stream<unified_packet>& data_out
) {
#pragma HLS INLINE off
#pragma HLS PIPELINE II=1

    // Cell state - stored in registers/BRAM
    static nodeleaf cell_node;
    static ap_uint<1> cell_initialized = 0;
    
    // Process incoming packets
    if (!data_in.empty()) {
        unified_packet pkt = data_in.read();
        
        if (pkt.is_node() && pkt.target_cell == cell_index) {
            // This node update is for us - load it into our cell state
            cell_node = pkt.node;
            cell_initialized = 1;
            
#ifndef __SYNTHESIS__
            std::cout << "Cell " << cell_index << " loaded node with key " 
                      << pkt.node.key << std::endl;
#endif
            // Node updates are consumed when they reach their target, not forwarded
            
        } else if (pkt.is_particle() && pkt.target_cell == cell_index && cell_initialized) {
            // This particle is for us to process
            
            // Compute distance from particle to node's center of mass
            double dx = cell_node.pos[0] - pkt.particle.pos[0];
            double dy = cell_node.pos[1] - pkt.particle.pos[1]; 
            double dz = cell_node.pos[2] - pkt.particle.pos[2];
            double rsq = dx*dx + dy*dy + dz*dz;
            double r = sqrt(rsq);
            
            // Check multipole acceptance criterion
            double node_size = 1.0 / (1 << cell_node.level);
            ap_uint<1> open_node = (r <= node_size/THETA) && !cell_node.is_leaf;
            
            if (open_node) {
                // Open the node - send particle to next cell
                unified_packet out_pkt = unified_packet::make_particle(pkt.particle, cell_index + 1);
                data_out.write(out_pkt);
                
#ifndef __SYNTHESIS__
                std::cout << "Cell " << cell_index << " opening node, sending to cell " 
                          << (cell_index + 1) << std::endl;
#endif
            } else {
                // Use this node for force calculation
                double rcubed = rsq * r;
                if (r > 0) {
                    double mass_double = cell_node.mass; // Convert ap_ufixed to double
                    pkt.particle.acc[0] += G * mass_double * dx / rcubed;
                    pkt.particle.acc[1] += G * mass_double * dy / rcubed;
                    pkt.particle.acc[2] += G * mass_double * dz / rcubed;
                }
                
                // Skip to next sibling (or end if this was the last)
                count_t next_target;
                if (cell_node.next_sibling != -1u) {
                    next_target = cell_index + cell_node.next_sibling;
                } else {
                    next_target = ARRAY_SIZE; // End processing
                }
                
                if (next_target < ARRAY_SIZE) {
                    unified_packet out_pkt = unified_packet::make_particle(pkt.particle, next_target);
                    data_out.write(out_pkt);
                }
                
#ifndef __SYNTHESIS__
                std::cout << "Cell " << cell_index << " computed force, next target: " 
                          << next_target << std::endl;
#endif
            }
            
        } else if (pkt.is_valid() && pkt.target_cell > cell_index) {
            // Pass packet along unchanged - it's meant for a downstream cell
            data_out.write(pkt);
            
        } else if (pkt.is_valid() && pkt.target_cell < cell_index) {
            // This shouldn't happen in normal operation
#ifndef __SYNTHESIS__
            std::cout << "Warning: Cell " << cell_index 
                      << " received packet for cell " << pkt.target_cell << std::endl;
#endif
        }
        // If packet is invalid (packet_type == 0), we just drop it
    }
}

// Systolic array of processing cells  
void systolic_tree_array(
    hls::stream<unified_packet>& data_in,
    hls::stream<unified_packet>& data_out
) {
#pragma HLS DATAFLOW

    // Create unified streams connecting the cells (much simpler!)
    hls::stream<unified_packet> data_streams[ARRAY_SIZE - 1];
    
#pragma HLS STREAM variable=data_streams depth=8

    // First cell takes from input stream
    processing_cell(0, data_in, data_streams[0]);
    
    // Middle cells connect through internal streams
    for (int i = 1; i < ARRAY_SIZE - 1; i++) {
#pragma HLS UNROLL
        processing_cell(i, data_streams[i-1], data_streams[i]);
    }
    
    // Last cell outputs to output stream
    processing_cell(ARRAY_SIZE - 1, data_streams[ARRAY_SIZE - 2], data_out);
}

// Stream merger that combines node loading and particle processing
void stream_controller(
    const ap_uint<512> *tree_memory,
    particle_t *particles,
    count_t num_nodes,
    count_t num_particles,
    hls::stream<unified_packet>& data_stream_out,
    hls::stream<unified_packet>& data_stream_in
) {
#pragma HLS DATAFLOW

    // Phase 1: Load all tree nodes first
    NODE_LOAD:
    for (count_t i = 0; i < num_nodes && i < ARRAY_SIZE; i++) {
#pragma HLS PIPELINE II=1
        
        // Read node from memory
        ap_uint<512> node_data = tree_memory[i];
        
        // Unpack the node data into nodeleaf struct
        // (This would need proper bit packing/unpacking logic based on your memory layout)
        nodeleaf node;
        node.key = node_data(29, 0);
        node.level = node_data(61, 30);
        // TODO: Complete unpacking for pos[3], mass, start_idx, num_particles, next_sibling, is_leaf, is_last
        // Example (you'll need to adjust bit ranges based on your actual packing):
        // node.pos[0] = node_data(93, 62);   // 32 bits for pos[0]
        // node.pos[1] = node_data(125, 94);  // 32 bits for pos[1] 
        // node.pos[2] = node_data(157, 126); // 32 bits for pos[2]
        // ... etc
        
        // Create node update packet
        unified_packet node_pkt = unified_packet::make_node(node, i);
        data_stream_out.write(node_pkt);
    }
    
    // Phase 2: Process all particles
    PARTICLE_PROCESS:
    for (count_t i = 0; i < num_particles; i++) {
#pragma HLS PIPELINE II=1
        
        // Send particle into the array
        unified_packet particle_pkt = unified_packet::make_particle(particles[i], 0);
        data_stream_out.write(particle_pkt);
        
        // Receive processed particle back
        unified_packet result_pkt = data_stream_in.read();
        if (result_pkt.is_particle()) {
            particles[i] = result_pkt.particle;
        }
    }
}

// Top-level kernel function
void walk_bhtree_kernel(
    const ap_uint<512> *tree,
    particle_t *particles,
    count_t num_particles,
    count_t num_nodes
) {
#pragma HLS INTERFACE m_axi port = particles bundle = gmem0 depth = 1024
#pragma HLS INTERFACE m_axi port = tree offset = slave bundle = gmem1 depth = 1024 max_widen_bitwidth = 512 max_write_burst_length = 64 num_write_outstanding = 16 latency = 64
#pragma HLS INTERFACE s_axilite port=particles bundle=control
#pragma HLS INTERFACE s_axilite port=num_particles bundle=control
#pragma HLS INTERFACE s_axilite port=num_nodes bundle=control
#pragma HLS INTERFACE s_axilite port=tree bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control

#pragma HLS DATAFLOW

    // Create unified communication streams - much cleaner!
    hls::stream<unified_packet> data_to_array("data_to_array");
    hls::stream<unified_packet> data_from_array("data_from_array");
    
#pragma HLS STREAM variable=data_to_array depth=64
#pragma HLS STREAM variable=data_from_array depth=64

    // Stream controller handles both node loading and particle processing
    stream_controller(tree, particles, num_nodes, num_particles, 
                     data_to_array, data_from_array);
    
    // Run the systolic array with unified streams
    systolic_tree_array(data_to_array, data_from_array);
}
