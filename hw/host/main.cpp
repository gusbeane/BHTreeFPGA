#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <string>
#include <fstream>
#include <random>

// Add OpenCL version macros before including OpenCL headers to target OpenCL 1.2 provided by Xilinx
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_ENABLE_EXCEPTIONS

// OpenCL headers
#include <CL/cl2.hpp>

// HLS headers
#include "../include/tree_config.h"
#include "../include/tree_types.h"

// Test particle structure
struct TestParticle {
    double x, y, z, mass;
};

// Helper function to create a test particle
particle_t create_test_particle(double x, double y, double z, double mass, 
                                uint32_t ph_key, int index) {
    particle_t p;
    p.pos[0] = pos_t(x);
    p.pos[1] = pos_t(y); 
    p.pos[2] = pos_t(z);
    p.mass = mass_t(mass);
    p.idx = count_t(index);
    p.key = phkey_t(ph_key & 0x3FFFFFFF); // Mask to 30 bits for phkey_t
    return p;
}

// Helper function to convert ap_uint<512> back to nodeleaf for analysis
nodeleaf convert_output_node(const ap_uint<512>& node_bits) {
    return *reinterpret_cast<const nodeleaf*>(&node_bits);
}

// Helper function to print a particle for debugging
void print_particle(const particle_t& p, int index) {
    std::cout << "Particle " << index << ": "
              << "pos=(" << double(p.pos[0]) << "," << double(p.pos[1]) << "," << double(p.pos[2]) << ") "
              << "mass=" << double(p.mass) << " "
              << "key=0x" << std::hex << int(p.key) << std::dec << " "
              << "idx=" << int(p.idx) << std::endl;
}

// Helper function to print a node for debugging  
void print_node(const nodeleaf& node, int index) {
    std::cout << "Node " << index << ": "
              << "level=" << int(node.level) << " "
              << "key=0x" << std::hex << int(node.key) << std::dec << " "
              << "npart=" << int(node.num_particles) << " "
              << "start_idx=" << int(node.start_idx) << " "
              << "pos=(" << double(node.pos[0]) << "," << double(node.pos[1]) << "," << double(node.pos[2]) << ") "
              << "mass=" << double(node.mass) << " "
              << "leaf=" << (node.is_leaf ? "true" : "false") << " "
              << "last=" << (node.is_last ? "true" : "false") << std::endl;
}

// Simple function to compute Peano-Hilbert key (simplified 3D version)
uint32_t simple_peano_hilbert_key(double x, double y, double z, int depth = 10) {
    // Convert to integer coordinates (0 to 2^depth - 1)
    uint32_t max_coord = (1U << depth) - 1;
    uint32_t ix = (uint32_t)(x * max_coord);
    uint32_t iy = (uint32_t)(y * max_coord);
    uint32_t iz = (uint32_t)(z * max_coord);
    
    // Very simplified PH key - just interleave bits for now
    uint32_t key = 0;
    for (int i = 0; i < depth; i++) {
        uint32_t bit_x = (ix >> i) & 1;
        uint32_t bit_y = (iy >> i) & 1;
        uint32_t bit_z = (iz >> i) & 1;
        key |= (bit_x << (3*i)) | (bit_y << (3*i + 1)) | (bit_z << (3*i + 2));
    }
    return key;
}

std::vector<char> load_file(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }
    
    file.seekg(0, std::ios::end);
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);
    
    std::vector<char> buffer(size);
    file.read(buffer.data(), size);
    return buffer;
}

int main(int argc, char* argv[]) {
    std::cout << "=== FPGA BHTree Host Application ===" << std::endl;
    
    std::string xclbin_file = "../build/bhtree/create_bhtree_kernel.xclbin";
    if (argc > 1) {
        xclbin_file = argv[1];
    }
    
    try {
        // Initialize OpenCL
        std::vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
        
        if (platforms.empty()) {
            throw std::runtime_error("No OpenCL platforms found");
        }
        
        cl::Platform platform;
        bool found_xilinx = false;
        for (const auto &p : platforms) {
          std::string platform_name = p.getInfo<CL_PLATFORM_NAME>();
          std::cout << "Found platform: " << platform_name << std::endl;
          if (platform_name.find("Xilinx") != std::string::npos) {
            platform = p;
            found_xilinx = true;
            break;
          }
        }

        if (!found_xilinx) {
          throw std::runtime_error("No Xilinx platform found");
        }
        std::cout << "Using platform: " << platform.getInfo<CL_PLATFORM_NAME>() << std::endl;
        
        std::vector<cl::Device> devices;
        platform.getDevices(CL_DEVICE_TYPE_ACCELERATOR, &devices);
        
        if (devices.empty()) {
            throw std::runtime_error("No accelerator devices found");
        }
        
        cl::Device device = devices[0];
        std::cout << "Using device: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
        
        cl::Context context(device);
        cl::CommandQueue queue(context, device, CL_QUEUE_PROFILING_ENABLE);
        
        // Load xclbin
        std::cout << "Loading xclbin file: " << xclbin_file << std::endl;
        std::vector<char> xclbin_data = load_file(xclbin_file);

        // Convert char vector to unsigned char vector for OpenCL
        std::vector<unsigned char> xclbin_binary(xclbin_data.begin(), xclbin_data.end());
        cl::Program::Binaries binaries = {xclbin_binary};
        std::vector<cl::Device> prog_devices = {device};
        cl::Program program(context, prog_devices, binaries);
        
        // Create kernel
        cl::Kernel kernel(program, "create_bhtree_kernel");
        
        // Prepare test data - 1 million random particles in the unit box
        const int NUM_PARTICLES = 1000000;
        const int TREE_SIZE = 3 * NUM_PARTICLES;

        std::vector<TestParticle> test_particles;
        test_particles.reserve(NUM_PARTICLES);

        std::mt19937 rng(42); // Fixed seed for reproducibility
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double mass = 1.0 / NUM_PARTICLES;

        for (int i = 0; i < NUM_PARTICLES; ++i) {
            double x = dist(rng);
            double y = dist(rng);
            double z = dist(rng);
            test_particles.push_back({x, y, z, mass});
        }
        
        // Compute PH keys and create sorted list
        std::vector<std::pair<uint32_t, int>> key_index_pairs;
        for (int i = 0; i < NUM_PARTICLES; i++) {
            uint32_t ph_key = simple_peano_hilbert_key(test_particles[i].x, 
                                                       test_particles[i].y, 
                                                       test_particles[i].z);
            key_index_pairs.push_back({ph_key, i});
        }
        
        // Sort by PH key
        std::sort(key_index_pairs.begin(), key_index_pairs.end());
        
        std::cout << "\nParticles sorted by PH keys:" << std::endl;
        std::vector<particle_t> particles(NUM_PARTICLES);
        for (int i = 0; i < NUM_PARTICLES; i++) {
            int orig_idx = key_index_pairs[i].second;
            uint32_t key = key_index_pairs[i].first;
            particles[i] = create_test_particle(test_particles[orig_idx].x,
                                               test_particles[orig_idx].y,
                                               test_particles[orig_idx].z,
                                               test_particles[orig_idx].mass,
                                               key, i);
            // print_particle(particles[i], i);
        }
        
        // Create OpenCL buffers
        cl::Buffer particles_buffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                   sizeof(particle_t) * NUM_PARTICLES, particles.data());
        
        cl::Buffer tree_buffer(context, CL_MEM_WRITE_ONLY,
                              sizeof(ap_uint<512>) * TREE_SIZE);
        
        // Set kernel arguments
        kernel.setArg(0, particles_buffer);
        kernel.setArg(1, tree_buffer);
        kernel.setArg(2, NUM_PARTICLES);
        
        // Execute kernel
        std::cout << "\nExecuting kernel on FPGA..." << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        
        cl::Event event;
        queue.enqueueTask(kernel, NULL, &event);
        event.wait();
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        
        std::cout << "Kernel execution time: " << duration.count() << " microseconds" << std::endl;
        std::cout << "Particles per microsecond: " << static_cast<float>(NUM_PARTICLES) / duration.count() << std::endl;
        
        // Read back results
        std::vector<ap_uint<512>> tree_output(TREE_SIZE);
        queue.enqueueReadBuffer(tree_buffer, CL_TRUE, 0, 
                               sizeof(ap_uint<512>) * TREE_SIZE, tree_output.data());
        
        // Analyze output
        int num_nodes = 0;
        std::cout << "\nOutput nodes:" << std::endl;
        for (int i = 0; i < TREE_SIZE; i++) {
            nodeleaf node = convert_output_node(tree_output[i]);
            
            // Stop if we encounter a node with zero mass and zero particles
            if (node.mass == 0.0 && node.num_particles == 0) {
                break;
            }
            
            num_nodes++;
            // print_node(node, i);
            if (node.is_last) break;
        }
        
        std::cout << "\nFPGA execution produced " << num_nodes << " nodes" << std::endl;
        
        // Basic sanity checks
        bool test_passed = true;
        
        if (num_nodes == 0) {
            std::cout << "ERROR: No output nodes generated!" << std::endl;
            test_passed = false;
        }
        
        // Check that masses are reasonable
        double total_mass = 0.0;
        for (int i = 0; i < num_nodes; i++) {
            nodeleaf node = convert_output_node(tree_output[i]);
            if(node.is_leaf) {
                total_mass += double(node.mass);
            }
        }
        
        std::cout << "Total mass in tree: " << total_mass << " (expected: 1.0)" << std::endl;
        if (std::abs(total_mass - 1.0) > 0.01) {
            std::cout << "ERROR: Total mass mismatch!" << std::endl;
            test_passed = false;
        }
        
        if (test_passed) {
            std::cout << "\n✅ FPGA execution PASSED!" << std::endl;
            return 0;
        } else {
            std::cout << "\n❌ FPGA execution FAILED!" << std::endl;
            return 1;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
