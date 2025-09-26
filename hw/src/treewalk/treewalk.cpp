#include "hls_math.h"
#include "tree_config.h"
#include "tree_types.h"

#ifndef __SYNTHESIS__
#include <thread>
#include <mutex>

  inline std::mutex& __stdout_mutex__() {
    static std::mutex m;
    return m;                 // one mutex per process for stdout
  }

#if defined(DEBUG) && (DEBUG)
#include <iostream>
#define SIM_PRINT(msg)                            \
    do {                                                    \
      std::lock_guard<std::mutex> lk(__stdout_mutex__());   \
      msg;                        \
    } while (0)
#else
#define SIM_PRINT(msg) do { } while (0)
#endif

#else
  #define SIM_PRINT(msg) do { } while (0)
#endif

// Simple work items
struct fk_work {
    particle_t particle;
    nodeleaf node;
    ap_uint<1> kill_signal;
};

bool eval_opening_criterion(double dx, double dy, double dz, double r2, double node_size, double thetasq) {
    #pragma HLS INLINE
    return r2 * thetasq < node_size;
}

// Pure force kernel - no memory access
void force_kernel(
    hls::stream<fk_work> &work_in,
    hls::stream<particle_t> &read_requests,
    pos_t thetasq, pos_t eps_sq)
{
    const pos_t node_sizes_sq[MAX_DEPTH] = {
        1.0 / (1 << 2),  1.0 / (1 << 4),  1.0 / (1 << 6),  1.0 / (1 << 8),
        1.0 / (1 << 10), 1.0 / (1 << 12), 1.0 / (1 << 14), 1.0 / (1 << 16),
        1.0 / (1 << 18), 1.0 / (1 << 20)};
    
    count_t num_particles_finished = 0;

    SIM_PRINT(std::cout << "FORCE_KERNEL STARTED" << std::endl);

    SIM_PRINT(std::cout << "node_sizes_sq[0]: " << node_sizes_sq[0] << std::endl);
    SIM_PRINT(std::cout << "node_sizes_sq[1]: " << node_sizes_sq[1] << std::endl);
    SIM_PRINT(std::cout << "node_sizes_sq[2]: " << node_sizes_sq[2] << std::endl);
    SIM_PRINT(std::cout << "node_sizes_sq[3]: " << node_sizes_sq[3] << std::endl);
    SIM_PRINT(std::cout << "node_sizes_sq[4]: " << node_sizes_sq[4] << std::endl);
    SIM_PRINT(std::cout << "node_sizes_sq[5]: " << node_sizes_sq[5] << std::endl);
    SIM_PRINT(std::cout << "node_sizes_sq[6]: " << node_sizes_sq[6] << std::endl);
    SIM_PRINT(std::cout << "node_sizes_sq[7]: " << node_sizes_sq[7] << std::endl);
    SIM_PRINT(std::cout << "node_sizes_sq[8]: " << node_sizes_sq[8] << std::endl);
    SIM_PRINT(std::cout << "node_sizes_sq[9]: " << node_sizes_sq[9] << std::endl);

    while(1) {
        #pragma HLS PIPELINE II=1
        
        // SIM_PRINT(std::cout << "+" << std::flush);
        
        fk_work work;
        if (!work_in.read_nb(work)) continue;

        if(work.kill_signal){
            particle_t p;
            p.valid = false;
            read_requests.write_nb(p);
            break;
        }
        
        particle_t p = work.particle;
        nodeleaf n = work.node;


        SIM_PRINT(std::cout << "p.acc= " << p.acc[0] << ", " << p.acc[1] << ", " << p.acc[2] << std::endl
        << "n.level: " << n.level << "n.start_idx: " << n.start_idx << std::endl);

        // Distance calculation
        double dx = double(n.pos[0]) - double(p.pos[0]);
        double dy = double(n.pos[1]) - double(p.pos[1]);
        double dz = double(n.pos[2]) - double(p.pos[2]);
        double r2 = dx * dx + dy * dy + dz * dz;

        SIM_PRINT(std::cout << "r2: " << r2 << " thetasq: " << thetasq << " node_sizes_sq[n.level-1]: " << node_sizes_sq[n.level-1] << std::endl);
        bool do_i_open = eval_opening_criterion(dx, dy, dz, r2, node_sizes_sq[n.level-1], thetasq);
        if (do_i_open && !n.is_leaf) {
            // Open node - go to first child
            p.next_tree_idx += 1;
        } else {
            // Accept node - accumulate force
            double r2_double = (double)r2;
            double rinv = hls::rsqrt(r2_double);
            double rinv2soft = 1.0 / (r2_double + (double)eps_sq);
            double rinv3 = rinv2soft * rinv;
            double force = (double)(n.mass) * rinv3;
            
            p.acc[0] += force * (double)dx;
            p.acc[1] += force * (double)dy;
            p.acc[2] += force * (double)dz;
            
            // Move to sibling
            if(n.next_sibling != -1u) {
                p.next_tree_idx += n.next_sibling;
            }
            else {
                p.next_tree_idx = -1u;
            }
        }
        read_requests.write_nb(p);
    }
    SIM_PRINT(std::cout << "FORCE_KERNEL FINISHED" << std::endl);
}

// Pure read kernel - only memory access
void read_kernel(
    hls::stream<particle_t> &requests,
    hls::stream<fk_work> &work_out,
    nodeleaf *tree)
{

    SIM_PRINT(std::cout << "READ_KERNEL STARTED" << std::endl);
    while(1) {
        #pragma HLS PIPELINE II=1
        
        particle_t p;
        if (!requests.read_nb(p)) continue;

        if(!p.valid) break;

        SIM_PRINT(std::cout << "read kernel got a request for particle " << p.idx << " and node" << p.next_tree_idx << std::endl);

        nodeleaf node;
        if(p.next_tree_idx != -1u) {
            node = tree[p.next_tree_idx];
        }

        // if p is -1u, it doesnt matter what node we send
        fk_work work = {p, node, false};

        work_out.write_nb(work);
    }
    SIM_PRINT(std::cout << "READ_KERNEL FINISHED" << std::endl);
}

// Work distributor - manages new particles and feedback loop
void work_distributor(
    hls::stream<particle_t> &new_particles,
    hls::stream<fk_work> &feedback,
    hls::stream<fk_work> &fk_queue, 
    hls::stream<particle_t> &output_particles, 
    nodeleaf first_node, const count_t NUM_PARTICLES)
{
    SIM_PRINT(std::cout << "WORK_DISTRIBUTOR STARTED" << std::endl);

    long num_particles_finished = 0;

    while(num_particles_finished < NUM_PARTICLES) {
        #pragma HLS PIPELINE II=1
        
        // SIM_PRINT(std::cout << "-" << std::flush);

        // Priority to feedback (continuing particles)
        fk_work work;
        if (feedback.read_nb(work)) {
            SIM_PRINT(std::cout << "read from feedback...");
            if(work.particle.next_tree_idx != -1u) {
                SIM_PRINT(std::cout << " writing to fk_queue..." << std::endl);
                fk_queue.write_nb(work);
            }
            else {
                SIM_PRINT(std::cout << "particle done, outputting" << std::endl);
                output_particles.write_nb(work.particle);
                num_particles_finished++;
                break;
            }
        }
        // Then new particles
        else {
            // SIM_PRINT(std::cout << "n" << std::flush);
            particle_t p;
            if (new_particles.read_nb(p)) {
                // Initialize particle
                p.next_tree_idx = 0;
                p.acc[0] = 0;
                p.acc[1] = 0;
                p.acc[2] = 0;
                p.valid = true;
                
                fk_work init_work = {p, first_node, false};
                SIM_PRINT(std::cout << "new particle writing to fk_queue..." << std::endl);
                fk_queue.write_nb(init_work);
            }
        }
    }

    SIM_PRINT(std::cout << "sending kill signal to force kernel" << std::endl);
    // now we need to send a kill signal to the force kernel
    fk_work kill_work = {particle_t(), nodeleaf(), true};
    fk_queue.write_nb(kill_work);

    SIM_PRINT(std::cout << "WORK_DISTRIBUTOR FINISHED" << std::endl);
}

void output_kernel(
    hls::stream<particle_t> &completed,
    hls::stream<particle_t> &particle_out,
    const count_t NUM_PARTICLES)
{
    SIM_PRINT(std::cout << "OUTPUT_KERNEL STARTED" << std::endl);
    // Write outputs
    long num_particles_written = 0;
    while(num_particles_written < NUM_PARTICLES) {
        #pragma HLS PIPELINE II=1
        particle_t p;
        if (completed.read_nb(p)) {
            particle_out.write_nb(p);
            num_particles_written++;
        }
    }
    SIM_PRINT(std::cout << "OUTPUT_KERNEL FINISHED" << std::endl);
}

// Main kernel
void treewalk_simple(
    hls::stream<particle_t> &particle_in,
    hls::stream<particle_t> &particle_out,
    nodeleaf *tree,
    pos_t theta, pos_t epsilon,
    const count_t NUM_PARTICLES)
{
#pragma HLS INTERFACE axis register port = particle_out depth=32
#pragma HLS INTERFACE axis register port = particle_in depth=32
#pragma HLS INTERFACE m_axi port = tree offset = slave bundle = gmem0 depth = 1024 max_widen_bitwidth = 512 max_write_burst_length = 64 num_write_outstanding = 16 latency = 64
#pragma HLS INTERFACE s_axilite port=NUM_PARTICLES bundle=control
#pragma HLS INTERFACE s_axilite port=tree bundle=control
#pragma HLS INTERFACE s_axilite port=theta bundle=control
#pragma HLS INTERFACE s_axilite port=epsilon bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control

    #pragma HLS DATAFLOW
    
    // Streams
    hls::stream<fk_work> fk_queue("fk_queue");
    #pragma HLS STREAM variable=fk_queue depth=32
    
    hls::stream<particle_t> rk_queue("rk_queue");
    #pragma HLS STREAM variable=rk_queue depth=32
    
    hls::stream<fk_work> feedback("feedback");
    #pragma HLS STREAM variable=feedback depth=32
    
    hls::stream<particle_t> completed("completed");
    #pragma HLS STREAM variable=completed depth=8
    
    pos_t thetasq = theta * theta;
    pos_t eps_sq = epsilon * epsilon;

    #ifndef __SYNTHESIS__
    std::thread work_distributor_thread(std::ref(work_distributor), std::ref(particle_in), std::ref(feedback), std::ref(fk_queue), std::ref(completed), tree[0], NUM_PARTICLES);
    std::thread force_kernel_thread(std::ref(force_kernel), std::ref(fk_queue), std::ref(rk_queue), thetasq, eps_sq);
    std::thread read_kernel_thread(std::ref(read_kernel), std::ref(rk_queue), std::ref(feedback), tree);
    std::thread output_kernel_thread(std::ref(output_kernel), std::ref(completed), std::ref(particle_out), NUM_PARTICLES);

    work_distributor_thread.join();
    force_kernel_thread.join();
    read_kernel_thread.join();
    output_kernel_thread.join();
    #else
    // Instantiate modules
    work_distributor(particle_in, feedback, fk_queue, tree[0], NUM_PARTICLES);
    force_kernel(fk_queue, rk_queue, thetasq, eps_sq);
    read_kernel(rk_queue, feedback, tree);
    output_kernel(completed, particle_out, NUM_PARTICLES);
    #endif
}