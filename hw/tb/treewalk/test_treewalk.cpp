#include "tree_config.h"
#include "tree_io.h"
#include "tree_types.h"
#include "treewalk.h"
#include "direct_sum.h"

int main() {
  // test tree walk with loaded tree data
  auto result = load_particles_and_tree("quasirandom_test_data");
  std::vector<particle_t> particles = result.first;
  std::vector<nodeleaf> tree = result.second;

  double h_part = particles[0].h;

  #if 1
  // For the testable version:
  particle_t input[1];
  particle_t output[1];

  // Copy first particle
  input[0].pos[0] = 0.2;
  input[0].pos[1] = 0.5;
  input[0].pos[2] = 0.5;
  input[0].mass = 1.0;
  input[0].idx = 0;
  input[0].next_tree_idx = 0;
  input[0].valid = true;
  input[0].acc[0] = 0; // Clear acceleration
  input[0].acc[1] = 0;
  input[0].acc[2] = 0;
  input[0].h = h_part;
  input[0].Nint_leaf = 0;
  input[0].Nint_node = 0;
  hls::stream<particle_t> particle_in("particle_in");
  for (int i = 0; i < 1; ++i) {
    particle_in.write(input[0]);
  }

  hls::stream<particle_t> particle_out("particle_out");

  // Run the testable kernel for a single particle
  treewalk_simple(particle_in, particle_out, tree.data(), pos_t(0.5),
                  pos_t(h_part), 1);

  // Read results
  count_t num_particles_received = 0;
  while (num_particles_received < 1) {
    particle_t p;
    if (particle_out.read_nb(p)) {
      output[num_particles_received] = p;
      num_particles_received++;
    }
  }

  // Print particle's position
  std::cout << "Particle 0 position: (" << output[0].pos[0] << ", "
            << output[0].pos[1] << ", " << output[0].pos[2] << ")" << std::endl;

  // Print result
  std::cout << "Treewalk force on particle 0: (" << output[0].acc[0] << ", "
            << output[0].acc[1] << ", " << output[0].acc[2] << ")" << std::endl;

  std::cout << "Nint_leaf: " << output[0].Nint_leaf << std::endl;
  std::cout << "Nint_node: " << output[0].Nint_node << std::endl;

  #endif

  // now try feeding in the first 200 particles
  const int NUM_PARTICLES = 200;
  std::vector<particle_t> particles_treeforce(NUM_PARTICLES);
  hls::stream<particle_t> particle_in_all("particle_in_all");

  for (int i = 0; i < NUM_PARTICLES; i++) {
    particle_in_all.write(particles[i]);
  }
  
  hls::stream<particle_t> particle_out_all("particle_out_all");
  treewalk_simple(particle_in_all, particle_out_all, tree.data(), pos_t(0.5),
                  pos_t(h_part), NUM_PARTICLES);
  
  count_t num_particles_received_all = 0;
  while (num_particles_received_all < NUM_PARTICLES) {
    particle_t p;
    if (particle_out_all.read_nb(p)) {
      particles_treeforce[num_particles_received_all] = p;
      num_particles_received_all++;
    }
  }

  // print acceleration of first particle
  std::cout << "acc for particle 0: " << particles_treeforce[0].acc[0] << ", " << particles_treeforce[0].acc[1] << ", " << particles_treeforce[0].acc[2] << std::endl;

  // Now compute direct sum force on first 200 particles
  std::vector<particle_t> particles_directsum(NUM_PARTICLES);
  std::vector<double> adiff(NUM_PARTICLES);
  for (int i = 0; i < NUM_PARTICLES; i++) {
    particles_directsum[i] = direct_sum_force(particles, particles[i], h_part);
    
    double adiff_i = pow(particles_directsum[i].acc[0] - particles_treeforce[i].acc[0], 2);
    adiff_i += pow(particles_directsum[i].acc[1] - particles_treeforce[i].acc[1], 2);
    adiff_i += pow(particles_directsum[i].acc[2] - particles_treeforce[i].acc[2], 2);
    adiff_i = std::sqrt(adiff_i);
    adiff[i] = adiff_i;
    // std::cout << "adiff[" << i << "]: " << adiff_i << std::endl;
    std::cout << "acc_x_tree[" << i << "]: " << particles_treeforce[i].acc[0] << " acc_x_direct[" << i << "]: " << particles_directsum[i].acc[0] << std::endl;
  }

  // print max adiff
  std::cout << "max adiff: " << *std::max_element(adiff.begin(), adiff.end()) << std::endl;

  // Now compute direct sum force
  // const int NUM_PARTICLES = 8000;

  // Create random number generator with fixed seed for reproducibility
  // std::mt19937 gen(42);
  // std::uniform_real_distribution<float> pos_dis(0.0f, 1.0f);

  double g = 1.22074408460575947536;
  double a1 = 1.0/g;
  double a2 = 1.0/(g*g);
  double a3 = 1.0/(g*g*g);

  std::vector<particle_t> particles_for_direct_sum(NUM_PARTICLES);
  particles_for_direct_sum[0].pos[0] = pos_t(0.5);
  particles_for_direct_sum[0].pos[1] = pos_t(0.5);
  particles_for_direct_sum[0].pos[2] = pos_t(0.5);

  double pos_dbl[3] = {0.5, 0.5, 0.5};

  for (int i = 1; i < NUM_PARTICLES; i++) {

    pos_dbl[0] = std::fmod(pos_dbl[0]+a1, 1.0);
    pos_dbl[1] = std::fmod(pos_dbl[1]+a2, 1.0);
    pos_dbl[2] = std::fmod(pos_dbl[2]+a3, 1.0);

    particles_for_direct_sum[i].pos[0] = pos_t(pos_dbl[0]);
    particles_for_direct_sum[i].pos[1] = pos_t(pos_dbl[1]);
    particles_for_direct_sum[i].pos[2] = pos_t(pos_dbl[2]);
  }
  
  for (int i = 0; i < NUM_PARTICLES; i++) {
    particles_for_direct_sum[i].mass = 1.0f / NUM_PARTICLES;
    particles_for_direct_sum[i].idx = i;
    particles_for_direct_sum[i].h = h_part;
  }

  double pos_for_direct_sum[3] = {0.2, 0.5, 0.5};
  particle_t p0;
  p0.pos[0] = pos_t(pos_for_direct_sum[0]);
  p0.pos[1] = pos_t(pos_for_direct_sum[1]);
  p0.pos[2] = pos_t(pos_for_direct_sum[2]);
  p0.mass = 1.0;
  p0.h = h_part;
  p0.idx = 0;

  p0 = direct_sum_force(particles_for_direct_sum, p0, h_part, true);

  std::cout << "Direct sum force on particle 0: (" << p0.acc[0] << ", "
            << p0.acc[1] << ", " << p0.acc[2] << ")" << std::endl;

  double acc_true[3] = {1.3676781978844061, 0.0, 0.0};
  std::cout << "relative difference in x: " << std::abs(p0.acc[0] - acc_true[0]) / acc_true[0] << std::endl;
  std::cout << "relative difference in y: " << std::abs(p0.acc[1] - acc_true[1]) / acc_true[0] << std::endl;
  std::cout << "relative difference in z: " << std::abs(p0.acc[2] - acc_true[2]) / acc_true[0] << std::endl;

  return 0;
}