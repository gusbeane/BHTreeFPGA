#ifndef TREEWALK_H
#define TREEWALK_H

#include "hls_stream.h"
#include "tree_config.h"
#include "tree_types.h"

void treewalk_simple(hls::stream<particle_t> &particle_in,
                     hls::stream<particle_t> &particle_out, nodeleaf *tree, pos_t theta, pos_t epsilon, const count_t NUM_PARTICLES);

// define softening parameters, following gadget-4 convention
#define SOFTFAC1 (32.0 / 3) /**< Coefficients for gravitational softening */
#define SOFTFAC2 32.0
#define SOFTFAC3 (-38.4)
#define SOFTFAC4 (-2.8)
#define SOFTFAC5 (16.0 / 3)
#define SOFTFAC6 6.4
#define SOFTFAC7 (-9.6)
#define SOFTFAC8 (64.0 / 3)
#define SOFTFAC9 (-48.0)
#define SOFTFAC10 38.4
#define SOFTFAC11 (-32.0 / 3)
#define SOFTFAC12 (-1.0 / 15)
#define SOFTFAC13 (-3.2)
#define SOFTFAC14 (1.0 / 15)
#define SOFTFAC15 (-16.0)
#define SOFTFAC16 9.6
#define SOFTFAC17 (-64.0 / 30)
#define SOFTFAC18 128.0
#define SOFTFAC19 (-115.2)
#define SOFTFAC20 (64.0 / 3)
#define SOFTFAC21 (-96.0)
#define SOFTFAC22 115.2
#define SOFTFAC23 (-128.0 / 3)
#define SOFTFAC24 (4.0 / 30)

#define SOFTFAC30 (32.0 / 3)
#define SOFTFAC31 (-576.0 / 5)
#define SOFTFAC32 (128.0)
#define SOFTFAC33 (-1152.0 / 5)
#define SOFTFAC34 (384.0)
#define SOFTFAC35 (2.0 * 384.0)

#define SOFTFAC40 (64.0 / 3)
#define SOFTFAC41 (2.0 / 15)
#define SOFTFAC42 (-96.0)
#define SOFTFAC43 (576.0 / 5)
#define SOFTFAC44 (-128.0 / 3)
#define SOFTFAC45 (-96.0)
#define SOFTFAC46 (-2.0 / 5)
#define SOFTFAC47 (1152.0 / 5)
#define SOFTFAC48 (-128.0)
#define SOFTFAC49 (8.0 / 5)
#define SOFTFAC50 (-256.0)
#define SOFTFAC51 (-8.0)

double softW2(double r, double rinv, double h) {
    if(r > h) {
        return - rinv * rinv;
    }
    else{
        double u = r / h;
        double u2 = u * u;
        double h_inv = 1.0 / h;
        double h2_inv = h_inv * h_inv;
        if(u < 0.5) {
            return - h2_inv * u * (static_cast<double>(SOFTFAC1) + u2 * (static_cast<double>(SOFTFAC2)*u + static_cast<double>(SOFTFAC3)));
        }
        else{
            double u3 = u2 * u;
            return - h2_inv * u * (static_cast<double>(SOFTFAC8) + static_cast<double>(SOFTFAC9) * u + static_cast<double>(SOFTFAC10) * u2 + static_cast<double>(SOFTFAC11) * u3 + static_cast<double>(SOFTFAC12) / u3);
        }
    }
}

#endif