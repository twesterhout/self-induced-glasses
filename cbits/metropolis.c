#include "metropolis.h"
#include <math.h>

static inline float read_spin(uint64_t const *const restrict x, HsInt const i) {
  int const bit = (x[i / 64] >> (i % 64)) & 1;
  return 2 * bit - 1;
}

static inline void flip_spin(uint64_t *const restrict x, HsInt const i) {
  x[i / 64] ^= ((uint64_t)1) << (i % 64);
}

static inline float
packed_dot_product_scalar(HsInt const n, float const *const restrict a,
                          uint64_t const *const restrict b) {
  double acc = 0.0;
  for (ptrdiff_t i = 0; i < n; ++i) {
    acc += (double)a[i] * (double)read_spin(b, i);
  }
  return (float)acc;
}

float energy_change_upon_flip(HsInt const number_bits,
                              float const *const restrict couplings,
                              uint64_t const *const restrict state,
                              HsInt const i) {
  float const overlap = packed_dot_product_scalar(
      number_bits, couplings + i * number_bits, state);
  float const s = read_spin(state, i);
  return 2 * overlap * s;
}

double run_one_sweep(HsInt const number_bits, HsInt const number_steps,
                     float const beta, float const couplings[],
                     HsInt const uniform_random_ints[],
                     float const uniform_random_floats[], uint64_t state[],
                     double *current_energy, float delta_energies[]) {
  (void)delta_energies; // unused
  HsInt number_accepted = 0;
  double e = *current_energy;
  for (HsInt step_index = 0; step_index < number_steps; ++step_index) {
    HsInt const i = uniform_random_ints[step_index];
    float const de = energy_change_upon_flip(number_bits, couplings, state, i);
    if (de > 0) {
      float const u = uniform_random_floats[step_index];
      if (u < expf(-beta * de)) {
        flip_spin(state, i);
        e += de;
        ++number_accepted;
      }
    } else { // when de <= 0, always accept
      flip_spin(state, i);
      e += de;
      ++number_accepted;
    }
  }

  *current_energy = e;
  return (double)number_accepted / (double)number_steps;
}
