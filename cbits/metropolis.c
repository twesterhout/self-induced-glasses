#include "metropolis.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern uint64_t xoshiro256plusplus_next(uint64_t s[4]);

static inline uint64_t
generate_one_uint64(uint64_t const range,
                    uint64_t xoshiro256plusplus_state[4]) {
  __uint128_t random64bit = xoshiro256plusplus_next(xoshiro256plusplus_state);
  __uint128_t multiresult = random64bit * range;
  uint64_t leftover = (uint64_t)multiresult;
  if (leftover < range) {
    uint64_t threshold = (-range) % range;
    while (leftover < threshold) {
      random64bit = xoshiro256plusplus_next(xoshiro256plusplus_state);
      multiresult = random64bit * range;
      leftover = (uint64_t)multiresult;
    }
  }
  return (uint64_t)(multiresult >> 64); // [0, range)
}

static inline double generate_one_double(uint64_t xoshiro256plusplus_state[4]) {
  uint64_t const random64bit =
      xoshiro256plusplus_next(xoshiro256plusplus_state);
  return (random64bit >> 11) * 0x1.0p-53;
}

void xoshiro256plusplus_many_HsInt(HsInt const count, HsInt const upper,
                                   HsInt buffer[],
                                   uint64_t xoshiro256plusplus_state[4]) {
  for (HsInt i = 0; i < count; ++i) {
    buffer[i] =
        (HsInt)generate_one_uint64((uint64_t)upper, xoshiro256plusplus_state);
  }
}

void xoshiro256plusplus_many_float(HsInt const count, float buffer[],
                                   uint64_t xoshiro256plusplus_state[4]) {
  float const float_unit = 0x1.0p-24f;
  HsInt i = 0;
  for (; i < count - 1; i += 2) {
    uint64_t const random64bit =
        xoshiro256plusplus_next(xoshiro256plusplus_state);
    int const a = (random64bit >> 16) & 0xFFFFFF;
    int const b = random64bit >> (24 + 16);
    buffer[i] = a * float_unit;
    buffer[i + 1] = b * float_unit;
  }
  if (i < count) {
    uint64_t const random64bit =
        xoshiro256plusplus_next(xoshiro256plusplus_state);
    int const a = (random64bit >> 16) & 0xFFFFFF;
    buffer[i] = a * float_unit;
  }
}

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
                              float const *const restrict field,
                              uint64_t const *const restrict state,
                              HsInt const i) {
  float const overlap = packed_dot_product_scalar(
      number_bits, couplings + i * number_bits, state);
  float const h = (field != NULL) ? field[i] : 0.0f;
  float const s = read_spin(state, i);
  return -2 * (2 * overlap + h) * s;
}

static void recompute_energy_changes_upon_flip(
    HsInt const number_bits, float const *const restrict couplings,
    uint64_t const *const restrict state, HsInt const i,
    float *const restrict energy_changes) {

  float const pre = -8 * read_spin(state, i);
  for (HsInt k = 0; k < i; ++k) {
    float const coupling = couplings[i * number_bits + k];
    float const sigma = read_spin(state, k);
    energy_changes[k] += pre * coupling * sigma;
  }
  energy_changes[i] *= -1;
  for (HsInt k = i + 1; k < number_bits; ++k) {
    float const coupling = couplings[i * number_bits + k];
    float const sigma = read_spin(state, k);
    energy_changes[k] += pre * coupling * sigma;
  }
}

float total_energy(HsInt const n, float const *const restrict couplings,
                   float const *const restrict field,
                   uint64_t const *const restrict state) {
  float e = 0.0f;
  for (HsInt k = 0; k < n; ++k) {
    const float t = packed_dot_product_scalar(n, couplings + k * n, state);
    e += read_spin(state, k) * t;
  }
  if (field != NULL) {
    e += packed_dot_product_scalar(n, field, state);
  }
  return e;
}

HsInt total_magnetization(HsInt const n, uint64_t const *const restrict state) {
  HsInt m = 0;
  for (HsInt k = 0; k < n; ++k) {
    m += (int)read_spin(state, k);
  }
  return m;
}

void add_spin_spin_correlation(HsInt const n,
                               uint64_t const *const restrict state,
                               float *const restrict out) {
  for (HsInt i = 0; i < n; ++i) {
    out[i * n + i] += 1.0f;

    float const si = read_spin(state, i);
    for (HsInt j = i + 1; j < n; ++j) {
      float const sj = read_spin(state, j);
      float const correlation = si * sj;
      out[i * n + j] += correlation;
      out[j * n + i] += correlation;
    }
  }
}

double run_one_sweep(HsInt const number_bits, HsInt const number_steps,
                     float const beta, float const couplings[],
                     float const field[], HsInt const uniform_random_ints[],
                     float const uniform_random_floats[], uint64_t state[],
                     double *current_energy, float delta_energies[]) {
  (void)field;
  HsInt number_accepted = 0;
  double e = *current_energy;

  for (HsInt step_index = 0; step_index < number_steps; ++step_index) {
    HsInt const i = uniform_random_ints[step_index];
    float const de = delta_energies[i];
    // float const de_ref =
    //     energy_change_upon_flip(number_bits, couplings, field, state, i);
    // if (de != de_ref) {
    //   fprintf(stderr, "de (%f) != de_ref (%f) for i=%li\n", de, de_ref, i);
    //   abort();
    // }
    float const u = -logf(uniform_random_floats[step_index]) / beta;
    if (de < u) {
      flip_spin(state, i);
      recompute_energy_changes_upon_flip(number_bits, couplings, state, i,
                                         delta_energies);
      e += de;
      ++number_accepted;
    }
  }

  *current_energy = e;
  return (double)number_accepted / (double)number_steps;
}
