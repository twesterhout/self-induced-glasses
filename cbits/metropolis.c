#include "metropolis.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <x86/avx2.h>
#include <time.h>

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

static inline simde__m256 unpack_byte(uint8_t const bits) {
  simde__m256i b1 = simde_mm256_set1_epi32((int32_t)bits); // broadcast bits
  simde__m256i m2 = simde_mm256_setr_epi32(1, 2, 4, 8, 0x10, 0x20, 0x40, 0x80);
  simde__m256i d1 =
      simde_mm256_and_si256(b1, m2); // isolate one bit in each dword
  simde__m256i mask = simde_mm256_cmpgt_epi32(
      d1, simde_mm256_setzero_si256()); // compare with 0
  return simde_mm256_castsi256_ps(mask);
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

static inline simde__m256 invert_mask(simde__m256 const mask,
                                      simde__m256i const minus_one) {
  return simde_mm256_castsi256_ps(
      simde_mm256_xor_si256(simde_mm256_castps_si256(mask), minus_one));
}

static inline simde__m256 new_energy_change(simde__m256 const change,
                                            simde__m256 const pre,
                                            simde__m256 const coupling,
                                            simde__m256 const mask,
                                            simde__m256 const not_mask) {
  return simde_mm256_add_ps(
      change,
      simde_mm256_mul_ps(
          pre, simde_mm256_sub_ps(simde_mm256_and_ps(coupling, mask),
                                  simde_mm256_and_ps(coupling, not_mask))));
}

static void recompute_energy_changes_upon_flip_simd(
    HsInt const number_bits, float const *const restrict couplings,
    uint64_t const *const restrict state, HsInt const i,
    float *const restrict energy_changes) {

  float const *const row = couplings + i * number_bits;
  float const pre = -8 * read_spin(state, i);
  simde__m256 const pre_v = simde_mm256_set1_ps(pre);
  simde__m256i const minus_one_v = simde_mm256_set1_epi32(-1);
  float const delta_energy = energy_changes[i];

  int const num_blocks = number_bits / 64;
  for (int block_idx = 0; block_idx < num_blocks; ++block_idx) {
    uint64_t word = state[block_idx];
    int offset = 64 * block_idx;

    // clang-format off
    simde__m256 const mask_v0 = unpack_byte(word & 0xFF);
    simde__m256 const mask_v1 = unpack_byte((word >> 8) & 0xFF);
    simde__m256 const mask_v2 = unpack_byte((word >> 16) & 0xFF);
    simde__m256 const mask_v3 = unpack_byte((word >> 24) & 0xFF);
    simde__m256 const mask_v4 = unpack_byte((word >> 32) & 0xFF);
    simde__m256 const mask_v5 = unpack_byte((word >> 40) & 0xFF);
    simde__m256 const mask_v6 = unpack_byte((word >> 48) & 0xFF);
    simde__m256 const mask_v7 = unpack_byte((word >> 56) & 0xFF);

    simde__m256 const not_mask_v0 = invert_mask(mask_v0, minus_one_v);
    simde__m256 const not_mask_v1 = invert_mask(mask_v1, minus_one_v);
    simde__m256 const not_mask_v2 = invert_mask(mask_v2, minus_one_v);
    simde__m256 const not_mask_v3 = invert_mask(mask_v3, minus_one_v);
    simde__m256 const not_mask_v4 = invert_mask(mask_v4, minus_one_v);
    simde__m256 const not_mask_v5 = invert_mask(mask_v5, minus_one_v);
    simde__m256 const not_mask_v6 = invert_mask(mask_v6, minus_one_v);
    simde__m256 const not_mask_v7 = invert_mask(mask_v7, minus_one_v);

    simde__m256 const coupling_v0 = simde_mm256_loadu_ps(row + offset + 8 * 0);
    simde__m256 const coupling_v1 = simde_mm256_loadu_ps(row + offset + 8 * 1);
    simde__m256 const coupling_v2 = simde_mm256_loadu_ps(row + offset + 8 * 2);
    simde__m256 const coupling_v3 = simde_mm256_loadu_ps(row + offset + 8 * 3);
    simde__m256 const coupling_v4 = simde_mm256_loadu_ps(row + offset + 8 * 4);
    simde__m256 const coupling_v5 = simde_mm256_loadu_ps(row + offset + 8 * 5);
    simde__m256 const coupling_v6 = simde_mm256_loadu_ps(row + offset + 8 * 6);
    simde__m256 const coupling_v7 = simde_mm256_loadu_ps(row + offset + 8 * 7);

    simde__m256 const change_v0 = simde_mm256_loadu_ps(energy_changes + offset + 8 * 0);
    simde__m256 const change_v1 = simde_mm256_loadu_ps(energy_changes + offset + 8 * 1);
    simde__m256 const change_v2 = simde_mm256_loadu_ps(energy_changes + offset + 8 * 2);
    simde__m256 const change_v3 = simde_mm256_loadu_ps(energy_changes + offset + 8 * 3);
    simde__m256 const change_v4 = simde_mm256_loadu_ps(energy_changes + offset + 8 * 4);
    simde__m256 const change_v5 = simde_mm256_loadu_ps(energy_changes + offset + 8 * 5);
    simde__m256 const change_v6 = simde_mm256_loadu_ps(energy_changes + offset + 8 * 6);
    simde__m256 const change_v7 = simde_mm256_loadu_ps(energy_changes + offset + 8 * 7);

    simde__m256 const new_change_v0 = new_energy_change(change_v0, pre_v, coupling_v0, mask_v0, not_mask_v0);
    simde__m256 const new_change_v1 = new_energy_change(change_v1, pre_v, coupling_v1, mask_v1, not_mask_v1);
    simde__m256 const new_change_v2 = new_energy_change(change_v2, pre_v, coupling_v2, mask_v2, not_mask_v2);
    simde__m256 const new_change_v3 = new_energy_change(change_v3, pre_v, coupling_v3, mask_v3, not_mask_v3);
    simde__m256 const new_change_v4 = new_energy_change(change_v4, pre_v, coupling_v4, mask_v4, not_mask_v4);
    simde__m256 const new_change_v5 = new_energy_change(change_v5, pre_v, coupling_v5, mask_v5, not_mask_v5);
    simde__m256 const new_change_v6 = new_energy_change(change_v6, pre_v, coupling_v6, mask_v6, not_mask_v6);
    simde__m256 const new_change_v7 = new_energy_change(change_v7, pre_v, coupling_v7, mask_v7, not_mask_v7);

    simde_mm256_storeu_ps(energy_changes + offset + 8 * 0, new_change_v0);
    simde_mm256_storeu_ps(energy_changes + offset + 8 * 1, new_change_v1);
    simde_mm256_storeu_ps(energy_changes + offset + 8 * 2, new_change_v2);
    simde_mm256_storeu_ps(energy_changes + offset + 8 * 3, new_change_v3);
    simde_mm256_storeu_ps(energy_changes + offset + 8 * 4, new_change_v4);
    simde_mm256_storeu_ps(energy_changes + offset + 8 * 5, new_change_v5);
    simde_mm256_storeu_ps(energy_changes + offset + 8 * 6, new_change_v6);
    simde_mm256_storeu_ps(energy_changes + offset + 8 * 7, new_change_v7);
    // clang-format on
  }

  for (int k = 64 * num_blocks; k < number_bits; ++k) {
    float const coupling = couplings[i * number_bits + k];
    float const sigma = read_spin(state, k);
    energy_changes[k] += pre * coupling * sigma;
  }

  energy_changes[i] = -delta_energy;
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

float replica_overlap(HsInt const n, uint64_t const *const restrict state1,
                      uint64_t const *const restrict state2) {
  float out = 0.0f;
  for (HsInt i = 0; i < n; ++i) {
    float const s1_i = read_spin(state1, i);
    float const s2_i = read_spin(state2, i);
    out += s1_i * s2_i;
  }
  return out;
}

float replica_overlap_squared(HsInt const n,
                              uint64_t const *const restrict state1,
                              uint64_t const *const restrict state2) {
  float out = 0.0f;
  for (HsInt i = 0; i < n; ++i) {
    float const s1_i = read_spin(state1, i);
    float const s2_i = read_spin(state2, i);
    for (HsInt j = i + 1; j < n; ++j) {
      float const s1_j = read_spin(state1, j);
      float const s2_j = read_spin(state2, j);
      out += s1_i * s1_j * s2_i * s2_j;
    }
  }
  return out;
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

void add_magnetization(HsInt n, uint64_t const *state, float *out) {
  for (HsInt i = 0; i < n; ++i) {
    out[i] += read_spin(state, i);
  }
}

double run_one_sweep(HsInt const number_bits, HsInt const number_steps,
                     float const beta, float const couplings[],
                     uint64_t state[], float delta_energies[],
                     uint64_t xoshiro256plusplus_state[4],
                     double *current_energy) {
  // struct timespec tp1;
  // clock_gettime(CLOCK_MONOTONIC, &tp1);

  HsInt number_accepted = 0;
  double e = *current_energy;

  int const block_size = 64;
  HsInt uniform_random_ints[block_size];
  float uniform_random_floats[block_size];

  for (HsInt step_index = 0; step_index < number_steps; ++step_index) {
    int const rand_index = step_index % block_size;
    if (rand_index == 0) {
      xoshiro256plusplus_many_HsInt(block_size, number_bits,
                                    uniform_random_ints,
                                    xoshiro256plusplus_state);
      xoshiro256plusplus_many_float(block_size, uniform_random_floats,
                                    xoshiro256plusplus_state);
    }

    HsInt const i = uniform_random_ints[rand_index];
    float const de = delta_energies[i];
    // float const de_ref =
    //     energy_change_upon_flip(number_bits, couplings, field, state, i);
    // if (de != de_ref) {
    //   fprintf(stderr, "de (%f) != de_ref (%f) for i=%li\n", de, de_ref, i);
    //   abort();
    // }
    float const u = -logf(uniform_random_floats[rand_index]) / beta;
    if (de < u) {
      flip_spin(state, i);
      recompute_energy_changes_upon_flip_simd(number_bits, couplings, state, i,
                                              delta_energies);
      e += de;
      ++number_accepted;
    }
  }

  *current_energy = e;

  // struct timespec tp2;
  // clock_gettime(CLOCK_MONOTONIC, &tp2);

  // fprintf(stderr, "run_one_sweep(number_steps=%zi, beta=%f) took %f\n",
  //         number_steps, (double)beta,
  //         (tp2.tv_sec - tp1.tv_sec) + (double)(tp2.tv_nsec - tp1.tv_nsec) / 1e9);
  return (double)number_accepted / (double)number_steps;
}
