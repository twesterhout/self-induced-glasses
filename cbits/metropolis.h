#pragma once

#include <HsFFI.h>
#include <stddef.h>
#include <stdint.h>

float energy_change_upon_flip(HsInt number_bits, float const couplings[],
                              float const field[], uint64_t const state[],
                              HsInt i);

float total_energy(HsInt n, float const *couplings, float const *field,
                   uint64_t const *state);
HsInt total_magnetization(HsInt n, uint64_t const *state);

void add_spin_spin_correlation(HsInt n, uint64_t const *state, float *out);
void add_magnetization(HsInt n, uint64_t const *state, float *out);

float replica_overlap(HsInt n, uint64_t const *state1, uint64_t const *state2);
float replica_overlap_squared(HsInt n, uint64_t const *state1,
                              uint64_t const *state2);

double run_one_sweep(HsInt number_bits, HsInt number_steps, float beta,
                     float const couplings[], uint64_t state[],
                     float delta_energies[],
                     uint64_t xoshiro256plusplus_state[4],
                     double *current_energy);

void xoshiro256plusplus_many_HsInt(HsInt count, HsInt upper, HsInt buffer[],
                                   uint64_t xoshiro256plusplus_state[4]);

void xoshiro256plusplus_many_float(HsInt count, float buffer[],
                                   uint64_t xoshiro256plusplus_state[4]);
