#pragma once

#include <HsFFI.h>
#include <stddef.h>
#include <stdint.h>

float total_energy(ptrdiff_t n, float const *couplings, uint64_t const *x);

float energy_change_upon_flip(HsInt n, float const *couplings,
                              uint64_t const *x, HsInt i);

void recompute_energy_changes(ptrdiff_t n, float const *couplings,
                              float *delta_energies, ptrdiff_t i,
                              float const *x);

ptrdiff_t hamming_weight(ptrdiff_t n, uint64_t const *x);

int autocorrelation(ptrdiff_t const n, float const *const signal,
                    float *const out);
