#pragma once

#include <stddef.h>

float total_energy(ptrdiff_t n, float const *couplings, float const *x);

float energy_change_upon_flip(ptrdiff_t n, float const *couplings,
                              float const *x, ptrdiff_t i);

void recompute_energy_changes(ptrdiff_t n, float const *couplings,
                              float *delta_energies, ptrdiff_t i,
                              float const *x);

int autocorrelation(ptrdiff_t const n, float const *const signal,
                    float *const out);
