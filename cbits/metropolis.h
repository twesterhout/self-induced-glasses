#pragma once

#include <HsFFI.h>
#include <stddef.h>
#include <stdint.h>

float energy_change_upon_flip(HsInt number_bits, float const *couplings,
                              uint64_t const *state, HsInt i);

double run_one_sweep(HsInt number_bits, HsInt number_steps, float beta,
                     float const couplings[], HsInt const uniform_random_ints[],
                     float const uniform_random_floats[], uint64_t state[],
                     double *current_energy, float delta_energies[]);
