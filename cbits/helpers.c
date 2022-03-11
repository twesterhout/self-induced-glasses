#include "helpers.h"
#include <cblas.h>
#include <complex.h>
#include <fftw3.h>
#include <stdlib.h>

// float effective_interaction_loop(sig_model const *model, int const px,
//                                  int const py, int const radius) {
//
// }

static inline float read_spin(uint64_t const *const restrict x,
                              ptrdiff_t const i) {
  int const bit = (x[i / 64] >> (i % 64)) & 1;
  return 2 * bit - 1;
}

static inline float packed_dot_product(ptrdiff_t const n,
                                       float const *const restrict a,
                                       uint64_t const *const restrict b) {
  float acc = 0.0f;
  for (ptrdiff_t i = 0; i < n; ++i) {
    acc += a[i] * read_spin(b, i);
  }
  return acc;
}

float energy_change_upon_flip(ptrdiff_t const n,
                              float const *const restrict couplings,
                              uint64_t const *const restrict x,
                              ptrdiff_t const i) {
  float const overlap = packed_dot_product(n, couplings + i * n, x);
  float const s = read_spin(x, i);
  return 2 * overlap * s;
}

float total_energy(ptrdiff_t const n, float const *const restrict couplings,
                   uint64_t const *const restrict x) {
  float e = 0.0f;
  for (ptrdiff_t k = 0; k < n; ++k) {
    const float t = packed_dot_product(n, couplings + k * n, x);
    e += read_spin(x, k) * t;
  }
  return -e;
}

ptrdiff_t hamming_weight(ptrdiff_t const n, uint64_t const *const restrict x) {
  ptrdiff_t m = 0;
  ptrdiff_t const number_words = n / 64;
  ptrdiff_t const rest = n % 64;
  for (ptrdiff_t k = 0; k < number_words; ++k) {
    m += __builtin_popcountl(x[k]);
  }
  if (rest != 0) {
    uint64_t const mask = (((uint64_t)1) << rest) - 1;
    uint64_t const t = x[number_words];
    m += __builtin_popcountl(t & mask);
  }
  return m;
}

void recompute_energy_changes(ptrdiff_t const n,
                              float const *const restrict couplings,
                              float *const restrict delta_energies,
                              ptrdiff_t const i,
                              float const *const restrict x) {
  float const pre = -8 * x[i];
  for (ptrdiff_t k = 0; k < n; ++k) {
    const float coupling = couplings[i * n + k];
    const float sigma = x[k];
    delta_energies[k] += pre * coupling * sigma;
  }
}

static void fft(ptrdiff_t const n, _Complex float const *in,
                _Complex float *out) {
  // NOTE: the cast is relatively safe because we're using FFTW_ESTIMATE which
  // should not touch the arrays
  fftwf_plan p = fftwf_plan_dft_1d(n, (_Complex float *)in, out, FFTW_FORWARD,
                                   FFTW_ESTIMATE);
  fftwf_execute(p);
  fftwf_destroy_plan(p);
}

static void ifft(ptrdiff_t const n, _Complex float const *in,
                 _Complex float *out) {
  // NOTE: the cast is relatively safe because we're using FFTW_ESTIMATE which
  // should not touch the arrays
  fftwf_plan p = fftwf_plan_dft_1d(n, (_Complex float *)in, out, FFTW_BACKWARD,
                                   FFTW_ESTIMATE);
  fftwf_execute(p);
  fftwf_destroy_plan(p);
}

int autocorrelation(ptrdiff_t const n, float const *const restrict signal,
                    float *const restrict out) {
  int status = 0;
  ptrdiff_t const extended_size = 2 * n;
  _Complex float *fft_in = fftwf_malloc(sizeof(_Complex float) * extended_size);
  if (fft_in == NULL) {
    status = -1;
    goto cleanup;
  }
  _Complex float *fft_out =
      fftwf_malloc(sizeof(_Complex float) * extended_size);
  if (fft_out == NULL) {
    status = -1;
    goto cleanup;
  }
  // fprintf(stderr, "Copying...\n");

  for (ptrdiff_t i = 0; i < n; ++i) {
    fft_in[i] = signal[i];
  }
  for (ptrdiff_t i = n; i < extended_size; ++i) {
    fft_in[i] = 0;
  }
  // fprintf(stderr, "FFT ...\n");
  fft(extended_size, fft_in, fft_out);

  // fprintf(stderr, "Squaring ...\n");
  for (ptrdiff_t i = 0; i < extended_size; ++i) {
    fft_in[i] = cabsf(fft_out[i]);
  }
  // fprintf(stderr, "iFFT ...\n");
  ifft(extended_size, fft_in, fft_out);

  // fprintf(stderr, "Copying ...\n");
  for (ptrdiff_t i = 0; i < n; ++i) {
    out[i] = crealf(fft_out[i]) / crealf(fft_out[0]);
  }

cleanup:
  if (fft_out != NULL) {
    fftwf_free(fft_out);
  }
  if (fft_in != NULL) {
    fftwf_free(fft_in);
  }
  return status;
}
