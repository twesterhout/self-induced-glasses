#include "helpers.h"
#include <cblas.h>
#include <complex.h>
#include <fftw3.h>
#include <stdlib.h>

// float effective_interaction_loop(sig_model const *model, int const px,
//                                  int const py, int const radius) {
//
// }

float energy_change_upon_flip(ptrdiff_t const n,
                              float const *const restrict couplings,
                              float const *const restrict x,
                              ptrdiff_t const i) {
  const float overlap = cblas_sdot(n, couplings + i * n, 1, x, 1);
  const float s = x[i];
  return 2 * overlap * s;
}

float total_energy(ptrdiff_t const n, float const *const restrict couplings,
                   float const *const restrict x) {
  float e = 0.0f;
  for (ptrdiff_t k = 0; k < n; ++k) {
    const float t = cblas_sdot(n, couplings + k * n, 1, x, 1);
    e += x[k] * t;
  }
  return -e;
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
