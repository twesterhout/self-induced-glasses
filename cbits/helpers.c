#include "helpers.h"
#include <cblas.h>
#include <complex.h>
#if USE_FFTW
#include <fftw3.h>
#else
#include "kiss_fft.h"
#endif
#include <stdlib.h>
#include <x86/avx2.h>

// float effective_interaction_loop(sig_model const *model, int const px,
//                                  int const py, int const radius) {
//
// }

static inline float read_spin(uint64_t const *const restrict x,
                              ptrdiff_t const i) {
  int const bit = (x[i / 64] >> (i % 64)) & 1;
  return 2 * bit - 1;
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

static inline float horizontal_add_128(simde__m128 const a) {
  simde__m128 t1 = simde_mm_movehl_ps(a, a);
  simde__m128 t2 = simde_mm_add_ps(a, t1);
  simde__m128 t3 = simde_mm_shuffle_ps(t2, t2, 1);
  simde__m128 t4 = simde_mm_add_ss(t2, t3);
  return simde_mm_cvtss_f32(t4);
}

static inline float horizontal_sum(simde__m256 const x) {
  return horizontal_add_128(simde_mm256_castps256_ps128(x)) +
         horizontal_add_128(simde_mm256_extractf128_ps(x, 1));
}

static inline float
packed_dot_product_scalar(ptrdiff_t const n, float const *const restrict a,
                          uint64_t const *const restrict b) {
  double acc = 0.0;
  for (ptrdiff_t i = 0; i < n; ++i) {
    acc += (double)a[i] * (double)read_spin(b, i);
  }
  return (float)acc;
}

static inline float packed_dot_product(ptrdiff_t const n,
                                       float const *const restrict a,
                                       uint64_t const *const restrict b) {
  // fprintf(stderr, "%s", "packed_dot_product...");
  simde__m256 sum_v = simde_mm256_setzero_ps();
#if USE_KAHAN
  simde__m256 err_v = simde_mm256_setzero_ps();
#endif
  ptrdiff_t const num_blocks = n / 64;
  for (ptrdiff_t block_idx = 0; block_idx < num_blocks; ++block_idx) {
    uint64_t word = b[block_idx];
    for (int k = 0; k < 8; ++k, word >>= 8) {
      uint8_t const byte = word & 0xFF;
      simde__m256 const mask = unpack_byte(byte);
      simde__m256 const not_mask =
          simde_mm256_castsi256_ps(simde_mm256_xor_si256(
              simde_mm256_castps_si256(mask), simde_mm256_set1_epi32(-1)));
      simde__m256 const a_v = simde_mm256_loadu_ps(a + 64 * block_idx + 8 * k);
#if USE_KAHAN
      simde__m256 t = simde_mm256_add_ps(err_v, simde_mm256_and_ps(a_v, mask));
      simde__m256 s = simde_mm256_add_ps(sum_v, t);
      simde__m256 z = simde_mm256_sub_ps(s, sum_v);
      err_v = simde_mm256_sub_ps(t, z);
      sum_v = s;

      t = simde_mm256_sub_ps(err_v, simde_mm256_and_ps(a_v, not_mask));
      s = simde_mm256_add_ps(sum_v, t);
      z = simde_mm256_sub_ps(s, sum_v);
      err_v = simde_mm256_sub_ps(t, z);
      sum_v = s;
#else
      sum_v = simde_mm256_add_ps(sum_v, simde_mm256_and_ps(a_v, mask));
      sum_v = simde_mm256_sub_ps(sum_v, simde_mm256_and_ps(a_v, not_mask));
#endif
    }
  }

  double acc = horizontal_sum(sum_v);
  for (ptrdiff_t i = 64 * num_blocks; i < n; ++i) {
    acc += (double)a[i] * (double)read_spin(b, i);
  }

#if CHECK_CORRECTNESS
  float const acc_ref = packed_dot_product_scalar(n, a, b);
  float const atol = 2e-6f;
  float const rtol = 1e-5f;
  if (fabs(acc - acc_ref) > atol + rtol * fabs(acc_ref)) {
    fprintf(stderr, "acc != acc_ref: %f vs. %f\n", acc, acc_ref);
    abort();
  }
#endif
  // fprintf(stderr, "%s\n", " done!");
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
  // fprintf(stderr, "%s", "total_energy...");
  float e = 0.0f;
  for (ptrdiff_t k = 0; k < n; ++k) {
    const float t = packed_dot_product(n, couplings + k * n, x);
    e += read_spin(x, k) * t;
  }
  // fprintf(stderr, "%s\n", " done!");
  return -e;
}

ptrdiff_t hamming_weight(ptrdiff_t const n, uint64_t const *const restrict x) {
  // fprintf(stderr, "%s", "hamming_weight...");
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
  // fprintf(stderr, "%s\n", " done!");
  return m;
}

void compute_structure_factor(HsInt const batch_size, HsInt const number_bits,
                              uint64_t const *const restrict input,
                              float *const restrict out) {
  memset(out, /*ch=*/0, /*count=*/number_bits * number_bits * sizeof(float));
  // Diagonal part
  for (ptrdiff_t i = 0; i < number_bits; ++i) {
    out[i * number_bits + i] = 1.0f;
  }
  // Off-diagonal part
  ptrdiff_t const number_words = (number_bits + 63) / 64;
  float const scale = 1.0f / batch_size;
  for (ptrdiff_t batch_idx = 0; batch_idx < batch_size; ++batch_idx) {
    uint64_t const *const x = input + batch_idx * number_words;
    for (ptrdiff_t i = 0; i < number_bits; ++i) {
      float const s_i = read_spin(x, i);
      for (ptrdiff_t j = i + 1; j < number_bits; ++j) {
        float const s_j = read_spin(x, j);
        float const term = scale * s_i * s_j;
        out[i * number_bits + j] += term;
        out[j * number_bits + i] += term;
      }
    }
  }
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

void extract_evolution(HsInt const batch_size, HsInt const number_words,
                       uint64_t const *const restrict input, HsInt const i,
                       float *const restrict out) {
  ptrdiff_t const word_index = i / 64;
  int const bit_index = i % 64;
  for (ptrdiff_t k = 0; k < batch_size; ++k) {
    uint64_t const word = input[word_index + k * number_words];
    int const bit = (word >> bit_index) & 1;
    out[k] = (float)(2 * bit - 1);
  }
}

void contiguous_axpy(HsInt const n, float const a,
                     float const *const restrict x, float *const restrict y) {
  cblas_saxpy(n, a, x, 1, y, 1);
}

static void fft(ptrdiff_t const n, _Complex float const *in,
                _Complex float *out) {
#if USE_FFTW
  // NOTE: the cast is relatively safe because we're using FFTW_ESTIMATE which
  // should not touch the arrays
  fftwf_plan p = fftwf_plan_dft_1d(n, (_Complex float *)in, out, FFTW_FORWARD,
                                   FFTW_ESTIMATE);
  fftwf_execute(p);
  fftwf_destroy_plan(p);
#else
  // fprintf(stderr, "Allocating...\n");
  kiss_fft_cfg cfg = kiss_fft_alloc(n, 0, NULL, NULL);
  // fprintf(stderr, "Running FFT...\n");
  kiss_fft(cfg, (kiss_fft_cpx const *)in, (kiss_fft_cpx *)out);
  // fprintf(stderr, "Freeing ...\n");
  kiss_fft_free(cfg);
#endif
}

static void ifft(ptrdiff_t const n, _Complex float const *in,
                 _Complex float *out) {
#if USE_FFTW
  // NOTE: the cast is relatively safe because we're using FFTW_ESTIMATE which
  // should not touch the arrays
  fftwf_plan p = fftwf_plan_dft_1d(n, (_Complex float *)in, out, FFTW_BACKWARD,
                                   FFTW_ESTIMATE);
  fftwf_execute(p);
  fftwf_destroy_plan(p);
#else
  // fprintf(stderr, "Allocating...\n");
  kiss_fft_cfg cfg = kiss_fft_alloc(n, 1, NULL, NULL);
  // fprintf(stderr, "Running iFFT...\n");
  kiss_fft(cfg, (kiss_fft_cpx const *)in, (kiss_fft_cpx *)out);
  // fprintf(stderr, "Freeing ...\n");
  kiss_fft_free(cfg);
#endif
}

int autocorrelation(ptrdiff_t const n, float const *const restrict signal,
                    float *const restrict out) {
  int status = 0;
  ptrdiff_t const extended_size = 2 * n;
  _Complex float *fft_in = NULL;
  _Complex float *fft_out = NULL;
  fft_in = malloc(sizeof(_Complex float) * extended_size);
  if (fft_in == NULL) {
    status = -1;
    goto cleanup;
  }
  fft_out = malloc(sizeof(_Complex float) * extended_size);
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
    free(fft_out);
  }
  if (fft_in != NULL) {
    free(fft_in);
  }
  return status;
}

static float compute_bits_overlap(ptrdiff_t const number_bits,
                                  uint64_t const *const restrict a,
                                  uint64_t const *const restrict b) {
  ptrdiff_t const num_blocks = number_bits / 64;
  ptrdiff_t const rest = number_bits - 64 * num_blocks;
  ptrdiff_t acc = 0;
  for (ptrdiff_t block_idx = 0; block_idx < num_blocks; ++block_idx) {
    uint64_t const word_a = a[block_idx];
    uint64_t const word_b = b[block_idx];
    acc += 64 - 2 * __builtin_popcountl(word_a ^ word_b);
  }
  if (rest != 0) {
    uint64_t const word_a = a[num_blocks];
    uint64_t const word_b = b[num_blocks];
    uint64_t const mask = (~(uint64_t)0) >> (64 - rest);
    acc += rest - 2 * __builtin_popcountl((word_a ^ word_b) & mask);
  }
  return (double)acc / (double)number_bits;
}

void two_point_autocorr_function(HsInt const n, HsInt const offset,
                                 uint64_t const *bits, HsInt const number_bits,
                                 float *out) {
  HsInt const number_words = (number_bits + 63) / 64;
  uint64_t const *a = bits + offset * number_words;
  for (ptrdiff_t i = 0; i < (n - offset); ++i) {
    out[i] = compute_bits_overlap(number_bits, a, a + i * number_words);
  }
}
