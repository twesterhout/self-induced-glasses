#include "fft.h"
#if USE_FFTW
#include <fftw3.h>
#else
#include "kiss_fft.h"
#endif

void fft(HsInt const n, _Complex float const *in, _Complex float *out) {
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

void ifft(HsInt const n, _Complex float const *in, _Complex float *out) {
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
