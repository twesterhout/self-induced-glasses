#pragma once

#include <HsFFI.h>

void fft(HsInt n, _Complex float const *in, _Complex float *out);
void ifft(HsInt n, _Complex float const *in, _Complex float *out);
