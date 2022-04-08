//Gets the first F frequencies for an STFT of order nfft,
//assuming a sample rate of fs.
//The frequencies are just linearly spaced from 0 to Nyquist (fs/2).
//If F>nfft/2+1, then the negative frequencies are included.

#include <stdio.h>
#include "codee_dsp.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int get_stft_freqs_s (float *Y, const size_t nfft, const size_t F, const float fs);
int get_stft_freqs_d (double *Y, const size_t nfft, const size_t F, const double fs);


int get_stft_freqs_s (float *Y, const size_t nfft, const size_t F, const float fs)
{
    if (F>nfft) { fprintf(stderr,"error in get_stft_freqs_s: F must be <= nfft\n"); return 1; }

    const size_t Nyq = nfft/2u + 1u;
    const float finc = fs/(float)nfft;

    *Y++ = 0.0f;
    for (size_t f=1u; f<F && f<Nyq; ++f, ++Y) { *Y = (float)f * finc; }
    for (size_t f=Nyq; f<F; ++f, ++Y) { *Y = -finc * (float)(nfft-f); }

    return 0;
}


int get_stft_freqs_d (double *Y, const size_t nfft, const size_t F, const double fs)
{
    if (F>nfft) { fprintf(stderr,"error in get_stft_freqs_d: F must be <= nfft\n"); return 1; }

    const size_t Nyq = nfft/2u + 1u;
    const double finc = fs/(double)nfft;

    *Y++ = 0.0;
    for (size_t f=1u; f<F && f<Nyq; ++f, ++Y) { *Y = (double)f * finc; }
    for (size_t f=Nyq; f<F; ++f, ++Y) { *Y = -finc * (double)(nfft-f); }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
