//Generates a delta-impulse signal (1-D vector of floats),
//with specified length (N), samp, and amplitude.
//The signal is amp at sample samp, and 0 elsewhere.

//This is identical to unit_impulse, except that amplitude
//can be any float, rather than just 1.0.

//For the complex case, the real part is the delta_impulse,
//and the imaginary part is left at 0.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int delta_impulse_s (float *Y, const size_t N, const size_t samp, const float amp);
int delta_impulse_d (double *Y, const size_t N, const size_t samp, const double amp);
int delta_impulse_c (float *Y, const size_t N, const size_t samp, const float amp);
int delta_impulse_z (double *Y, const size_t N, const size_t samp, const double amp);


int delta_impulse_s (float *Y, const size_t N, const size_t samp, const float amp)
{
    if (samp+1u>=N) { fprintf(stderr, "error in delta_impulse_s: samp (delay) must be less than N-1\n"); return 1; }

    for (size_t n=samp; n>0u; --n, ++Y) { *Y = 0.0f; }
    *Y++ = amp;
    for (size_t n=N-samp-1u; n>0u; --n, ++Y) { *Y = 0.0f; }

    return 0;
}


int delta_impulse_d (double *Y, const size_t N, const size_t samp, const double amp)
{
    if (samp+1u>=N) { fprintf(stderr, "error in delta_impulse_d: samp (delay) must be less than N-1\n"); return 1; }

    for (size_t n=samp; n>0u; --n, ++Y) { *Y = 0.0; }
    *Y++ = amp;
    for (size_t n=N-samp-1u; n>0u; --n, ++Y) { *Y = 0.0; }

    return 0;
}


int delta_impulse_c (float *Y, const size_t N, const size_t samp, const float amp)
{
    if (samp+1u>=N) { fprintf(stderr, "error in delta_impulse_c: samp (delay) must be less than N-1\n"); return 1; }

    for (size_t n=2u*samp; n>0u; --n, ++Y) { *Y = 0.0f; }
    *Y++ = amp; *Y++ = 0.0f;
    for (size_t n=2u*(N-samp-1u); n>0u; --n, ++Y) { *Y = 0.0f; }

    return 0;
}


int delta_impulse_z (double *Y, const size_t N, const size_t samp, const double amp)
{
    if (samp+1u>=N) { fprintf(stderr, "error in delta_impulse_z: samp (delay) must be less than N-1\n"); return 1; }

    for (size_t n=2u*samp; n>0u; --n, ++Y) { *Y = 0.0; }
    *Y++ = amp; *Y++ = 0.0;
    for (size_t n=2u*(N-samp-1u); n>0u; --n, ++Y) { *Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
