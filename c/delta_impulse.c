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

int delta_impulse_s (float *Y, const size_t N, const size_t tau, const float amp);
int delta_impulse_d (double *Y, const size_t N, const size_t tau, const float amp);
int delta_impulse_c (float *Y, const size_t N, const size_t tau, const float amp);
int delta_impulse_z (double *Y, const size_t N, const size_t tau, const float amp);


int delta_impulse_s (float *Y, const size_t N, const size_t tau, const float amp)
{
    if (tau>=N) { fprintf(stderr, "error in delta_impulse_s: tau (delay) must be less than N\n"); return 1; }

    for (size_t n=0u; n<tau; ++n, ++Y) { *Y = 0.0f; }
    *Y++ = amp;
    for (size_t n=tau+1; n<N; ++n, ++Y) { *Y = 0.0f; }

    return 0;
}


int delta_impulse_d (double *Y, const size_t N, const size_t tau, const float amp)
{
    if (tau>=N) { fprintf(stderr, "error in delta_impulse_d: tau (delay) must be less than N\n"); return 1; }

    for (size_t n=0u; n<tau; ++n, ++Y) { *Y = 0.0; }
    *Y++ = amp;
    for (size_t n=tau+1; n<N; ++n, ++Y) { *Y = 0.0; }

    return 0;
}


int delta_impulse_c (float *Y, const size_t N, const size_t tau, const float amp)
{
    if (tau>=N) { fprintf(stderr, "error in delta_impulse_c: tau (delay) must be less than N\n"); return 1; }

    for (size_t n=0u; n<2u*tau; ++n, ++Y) { *Y = 0.0f; }
    *Y++ = amp;
    for (size_t n=2u*tau+1; n<2u*N; ++n, ++Y) { *Y = 0.0f; }

    return 0;
}


int delta_impulse_z (double *Y, const size_t N, const size_t tau, const float amp)
{
    if (tau>=N) { fprintf(stderr, "error in delta_impulse_z: tau (delay) must be less than N\n"); return 1; }

    for (size_t n=0u; n<2u*tau; ++n, ++Y) { *Y = 0.0; }
    *Y++ = amp;
    for (size_t n=2u*tau+1; n<2u*N; ++n, ++Y) { *Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
