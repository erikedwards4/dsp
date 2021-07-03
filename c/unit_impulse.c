//Generates a unit-impulse signal (1-D vector of floats),
//with specified length (N), delay (tau), and amplitude.
//The signal is amp at sample tau, and 0 elsewhere.

//This is identical to delta_impulse, except that the
//amplitude is fixed at 1.0 (since "unit").

//For the complex case, the real part is the unit_impulse,
//and the imaginary part is left at 0.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int unit_impulse_s (float *Y, const size_t N, const size_t tau);
int unit_impulse_d (double *Y, const size_t N, const size_t tau);
int unit_impulse_c (float *Y, const size_t N, const size_t tau);
int unit_impulse_z (double *Y, const size_t N, const size_t tau);


int unit_impulse_s (float *Y, const size_t N, const size_t tau)
{
    if (tau>=N) { fprintf(stderr, "error in unit_impulse_s: tau (delay) must be less than N\n"); return 1; }

    for (size_t n=0u; n<tau; ++n, ++Y) { *Y = 0.0f; }
    *Y++ = 1.0f;
    for (size_t n=tau+1; n<N; ++n, ++Y) { *Y = 0.0f; }

    return 0;
}


int unit_impulse_d (double *Y, const size_t N, const size_t tau)
{
    if (tau>=N) { fprintf(stderr, "error in unit_impulse_d: tau (delay) must be less than N\n"); return 1; }

    for (size_t n=0u; n<tau; ++n, ++Y) { *Y = 0.0; }
    *Y++ = 1.0;
    for (size_t n=tau+1; n<N; ++n, ++Y) { *Y = 0.0; }

    return 0;
}


int unit_impulse_c (float *Y, const size_t N, const size_t tau)
{
    if (tau>=N) { fprintf(stderr, "error in unit_impulse_c: tau (delay) must be less than N\n"); return 1; }

    for (size_t n=0u; n<2u*tau; ++n, ++Y) { *Y = 0.0f; }
    *Y++ = 1.0f;
    for (size_t n=2u*tau+1; n<2u*N; ++n, ++Y) { *Y = 0.0f; }

    return 0;
}


int unit_impulse_z (double *Y, const size_t N, const size_t tau)
{
    if (tau>=N) { fprintf(stderr, "error in unit_impulse_z: tau (delay) must be less than N\n"); return 1; }

    for (size_t n=0u; n<2u*tau; ++n, ++Y) { *Y = 0.0; }
    *Y++ = 1.0;
    for (size_t n=2u*tau+1; n<2u*N; ++n, ++Y) { *Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
