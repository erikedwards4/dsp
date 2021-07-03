//Generates a rectangular-pulse signal (1-D vector of floats),
//with specified length (N), delay (tau), width, and amplitude.
//The signal is amp from sample tau to tau+w, and 0 elsewhere.

//For a width of 1, this is identical to delta_impulse.

//For the complex case, the real part is the rectangular_pulse,
//and the imaginary part is left at 0.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int rectangular_pulse_s (float *Y, const size_t N, const size_t tau, const size_t width, const float amp);
int rectangular_pulse_d (double *Y, const size_t N, const size_t tau, const size_t width, const float amp);
int rectangular_pulse_c (float *Y, const size_t N, const size_t tau, const size_t width, const float amp);
int rectangular_pulse_z (double *Y, const size_t N, const size_t tau, const size_t width, const float amp);


int rectangular_pulse_s (float *Y, const size_t N, const size_t tau, const size_t width, const float amp)
{
    if (tau>=N) { fprintf(stderr, "error in rectangular_pulse_s: tau (delay) must be less than N\n"); return 1; }

    for (size_t n=0u; n<tau; ++n, ++Y) { *Y = 0.0f; }
    for (size_t n=tau+1; n<tau+width && n<N; ++n, ++Y) { *Y = amp; }
    for (size_t n=tau+width; n<N; ++n, ++Y) { *Y = 0.0f; }

    return 0;
}


int rectangular_pulse_d (double *Y, const size_t N, const size_t tau, const size_t width, const float amp)
{
    if (tau>=N) { fprintf(stderr, "error in rectangular_pulse_d: tau (delay) must be less than N\n"); return 1; }

    for (size_t n=0u; n<tau; ++n, ++Y) { *Y = 0.0; }
    for (size_t n=tau+1; n<tau+width; ++n, ++Y) { *Y = amp; }
    for (size_t n=tau+width; n<N; ++n, ++Y) { *Y = 0.0; }

    return 0;
}


int rectangular_pulse_c (float *Y, const size_t N, const size_t tau, const size_t width, const float amp)
{
    if (tau>=N) { fprintf(stderr, "error in rectangular_pulse_c: tau (delay) must be less than N\n"); return 1; }

    for (size_t n=0u; n<2u*tau; ++n, ++Y) { *Y = 0.0f; }
    for (size_t n=tau+1; n<tau+width; ++n) { *Y++ = amp; *Y++ = 0.0f; }
    for (size_t n=2u*(tau+width); n<2u*N; ++n, ++Y) { *Y = 0.0f; }

    return 0;
}


int rectangular_pulse_z (double *Y, const size_t N, const size_t tau, const size_t width, const float amp)
{
    if (tau>=N) { fprintf(stderr, "error in rectangular_pulse_z: tau (delay) must be less than N\n"); return 1; }

    for (size_t n=0u; n<2u*tau; ++n, ++Y) { *Y = 0.0; }
    for (size_t n=tau+1; n<tau+width; ++n) { *Y++ = amp; *Y++ = 0.0; }
    for (size_t n=2u*(tau+width); n<2u*N; ++n, ++Y) { *Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
