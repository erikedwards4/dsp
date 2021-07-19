//Generates a rectangular-pulse signal (1-D vector of floats),
//with specified length (N), samp, width, and amplitude.
//The signal is amp from sample samp to samp+w, and 0 elsewhere.

//For a width of 1, this is identical to delta_impulse.

//For the complex case, the real part is the rect_pulse,
//and the imaginary part is left at 0.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int rect_pulse_s (float *Y, const size_t N, const size_t samp, const size_t width, const float amp);
int rect_pulse_d (double *Y, const size_t N, const size_t samp, const size_t width, const double amp);
int rect_pulse_c (float *Y, const size_t N, const size_t samp, const size_t width, const float amp);
int rect_pulse_z (double *Y, const size_t N, const size_t samp, const size_t width, const double amp);


int rect_pulse_s (float *Y, const size_t N, const size_t samp, const size_t width, const float amp)
{
    if (samp>=N) { fprintf(stderr, "error in rect_pulse_s: samp (delay) must be less than N\n"); return 1; }

    for (size_t n=0u; n<samp; ++n, ++Y) { *Y = 0.0f; }
    for (size_t n=samp; n<samp+width && n<N; ++n, ++Y) { *Y = amp; }
    for (size_t n=samp+width; n<N; ++n, ++Y) { *Y = 0.0f; }

    return 0;
}


int rect_pulse_d (double *Y, const size_t N, const size_t samp, const size_t width, const double amp)
{
    if (samp>=N) { fprintf(stderr, "error in rect_pulse_d: samp (delay) must be less than N\n"); return 1; }

    for (size_t n=0u; n<samp; ++n, ++Y) { *Y = 0.0; }
    for (size_t n=samp; n<samp+width && n<N; ++n, ++Y) { *Y = amp; }
    for (size_t n=samp+width; n<N; ++n, ++Y) { *Y = 0.0; }

    return 0;
}


int rect_pulse_c (float *Y, const size_t N, const size_t samp, const size_t width, const float amp)
{
    if (samp>=N) { fprintf(stderr, "error in rect_pulse_c: samp (delay) must be less than N\n"); return 1; }

    for (size_t n=0u; n<2u*samp; ++n, ++Y) { *Y = 0.0f; }
    for (size_t n=samp; n<samp+width && n<N; ++n) { *Y++ = amp; *Y++ = 0.0f; }
    for (size_t n=2u*(samp+width); n<2u*N; ++n, ++Y) { *Y = 0.0f; }

    return 0;
}


int rect_pulse_z (double *Y, const size_t N, const size_t samp, const size_t width, const double amp)
{
    if (samp>=N) { fprintf(stderr, "error in rect_pulse_z: samp (delay) must be less than N\n"); return 1; }

    for (size_t n=0u; n<2u*samp; ++n, ++Y) { *Y = 0.0; }
    for (size_t n=samp; n<samp+width && n<N; ++n) { *Y++ = amp; *Y++ = 0.0; }
    for (size_t n=2u*(samp+width); n<2u*N; ++n, ++Y) { *Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
