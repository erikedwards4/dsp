//Generates a sinewave signal (1-D vector of floats),
//with specified length (N), amp, freq, and phase.
//Sine phase is used, such that signal is 0 at phase 0.

//For the complex case, the imaginary part is the
//corresponding cosinewave signal, thus giving a cisoid
//(but see cosinewave for other phase convention).

#include <stdio.h>
#include <float.h>
#include <math.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int sinewave_s (float *Y, const size_t N, const float amp, const float frq, const float phs);
int sinewave_d (double *Y, const size_t N, const double amp, const double frq, const double phs);
int sinewave_c (float *Y, const size_t N, const float amp, const float frq, const float phs);
int sinewave_z (double *Y, const size_t N, const double amp, const double frq, const double phs);


int sinewave_s (float *Y, const size_t N, const float amp, const float frq, const float phs)
{
    if (amp<0.0f) { fprintf(stderr, "error in sinewave_s: amp must be nonnegative\n"); return 1; }
    if (frq<FLT_EPSILON) { fprintf(stderr, "error in sinewave_s: freq must be positive\n"); return 1; }

    if (amp<FLT_EPSILON)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0f; }
    }
    else
    {
        const float f2pi = (float)(2.0*M_PI) * frq;

        for (size_t n=0u; n<N; ++n, ++Y)
        {
            *Y = amp * sinf(fmaf((float)n,f2pi,phs));
        }
    }

    return 0;
}


int sinewave_d (double *Y, const size_t N, const double amp, const double frq, const double phs)
{
    if (amp<0.0) { fprintf(stderr, "error in sinewave_d: amp must be nonnegative\n"); return 1; }
    if (frq<DBL_EPSILON) { fprintf(stderr, "error in sinewave_s: freq must be positive\n"); return 1; }

    if (amp<DBL_EPSILON)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0; }
    }
    else
    {
        const double f2pi = 2.0 * M_PI * frq;

        for (size_t n=0u; n<N; ++n, ++Y)
        {
            *Y = amp * sin(fma((double)n,f2pi,phs));
        }
    }

    return 0;
}


int sinewave_c (float *Y, const size_t N, const float amp, const float frq, const float phs)
{
    if (amp<0.0f) { fprintf(stderr, "error in sinewave_c: amp must be nonnegative\n"); return 1; }
    if (frq<FLT_EPSILON) { fprintf(stderr, "error in sinewave_s: freq must be positive\n"); return 1; }

    const float f2pi = (float)(2.0*M_PI) * frq;

    for (size_t n=0u; n<N; ++n)
    {
        *Y++ = amp * sinf(fmaf((float)n,f2pi,phs));
        *Y++ = amp * cosf(fmaf((float)n,f2pi,phs));
    }

    return 0;
}


int sinewave_z (double *Y, const size_t N, const double amp, const double frq, const double phs)
{
    if (amp<0.0) { fprintf(stderr, "error in sinewave_z: amp must be nonnegative\n"); return 1; }
    if (frq<DBL_EPSILON) { fprintf(stderr, "error in sinewave_z: freq must be positive\n"); return 1; }

    if (amp<DBL_EPSILON)
    {
        for (size_t n=0u; n<2u*N; ++n, ++Y) { *Y = 0.0; }
    }
    else
    {
        const double f2pi = 2.0 * M_PI * frq;

        for (size_t n=0u; n<N; ++n)
        {
            *Y++ = amp * sin(fma((double)n,f2pi,phs));
            *Y++ = amp * cos(fma((double)n,f2pi,phs));
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
