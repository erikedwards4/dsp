//Generates a cosinewave signal (1-D vector of floats),
//with specified length (N), amp, freq, and phase.
//Cosine phase is used, such that the signal is at its max at phase 0.
//Otherwise, this is identical to sinewave.

//For the complex case, the real part is the cosine-phase wave,
//and the imaginary part is the corresponding sine-phase wave,
//thus making a cisoid.

#include <stdio.h>
#include <float.h>
#include <math.h>
#include "codee_dsp.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int cosinewave_s (float *Y, const size_t N, const float amp, const float frq, const float phs)
{
    if (amp<0.0f) { fprintf(stderr, "error in cosinewave_s: amp must be nonnegative\n"); return 1; }
    if (frq<FLT_EPSILON) { fprintf(stderr, "error in cosinewave_s: freq must be positive\n"); return 1; }

    if (N==0u) {}
    else if (amp<FLT_EPSILON)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (amp==1.0f)
    {
        const float f2pi = (float)(2.0*M_PI) * frq;
        for (size_t n=0u; n<N; ++n, ++Y)
        {
            *Y = cosf(fmaf((float)n,f2pi,phs));
        }
    }
    else
    {
        const float f2pi = (float)(2.0*M_PI) * frq;
        for (size_t n=0u; n<N; ++n, ++Y)
        {
            *Y = amp * cosf(fmaf((float)n,f2pi,phs));
        }
    }

    return 0;
}


int cosinewave_d (double *Y, const size_t N, const double amp, const double frq, const double phs)
{
    if (amp<0.0) { fprintf(stderr, "error in cosinewave_d: amp must be nonnegative\n"); return 1; }
    if (frq<DBL_EPSILON) { fprintf(stderr, "error in cosinewave_d: freq must be positive\n"); return 1; }

    if (N==0u) {}
    else if (amp<DBL_EPSILON)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (amp==1.0)
    {
        const double f2pi = 2.0 * M_PI * frq;
        for (size_t n=0u; n<N; ++n, ++Y)
        {
            *Y = cos(fma((double)n,f2pi,phs));
        }
    }
    else
    {
        const double f2pi = 2.0 * M_PI * frq;
        for (size_t n=0u; n<N; ++n, ++Y)
        {
            *Y = amp * cos(fma((double)n,f2pi,phs));
        }
    }

    return 0;
}


int cosinewave_c (float *Y, const size_t N, const float amp, const float frq, const float phs)
{
    if (amp<0.0f) { fprintf(stderr, "error in cosinewave_c: amp must be nonnegative\n"); return 1; }
    if (frq<FLT_EPSILON) { fprintf(stderr, "error in cosinewave_c: freq must be positive\n"); return 1; }

    if (N==0u) {}
    else if (amp<FLT_EPSILON)
    {
        for (size_t n=2u*N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (amp==1.0f)
    {
        const float f2pi = (float)(2.0*M_PI) * frq;
        for (size_t n=0u; n<N; ++n)
        {
            *Y++ = cosf(fmaf((float)n,f2pi,phs));
            *Y++ = sinf(fmaf((float)n,f2pi,phs));
        }
    }
    else
    {
        const float f2pi = (float)(2.0*M_PI) * frq;
        for (size_t n=0u; n<N; ++n)
        {
            *Y++ = amp * cosf(fmaf((float)n,f2pi,phs));
            *Y++ = amp * sinf(fmaf((float)n,f2pi,phs));
        }
    }

    return 0;
}


int cosinewave_z (double *Y, const size_t N, const double amp, const double frq, const double phs)
{
    if (amp<0.0) { fprintf(stderr, "error in cosinewave_z: amp must be nonnegative\n"); return 1; }
    if (frq<DBL_EPSILON) { fprintf(stderr, "error in cosinewave_z: freq must be positive\n"); return 1; }

    if (N==0u) {}
    else if (amp<DBL_EPSILON)
    {
        for (size_t n=2u*N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (amp==1.0)
    {
        const double f2pi = 2.0 * M_PI * frq;
        for (size_t n=0u; n<N; ++n)
        {
            *Y++ = cos(fma((double)n,f2pi,phs));
            *Y++ = sin(fma((double)n,f2pi,phs));
        }
    }
    else
    {
        const double f2pi = 2.0 * M_PI * frq;
        for (size_t n=0u; n<N; ++n)
        {
            *Y++ = amp * cos(fma((double)n,f2pi,phs));
            *Y++ = amp * sin(fma((double)n,f2pi,phs));
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
