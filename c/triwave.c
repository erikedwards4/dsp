//Generates a triangle-wave signal (1-D vector of floats),
//with specified length (N), amp, freq, and phase.
//Sine phase is used, such that signal is 0 at phase 0.

//For the complex case, the real part is the cosine-phase triwave,
//and the imaginary part is the corresponding sine-phase triwave.

#include <stdio.h>
#include <float.h>
#include <math.h>
#include "codee_dsp.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifndef M_PIf
    #define M_PIf 3.14159265358979323846f
#endif

#ifndef M_PI_2
    #define M_PI_2 1.57079632679489661923
#endif

#ifndef M_PI_2f
    #define M_PI_2f 1.57079632679489661923f
#endif

#ifndef M_2_PI
    #define M_2_PI 0.636619772367581343076
#endif

#ifndef M_2_PIf
    #define M_2_PIf 0.636619772367581343076f
#endif


#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int triwave_s (float *Y, const size_t N, const float amp, const float frq, const float phs)
{
    if (amp<0.0f) { fprintf(stderr, "error in triwave_s: amp must be nonnegative\n"); return 1; }
    if (frq<FLT_EPSILON) { fprintf(stderr, "error in triwave_s: freq must be positive\n"); return 1; }

    if (N==0u) {}
    else if (amp<FLT_EPSILON)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (amp==1.0f)
    {
        const float M_2PI = (float)(2.0*M_PI);
        const float M_3PI_2 = (float)(3.0*M_PI/2.0);
        const float f2pi = frq * M_2PI;
        float arg;

        for (size_t n=0u; n<N; ++n, ++Y)
        {
            arg = fmodf(fmaf((float)n,f2pi,phs),M_2PI);
            if (arg>M_3PI_2 && arg<M_2PI) { *Y = M_2_PIf*arg - 4.0f; }
            else if (arg>M_PI_2f) { *Y = 2.0f - M_2_PIf*arg; }
            else if (arg>0.0f) { *Y = M_2_PIf * arg; }
            else { *Y = 0.0f; }
        }
    }
    else
    {
        const float M_2PI = (float)(2.0*M_PI);
        const float M_3PI_2 = (float)(3.0*M_PI/2.0);
        const float f2pi = frq * M_2PI;
        const float a2pi = amp * M_2_PIf;
        float arg;

        for (size_t n=0u; n<N; ++n, ++Y)
        {
            arg = fmodf(fmaf((float)n,f2pi,phs),M_2PI);
            if (arg>M_3PI_2 && arg<M_2PI) { *Y = a2pi*arg - 4.0f*amp; }
            else if (arg>M_PI_2f) { *Y = 2.0f*amp - a2pi*arg; }
            else if (arg>0.0f) { *Y = a2pi * arg; }
            else { *Y = 0.0f; }
        }
    }

    return 0;
}


int triwave_d (double *Y, const size_t N, const double amp, const double frq, const double phs)
{
    if (amp<0.0) { fprintf(stderr, "error in triwave_d: amp must be nonnegative\n"); return 1; }
    if (frq<DBL_EPSILON) { fprintf(stderr, "error in triwave_d: freq must be positive\n"); return 1; }

    if (N==0u) {}
    else if (amp<DBL_EPSILON)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (amp==1.0)
    {
        const double M_2PI = 2.0 * M_PI;
        const double M_3PI_2 = 3.0 * M_PI / 2.0;
        const double f2pi = M_2PI * frq;
        double arg;

        for (size_t n=0u; n<N; ++n, ++Y)
        {
            arg = fmod(fma((double)n,f2pi,phs),M_2PI);
            if (arg>M_3PI_2 && arg<M_2PI) { *Y = M_2_PI*arg - 4.0; }
            else if (arg>M_PI_2) { *Y = 2.0 - M_2_PI*arg; }
            else if (arg>0.0) { *Y = M_2_PI * arg; }
            else { *Y = 0.0; }
        }
    }
    else
    {
        const double M_2PI = 2.0 * M_PI;
        const double M_3PI_2 = 3.0 * M_PI / 2.0;
        const double f2pi = M_2PI * frq;
        const double a2pi = amp * M_2_PI;
        double arg;

        for (size_t n=0u; n<N; ++n, ++Y)
        {
            arg = fmod(fma((double)n,f2pi,phs),M_2PI);
            if (arg>M_3PI_2 && arg<M_2PI) { *Y = a2pi*arg - 4.0*amp; }
            else if (arg>M_PI_2) { *Y = 2.0*amp - a2pi*arg; }
            else if (arg>0.0) { *Y = a2pi * arg; }
            else { *Y = 0.0; }
        }
    }

    return 0;
}


int triwave_c (float *Y, const size_t N, const float amp, const float frq, const float phs)
{
    if (amp<0.0f) { fprintf(stderr, "error in triwave_c: amp must be nonnegative\n"); return 1; }
    if (frq<FLT_EPSILON) { fprintf(stderr, "error in triwave_c: freq must be positive\n"); return 1; }

    if (N==0u) {}
    else if (amp<FLT_EPSILON)
    {
        for (size_t n=2u*N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else
    {
        const float M_2PI = (float)(2.0*M_PI);
        const float M_3PI_2 = (float)(3.0*M_PI/2.0);
        const float f2pi = frq * M_2PI;
        const float a2pi = amp * M_2_PIf;
        float arg;

        for (size_t n=0u; n<N; ++n)
        {
            arg = fmodf(fmaf((float)n,f2pi,phs),M_2PI);

            if (arg>M_PIf) { *Y++ = a2pi*arg - 3.0f*amp; }
            else if (arg>0.0f) { *Y++ = amp - a2pi*arg; }
            else { *Y++ = amp; }

            if (arg>M_3PI_2 && arg<M_2PI) { *Y++ = a2pi*arg - 4.0f*amp; }
            else if (arg>M_PI_2f) { *Y++ = 2.0f*amp - a2pi*arg; }
            else if (arg>0.0f) { *Y++ = a2pi * arg; }
            else { *Y++ = 0.0f; }
        }
    }

    return 0;
}


int triwave_z (double *Y, const size_t N, const double amp, const double frq, const double phs)
{
    if (amp<0.0) { fprintf(stderr, "error in triwave_z: amp must be nonnegative\n"); return 1; }
    if (frq<DBL_EPSILON) { fprintf(stderr, "error in triwave_z: freq must be positive\n"); return 1; }

    if (N==0u) {}
    else if (amp<DBL_EPSILON)
    {
        for (size_t n=2u*N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else
    {
        const double M_2PI = 2.0 * M_PI;
        const double M_3PI_2 = 3.0 * M_PI / 2.0;
        const double f2pi = frq * M_2PI;
        const double a2pi = amp * M_2_PI;
        double arg;

        for (size_t n=0u; n<N; ++n)
        {
            arg = fmod(fma((double)n,f2pi,phs),M_2PI);

            if (arg>M_PI) { *Y++ = a2pi*arg - 3.0*amp; }
            else if (arg>0.0) { *Y++ = amp - a2pi*arg; }
            else { *Y++ = amp; }

            if (arg>M_3PI_2 && arg<M_2PI) { *Y++ = a2pi*arg - 4.0*amp; }
            else if (arg>M_PI_2) { *Y++ = 2.0*amp - a2pi*arg; }
            else if (arg>0.0) { *Y++ = a2pi * arg; }
            else { *Y++ = 0.0; }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
