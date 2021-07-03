//Generates a sawtooth-wave signal (1-D vector of floats),
//with specified length (N), amp, freq, and phase.
//Sine phase is used, such that signal is 0 at phase 0.

#include <stdio.h>
#include <float.h>
#include <math.h>

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

#ifndef M_1_PI
    #define M_1_PI 0.318309886183790671538
#endif

#ifndef M_1_PIf
    #define M_1_PIf 0.318309886183790671538f
#endif


#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int sawwave_s (float *Y, const size_t N, const float amp, const float frq, const float phs);
int sawwave_d (double *Y, const size_t N, const double amp, const double frq, const double phs);
int sawwave_c (float *Y, const size_t N, const float amp, const float frq, const float phs);
int sawwave_z (double *Y, const size_t N, const double amp, const double frq, const double phs);


int sawwave_s (float *Y, const size_t N, const float amp, const float frq, const float phs)
{
    if (amp<0.0f) { fprintf(stderr, "error in sawwave_s: amp must be nonnegative\n"); return 1; }
    if (frq<FLT_EPSILON) { fprintf(stderr, "error in sawwave_s: freq must be positive\n"); return 1; }

    if (amp<FLT_EPSILON)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0f; }
    }
    else if (amp==1.0f)
    {
        const float M_2PI = 2.0f * M_PIf;
        const float f2pi = frq * M_2PI;
        float arg;

        for (size_t n=0; n<N; ++n, ++Y)
        {
            arg = fmodf(fmaf((float)n,f2pi,phs),M_2PI);
            if (arg>M_PIf && arg<M_2PI) { *Y = M_1_PIf*arg - 2.0f; }
            else if (arg>0.0f && arg<M_PIf) { *Y = M_1_PIf * arg; }
            else { *Y = 0.0f; }
        }
    }
    else
    {
        const float M_2PI = 2.0f * M_PIf;
        const float f2pi = frq * M_2PI;
        const float a2pi = amp * M_1_PIf;
        float arg;

        for (size_t n=0; n<N; ++n, ++Y)
        {
            arg = fmodf(fmaf((float)n,f2pi,phs),M_2PI);
            if (arg>M_PIf && arg<M_2PI) { *Y = a2pi*arg - 2.0f*amp; }
            else if (arg>0.0f && arg<M_PIf) { *Y = a2pi * arg; }
            else { *Y = 0.0f; }
        }
    }

    return 0;
}


int sawwave_d (double *Y, const size_t N, const double amp, const double frq, const double phs)
{
    if (amp<0.0) { fprintf(stderr, "error in sawwave_d: amp must be nonnegative\n"); return 1; }
    if (frq<DBL_EPSILON) { fprintf(stderr, "error in sawwave_d: freq must be positive\n"); return 1; }

    if (amp<DBL_EPSILON)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0; }
    }
    else if (amp==1.0)
    {
        const double M_2PI = 2.0 * M_PI;
        const double f2pi = M_2PI * frq;
        double arg;

        for (size_t n=0; n<N; ++n, ++Y)
        {
            arg = fmod(fma((double)n,f2pi,phs),M_2PI);
            if (arg>M_PI && arg<M_2PI) { *Y = M_1_PI*arg - 2.0; }
            else if (arg>0.0 && arg<M_PI) { *Y = M_1_PI * arg; }
            else { *Y = 0.0; }
        }
    }
    else
    {
        const double M_2PI = 2.0 * M_PI;
        const double f2pi = frq * M_2PI;
        const double a2pi = amp * M_1_PI;
        double arg;

        for (size_t n=0; n<N; ++n, ++Y)
        {
            arg = fmod(fma((double)n,f2pi,phs),M_2PI);
            if (arg>M_PI && arg<M_2PI) { *Y = a2pi*arg - 2.0*amp; }
            else if (arg>0.0 && arg<M_PI) { *Y = a2pi * arg; }
            else { *Y = 0.0; }
        }
    }

    return 0;
}


int sawwave_c (float *Y, const size_t N, const float amp, const float frq, const float phs)
{
    if (amp<0.0f) { fprintf(stderr, "error in sawwave_c: amp must be nonnegative\n"); return 1; }
    if (frq<FLT_EPSILON) { fprintf(stderr, "error in sawwave_c: freq must be positive\n"); return 1; }

    if (amp<FLT_EPSILON)
    {
        for (size_t n=0u; n<2u*N; ++n, ++Y) { *Y = 0.0f; }
    }
    else
    {
        const float M_2PI = 2.0f * M_PIf;
        const float M_3PI_2 = 3.0f * M_PIf / 2.0f;
        const float f2pi = frq * M_2PI;
        const float a2pi = amp * M_1_PIf;
        float arg;

        for (size_t n=0; n<N; ++n)
        {
            arg = fmodf(fmaf((float)n,f2pi,phs),M_2PI);

            if (arg<M_PI_2f) { *Y++ = a2pi*arg + 0.5f; }
            else if (arg>M_PI_2f && arg<M_3PI_2) { *Y++ = a2pi*arg - 1.5f; }
            else if (arg>M_3PI_2) { *Y++ = a2pi*arg - 2.5f; }
            else { *Y++ = 0.0f; }

            if (arg>M_PIf) { *Y++ = a2pi*arg - 2.0f*amp; }
            else if (arg>0.0f && arg<M_PIf) { *Y++ = a2pi * arg; }
            else { *Y++ = 0.0f; }
        }
    }

    return 0;
}


int sawwave_z (double *Y, const size_t N, const double amp, const double frq, const double phs)
{
    if (amp<0.0) { fprintf(stderr, "error in sawwave_z: amp must be nonnegative\n"); return 1; }
    if (frq<DBL_EPSILON) { fprintf(stderr, "error in sawwave_z: freq must be positive\n"); return 1; }

    if (amp<DBL_EPSILON)
    {
        for (size_t n=0u; n<2u*N; ++n, ++Y) { *Y = 0.0; }
    }
    else
    {
        const double M_2PI = 2.0 * M_PI;
        const double M_3PI_2 = 3.0 * M_PI / 2.0;
        const double f2pi = frq * M_2PI;
        const double a2pi = amp * M_1_PI;
        double arg;

        for (size_t n=0; n<N; ++n)
        {
            arg = fmod(fma((double)n,f2pi,phs),M_2PI);

            if (arg<M_PI_2) { *Y++ = a2pi*arg + 0.5; }
            else if (arg>M_PI_2 && arg<M_3PI_2) { *Y++ = a2pi*arg - 1.5; }
            else if (arg>M_3PI_2) { *Y++ = a2pi*arg - 2.5; }
            else { *Y++ = 0.0; }

            if (arg>M_PI) { *Y++ = a2pi*arg - 2.0*amp; }
            else if (arg>0.0 && arg<M_PI) { *Y++ = a2pi * arg; }
            else { *Y++ = 0.0; }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
