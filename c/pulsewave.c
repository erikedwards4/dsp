//Generates a pulse-wave signal (1-D vector of floats),
//with specified length (N), amp, freq, phase, and duty-cycle.
//Cosine phase is used, such that the signal is 1 at phase 0.

//For the complex case, the real part is the cosine-phase pulsewave,
//and the imaginary part is the corresponding sine-phase pulsewave.

#include <stdio.h>
#include <float.h>
#include <math.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
    #define M_PI_2 1.57079632679489661923
#endif

#ifndef M_PIf
    #define M_PIf 3.14159265358979323846f
#endif

#ifndef M_PI_2f
    #define M_PI_2f 1.57079632679489661923f
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int pulsewave_s (float *Y, const size_t N, const float amp, const float frq, const float phs, const float dty);
int pulsewave_d (double *Y, const size_t N, const double amp, const double frq, const double phs, const double dty);
int pulsewave_c (float *Y, const size_t N, const float amp, const float frq, const float phs, const float dty);
int pulsewave_z (double *Y, const size_t N, const double amp, const double frq, const double phs, const double dty);


int pulsewave_s (float *Y, const size_t N, const float amp, const float frq, const float phs, const float dty)
{
    if (amp<0.0f) { fprintf(stderr, "error in pulsewave_s: amp must be nonnegative\n"); return 1; }
    if (frq<FLT_EPSILON) { fprintf(stderr, "error in pulsewave_s: freq must be positive\n"); return 1; }

    if (amp<FLT_EPSILON)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0f; }
    }
    else
    {
        const float M_2PI = (float)(2.0*M_PI);
        const float f2pi = frq * M_2PI;
        float arg;

        for (size_t n=0; n<N; ++n, ++Y)
        {
            arg = fmodf(fmaf((float)n,f2pi,phs),M_2PI);
            if (arg>M_PIf && arg<M_2PI) { *Y = -amp; }
            else if (arg>0.0f && arg<M_PIf) { *Y = amp; }
            else { *Y = 0.0f; }
        }
    }

    return 0;
}


int pulsewave_d (double *Y, const size_t N, const double amp, const double frq, const double phs, const double dty)
{
    if (amp<0.0) { fprintf(stderr, "error in pulsewave_d: amp must be nonnegative\n"); return 1; }
    if (frq<DBL_EPSILON) { fprintf(stderr, "error in pulsewave_d: freq must be positive\n"); return 1; }

    if (amp<DBL_EPSILON)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0; }
    }
    else
    {
        const double M_2PI = 2.0 * M_PI;
        const double f2pi = frq * M_2PI;
        double arg;

        for (size_t n=0; n<N; ++n, ++Y)
        {
            arg = fmod(fma((double)n,f2pi,phs),M_2PI);
            if (arg>M_PI && arg<M_2PI) { *Y = -amp; }
            else if (arg>0.0 && arg<M_PI) { *Y = amp; }
            else { *Y = 0.0; }
        }
    }

    return 0;
}


int pulsewave_c (float *Y, const size_t N, const float amp, const float frq, const float phs, const float dty)
{
    if (amp<0.0f) { fprintf(stderr, "error in pulsewave_c: amp must be nonnegative\n"); return 1; }
    if (frq<FLT_EPSILON) { fprintf(stderr, "error in pulsewave_c: freq must be positive\n"); return 1; }

    if (amp<FLT_EPSILON)
    {
        for (size_t n=0u; n<2u*N; ++n, ++Y) { *Y = 0.0f; }
    }
    else
    {
        const float M_2PI = (float)(2.0*M_PI);
        const float M_3PI_2 = (float)(M_PI+M_PI_2);
        const float f2pi = frq * M_2PI;
        float arg;

        for (size_t n=0; n<N; ++n)
        {
            arg = fmodf(fmaf((float)n,f2pi,phs),M_2PI);

            if (arg<M_PI_2f || arg>M_3PI_2) { *Y++ = amp; }
            else if (arg>M_PI_2f && arg<M_3PI_2) { *Y++ = -amp; }
            else { *Y++ = 0.0f; }

            if (arg>M_PIf && arg<M_2PI) { *Y++ = -amp; }
            else if (arg>0.0f && arg<M_PIf) { *Y++ = amp; }
            else { *Y++ = 0.0f; }
        }
    }

    return 0;
}


int pulsewave_z (double *Y, const size_t N, const double amp, const double frq, const double phs, const double dty)
{
    if (amp<0.0) { fprintf(stderr, "error in pulsewave_z: amp must be nonnegative\n"); return 1; }
    if (frq<DBL_EPSILON) { fprintf(stderr, "error in pulsewave_z: freq must be positive\n"); return 1; }

    if (amp<DBL_EPSILON)
    {
        for (size_t n=0u; n<2u*N; ++n, ++Y) { *Y = 0.0; }
    }
    else
    {
        const double M_2PI = 2.0 * M_PI;
        const double M_3PI_2 = M_PI + M_PI_2;
        const double f2pi = frq * M_2PI;
        double arg;

        for (size_t n=0; n<N; ++n)
        {
            arg = fmod(fma((double)n,f2pi,phs),M_2PI);

            if (arg<M_PI_2 || arg>M_3PI_2) { *Y++ = amp; }
            else if (arg>M_PI_2 && arg<M_3PI_2) { *Y++ = -amp; }
            else { *Y++ = 0.0; }

            if (arg>M_PI && arg<M_2PI) { *Y++ = -amp; }
            else if (arg>0.0 && arg<M_PI) { *Y++ = amp; }
            else { *Y++ = 0.0; }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
