//Generates a pulse-wave signal (1-D vector of floats),
//with specified length (N), amp, freq, phase, and duty-cycle,
//where the duty-cycle is a float in [0.0 1.0].
//Cosine phase is used, such that the signal is 1 at phase 0.

//Negative amplitudes are allowed, although that stretches the
//definition of a "pulse-wave".

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
    if (frq<FLT_EPSILON) { fprintf(stderr, "error in pulsewave_s: freq must be positive\n"); return 1; }
    if (dty<0.0f || dty>1.0f) { fprintf(stderr, "error in pulsewave_s: dty must be in [0 1]\n"); return 1; }

    if (N==0u) {}
    else if (dty<FLT_EPSILON || (amp<FLT_EPSILON && amp>-FLT_EPSILON))
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (dty>1.0f-FLT_EPSILON)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = amp; }
    }
    else
    {
        const float M_2PI = (float)(2.0*M_PI);
        const float f2pi = frq * M_2PI;
        const float dtypi = dty * (float)(M_PI);
        float arg;

        for (size_t n=0u; n<N; ++n, ++Y)
        {
            arg = fmodf(fmaf((float)n,f2pi,phs),M_2PI);
            if (arg>dtypi && arg<M_2PI-dtypi) { *Y = 0.0f; }
            else { *Y = amp; }
        }
    }

    return 0;
}


int pulsewave_d (double *Y, const size_t N, const double amp, const double frq, const double phs, const double dty)
{
    if (frq<DBL_EPSILON) { fprintf(stderr, "error in pulsewave_d: freq must be positive\n"); return 1; }
    if (dty<0.0 || dty>1.0) { fprintf(stderr, "error in pulsewave_d: dty must be in [0 1]\n"); return 1; }

    if (N==0u) {}
    else if (dty<DBL_EPSILON || (amp<DBL_EPSILON && amp>-DBL_EPSILON))
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (dty>1.0-DBL_EPSILON)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = amp; }
    }
    else
    {
        const double M_2PI = 2.0 * M_PI;
        const double f2pi = frq * M_2PI;
        const double dtypi = dty * M_PI;
        double arg;

        for (size_t n=0u; n<N; ++n, ++Y)
        {
            arg = fmod(fma((double)n,f2pi,phs),M_2PI);
            if (arg>dtypi && arg<M_2PI-dtypi) { *Y = 0.0; }
            else { *Y = amp; }
        }
    }

    return 0;
}


int pulsewave_c (float *Y, const size_t N, const float amp, const float frq, const float phs, const float dty)
{
    if (amp<0.0f) { fprintf(stderr, "error in pulsewave_c: amp must be nonnegative\n"); return 1; }
    if (frq<FLT_EPSILON) { fprintf(stderr, "error in pulsewave_c: freq must be positive\n"); return 1; }
    if (dty<0.0f || dty>1.0f) { fprintf(stderr, "error in pulsewave_c: dty must be in [0 1]\n"); return 1; }

    if (N==0u) {}
    else if (dty<FLT_EPSILON || (amp<FLT_EPSILON && amp>-FLT_EPSILON))
    {
        for (size_t n=2u*N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (dty>1.0f-FLT_EPSILON)
    {
        for (size_t n=2u*N; n>0u; --n, ++Y) { *Y = amp; }
    }
    else
    {
        const float M_2PI = (float)(2.0*M_PI);
        const float f2pi = frq * M_2PI;
        const float dtypi = dty * (float)(M_PI);
        float arg;

        for (size_t n=0u; n<N; ++n, ++Y)
        {
            arg = fmodf(fmaf((float)n,f2pi,phs),M_2PI);
            if (arg>dtypi && arg<M_2PI-dtypi) { *Y++ = 0.0f; }
            else { *Y++ = amp; }
            arg = fmodf((float)(M_PI_2)-arg,M_2PI);
            if (arg>dtypi && arg<M_2PI-dtypi) { *Y++ = 0.0f; }
            else { *Y++ = amp; }
        }
    }

    return 0;
}


int pulsewave_z (double *Y, const size_t N, const double amp, const double frq, const double phs, const double dty)
{
    if (amp<0.0) { fprintf(stderr, "error in pulsewave_z: amp must be nonnegative\n"); return 1; }
    if (frq<DBL_EPSILON) { fprintf(stderr, "error in pulsewave_z: freq must be positive\n"); return 1; }
    if (dty<0.0 || dty>1.0) { fprintf(stderr, "error in pulsewave_z: dty must be in [0 1]\n"); return 1; }

    if (N==0u) {}
    else if (dty<DBL_EPSILON || (amp<DBL_EPSILON && amp>-DBL_EPSILON))
    {
        for (size_t n=2u*N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (dty>1.0-DBL_EPSILON)
    {
        for (size_t n=2u*N; n>0u; --n, ++Y) { *Y = amp; }
    }
    else
    {
        const double M_2PI = 2.0 * M_PI;
        const double f2pi = frq * M_2PI;
        const double dtypi = dty * M_PI;
        double arg;

        for (size_t n=0u; n<N; ++n)
        {
            arg = fmod(fma((double)n,f2pi,phs),M_2PI);
            if (arg>dtypi && arg<M_2PI-dtypi) { *Y++ = 0.0; }
            else { *Y++ = amp; }
            arg = fmod(M_PI_2-arg,M_2PI);
            if (arg>dtypi && arg<M_2PI-dtypi) { *Y++ = 0.0; }
            else { *Y++ = amp; }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
