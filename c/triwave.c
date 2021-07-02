//Generates a triangle-wave signal (1-D vector of floats),
//with specified length (N), amp, freq, and phase.
//Sine phase is used, such that signal is 0 at phase 0.

//For the complex case, the imaginary part is the
//corresponding cosine-phase signal.

#include <stdio.h>
#include <math.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
    #define M_PI_2 1.57079632679489661923
#endif

#ifndef M_PI_4
    #define M_PI_4 0.785398163397448309616
#endif


#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int triwave_s (float *Y, const size_t N, const float amp, const float frq, const float phs);
int triwave_d (double *Y, const size_t N, const double amp, const double frq, const double phs);
int triwave_c (float *Y, const size_t N, const float amp, const float frq, const float phs);
int triwave_z (double *Y, const size_t N, const double amp, const double frq, const double phs);


int triwave_s (float *Y, const size_t N, const float amp, const float frq, const float phs)
{
    if (amp<0.0f) { fprintf(stderr, "error in triwave_s: amp must be nonnegative\n"); return 1; }

    const float M_2PI = (float)(2.0*M_PI);
    const float f2pi = M_2PI * frq;
    float arg;

    for (size_t n=0; n<N; ++n, ++Y)
    {
        arg = fmodf(fmaf((float)n,f2pi,phs),M_2PI);
        if (arg>(float)M_PI) { *Y = -amp; }
        else if (arg>0.0f && arg<(float)M_PI) { *Y = amp; }
        //if (arg<0.0f) { *Y = -amp; }
        //else if (arg>0.0f) { *Y = amp; }
        else { *Y = 0.0f; }
    }

    return 0;
}


int triwave_d (double *Y, const size_t N, const double amp, const double frq, const double phs)
{
    if (amp<0.0) { fprintf(stderr, "error in triwave_d: amp must be nonnegative\n"); return 1; }

    const double M_2PI = 2.0 * M_PI;
    const double f2pi = M_2PI * frq;
    double arg;

    for (size_t n=0; n<N; ++n, ++Y)
    {
        arg = fmod(fma((double)n,f2pi,phs),M_2PI);
        if (arg>M_PI) { *Y = -amp; }
        else if (arg>0.0 && arg<M_PI) { *Y = amp; }
        else { *Y = 0.0; }
    }

    return 0;
}


int triwave_c (float *Y, const size_t N, const float amp, const float frq, const float phs)
{
    if (amp<0.0f) { fprintf(stderr, "error in triwave_c: amp must be nonnegative\n"); return 1; }

    const float M_2PI = (float)(2.0*M_PI);
    const float M_3PI_2 = (float)(M_PI+M_PI_2);
    const float f2pi = M_2PI * frq;
    float arg;

    for (size_t n=0; n<N; ++n)
    {
        arg = fmodf(fmaf((float)n,f2pi,phs),M_2PI);

        if (arg>(float)M_PI) { *Y++ = -amp; }
        else if (arg>0.0f && arg<(float)M_PI) { *Y++ = amp; }
        else { *Y++ = 0.0f; }

        if (arg<(float)M_PI_2) { *Y++ = amp; }
        else if (arg>(float)M_PI_2 && arg<M_3PI_2) { *Y++ = -amp; }
        else { *Y++ = 0.0f; }
    }

    return 0;
}


int triwave_z (double *Y, const size_t N, const double amp, const double frq, const double phs)
{
    if (amp<0.0) { fprintf(stderr, "error in triwave_z: amp must be nonnegative\n"); return 1; }

    const double M_2PI = 2.0 * M_PI;
    const double M_3PI_2 = M_PI + M_PI_2;
    const double f2pi = M_2PI * frq;
    double arg;

    for (size_t n=0; n<N; ++n)
    {
        arg = fmod(fma((double)n,f2pi,phs),M_2PI);

        if (arg>M_PI) { *Y++ = -amp; }
        else if (arg>0.0 && arg<M_PI) { *Y++ = amp; }
        else { *Y++ = 0.0; }

        if (arg<M_PI_2) { *Y++ = amp; }
        else if (arg>M_PI_2 && arg<M_3PI_2) { *Y++ = -amp; }
        else { *Y++ = 0.0; }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
