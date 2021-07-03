//Gets random floats from normal distributions with mean=0 and std,
//takes fft, multiplies by 1/sqrt(f) amplitude, and takes ifft.
//Although this requires an extra fft compared to Octave,
//it requires half as many random numbers, and allows control of std.

//This uses modified code from PCG randoms minimal C library,
//but remains stand-alone (no install of PCG required).
//I use the textbook cos and sin formulae to make Gaussian.

#include <stdio.h>
#include <stdint.h>
#include <float.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int pink_s (float *Y, const size_t N, const float std, const char zmn);
int pink_d (double *Y, const size_t N, const double std, const char zmn);
int pink_c (float *Y, const size_t N, const float std, const char zmn);
int pink_z (double *Y, const size_t N, const double std, const char zmn);


int pink_s (float *Y, const size_t N, const float std, const char zmn)
{
    if (std<0.0f) { fprintf(stderr, "error in pink_s: std must be nonnegative\n"); return 1; }

    if (N==0u) {}
    else if (std<FLT_EPSILON)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0f; }
    }
    else
    {
        const float M_2PI = (float)(2.0*M_PI);
        float u1, u2, R;
        uint32_t r, xorshifted, rot;
        uint64_t state = 0u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
        struct timespec ts;

        //Init random num generator
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in pink_s: timespec_get.\n"); perror("timespec_get"); return 1; }
	    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

        //Generate
        if (std==1.0f)
        {
            for (size_t n=0u; n<N-1u; n+=2u)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((float)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((float)r,-32);
                R = sqrtf(-2.0f*logf(1.0f-u1));
                *Y++ = R * cosf(M_2PI*u2);
                *Y++ = R * sinf(M_2PI*u2);
            }
            if (N%2==1)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((float)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((float)r,-32);
                R = sqrtf(-2.0f*logf(1.0f-u1));
                *Y++ = R * cosf(M_2PI*u2);
            }
        }
        else
        {
            for (size_t n=0u; n<N-1u; n+=2u)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((float)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((float)r,-32);
                R = std * sqrtf(-2.0f*logf(1.0f-u1));
                *Y++ = R * cosf(M_2PI*u2);
                *Y++ = R * sinf(M_2PI*u2);
            }
            if (N%2==1)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((float)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((float)r,-32);
                R = std * sqrtf(-2.0f*logf(1.0f-u1));
                *Y++ = R * cosf(M_2PI*u2);
            }
        }
    }

    if (zmn)
    {
        float sm = 0.0f;
        for (size_t n=0; n<N; ++n) { sm += *--Y; }
        sm /= (float)N;
        for (size_t n=0; n<N; ++n, ++Y) { *Y -= sm; }
    }

    return 0;
}


int pink_d (double *Y, const size_t N, const double std, const char zmn)
{
    if (std<0.0) { fprintf(stderr, "error in pink_d: std must be nonnegative\n"); return 1; }

    if (N==0u) {}
    else if (std<DBL_EPSILON)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0; }
    }
    else
    {
        const double M_2PI = 2.0*M_PI;
        double u1, u2, R;
        uint32_t r, xorshifted, rot;
        uint64_t state = 0u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
        struct timespec ts;

        //Init random num generator
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in pink_d: timespec_get.\n"); perror("timespec_get"); return 1; }
	    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

        //Generate
        if (std==1.0)
        {
            for (size_t n=0u; n<N-1u; n+=2u)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((double)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((double)r,-32);
                R = sqrt(-2.0*log(1.0-u1));
                *Y++ = R * cos(M_2PI*u2);
                *Y++ = R * sin(M_2PI*u2);
            }
            if (N%2==1)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((double)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((double)r,-32);
                R = sqrt(-2.0*log(1.0-u1));
                *Y++ = R * cos(M_2PI*u2);
            }
        }
        else
        {
            for (size_t n=0u; n<N-1u; n+=2u)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((double)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((double)r,-32);
                R = std * sqrt(-2.0*log(1.0-u1));
                *Y++ = R * cos(M_2PI*u2);
                *Y++ = R * sin(M_2PI*u2);
            }
            if (N%2==1)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((double)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((double)r,-32);
                R = std * sqrt(-2.0*log(1.0-u1));
                *Y++ = R * cos(M_2PI*u2);
            }
        }
    }

    if (zmn)
    {
        double sm = 0.0;
        for (size_t n=0; n<N; ++n) { sm += *--Y; }
        sm /= (double)N;
        for (size_t n=0; n<N; ++n, ++Y) { *Y -= sm; }
    }

    return 0;
}


int pink_c (float *Y, const size_t N, const float std, const char zmn)
{
    if (std<0.0f) { fprintf(stderr, "error in pink_s: std must be nonnegative\n"); return 1; }

    if (N==0u) {}
    else if (std<FLT_EPSILON)
    {
        for (size_t n=0u; n<2u*N; ++n, ++Y) { *Y = 0.0f; }
    }
    else
    {
        const float M_2PI = (float)(2.0*M_PI);
        float u1, u2, R;
        uint32_t r, xorshifted, rot;
        uint64_t state = 0u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
        struct timespec ts;

        //Init random num generator
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in pink_c: timespec_get.\n"); perror("timespec_get"); return 1; }
	    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

        //Generate
        if (std==1.0f)
        {
            for (size_t n=0u; n<N; ++n)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((float)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((float)r,-32);
                R = sqrtf(-2.0f*logf(1.0f-u1));
                *Y++ = R * cosf(M_2PI*u2);
                *Y++ = R * sinf(M_2PI*u2);
            }
        }
        else
        {
            for (size_t n=0u; n<N; ++n)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((float)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((float)r,-32);
                R = std * sqrtf(-2.0f*logf(1.0f-u1));
                *Y++ = R * cosf(M_2PI*u2);
                *Y++ = R * sinf(M_2PI*u2);
            }
        }
    }

    if (zmn)
    {
        float smr=0.0f, smi=0.0f;
        for (size_t n=0; n<N; ++n) { smi += *--Y; smr += *--Y; }
        smr /= (float)N; smi /= (float)N;
        for (size_t n=0; n<N; ++n) { *Y++ -= smr; *Y++ -= smi; }
    }

    return 0;
}


int pink_z (double *Y, const size_t N, const double std, const char zmn)
{
    if (std<0.0) { fprintf(stderr, "error in pink_z: std must be nonnegative\n"); return 1; }

    if (N==0u) {}
    else if (std<DBL_EPSILON)
    {
        for (size_t n=0u; n<2u*N; ++n, ++Y) { *Y = 0.0; }
    }
    else
    {
        const double M_2PI = 2.0*M_PI;
        double u1, u2, R;
        uint32_t r, xorshifted, rot;
        uint64_t state = 0u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
        struct timespec ts;

        //Init random num generator
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in pink_z: timespec_get.\n"); perror("timespec_get"); return 1; }
	    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

        //Generate
        if (std==1.0)
        {
            for (size_t n=0u; n<N; ++n)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((double)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((double)r,-32);
                R = sqrt(-2.0*log(1.0-u1));
                *Y++ = R * cos(M_2PI*u2);
                *Y++ = R * sin(M_2PI*u2);
            }
        }
        else
        {
            for (size_t n=0u; n<N; ++n)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((double)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((double)r,-32);
                R = std * sqrt(-2.0*log(1.0-u1));
                *Y++ = R * cos(M_2PI*u2);
                *Y++ = R * sin(M_2PI*u2);
            }
        }
    }

    if (zmn)
    {
        double smr=0.0, smi = 0.0;
        for (size_t n=0; n<N; ++n) { smi += *--Y; smr += *--Y; }
        smr /= (double)N; smi /= (double)N;
        for (size_t n=0; n<N; ++n) { *Y++ -= smr; *Y++ -= smi; }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
