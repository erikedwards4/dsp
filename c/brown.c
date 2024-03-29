//Gets random floats from a normal distribution with mean=0 and std,
//and then integrates (cumsum) to generate the brown noise.

//This uses modified code from PCG randoms minimal C library,
//but remains stand-alone (no install of PCG required).
//I use the textbook cos and sin formulae to make Gaussian.

//Since numerical overflow is very possible for long signals
//I use a sign reversal of the current random float 
//to prevent this (only if the resulting sum would overflow).

//Brown noise is also called "Brownian" noise (as in Brownian motion),
//"random-walk" noise (since generated by a random-walk process),
//and "red" noise (since pink is between white and red).

#include <stdio.h>
#include <stdint.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include "codee_dsp.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int brown_s (float *Y, const size_t N, const float std, const int zmn)
{
    if (std<0.0f) { fprintf(stderr, "error in brown_s: std must be nonnegative\n"); return 1; }

    if (N==0u) {}
    else if (std<FLT_EPSILON)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else
    {
        const float M_2PI = (float)(2.0*M_PI);
        float u1, u2, R, sm = 0.0f;
        uint32_t r, xorshifted, rot;
        uint64_t state = 0u;
        const uint64_t mul = 6364136223846793005u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
        struct timespec ts;

        //Init random num generator
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in brown_s: timespec_get.\n"); perror("timespec_get"); return 1; }
	    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

        //Generate
        if (std==1.0f)
        {
            for (size_t n=N-1u; n>0u; n-=2u)
            {
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((float)r,-32);
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((float)r,-32);
                R = sqrtf(-2.0f*logf(1.0f-u1));
                //sm += R * cosf(M_2PI*u2); *Y++ = sm;
                //sm += R * sinf(M_2PI*u2); *Y++ = sm;
                //Idea to fix overflow
                sm += R * cosf(M_2PI*u2);
                if (isnan(sm) || !isfinite(sm)) { sm = *(Y-1) - R*cosf(M_2PI*u2); }
                *Y++ = sm;
                sm += R * sinf(M_2PI*u2);
                if (isnan(sm) || !isfinite(sm)) { sm = *(Y-1) - R*sinf(M_2PI*u2); }
                *Y++ = sm;
            }
            if (N%2u==1u)
            {
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((float)r,-32);
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((float)r,-32);
                R = sqrtf(-2.0f*logf(1.0f-u1));
                //sm += R * cosf(M_2PI*u2); *Y++ = sm;
                //Idea to fix overflow
                sm += R * cosf(M_2PI*u2);
                if (isnan(sm) || !isfinite(sm)) { sm = *(Y-1) - R*cosf(M_2PI*u2); }
                *Y++ = sm;
            }
        }
        else
        {
            for (size_t n=N-1u; n>0u; n-=2u)
            {
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((float)r,-32);
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((float)r,-32);
                R = std * sqrtf(-2.0f*logf(1.0f-u1));
                //sm += R * cosf(M_2PI*u2); *Y++ = sm;
                //sm += R * sinf(M_2PI*u2); *Y++ = sm;
                //Idea to fix overflow
                sm += R * cosf(M_2PI*u2);
                if (isnan(sm) || !isfinite(sm)) { sm = *(Y-1) - R*cosf(M_2PI*u2); }
                *Y++ = sm;
                sm += R * sinf(M_2PI*u2);
                if (isnan(sm) || !isfinite(sm)) { sm = *(Y-1) - R*sinf(M_2PI*u2); }
                *Y++ = sm;
            }
            if (N%2u==1u)
            {
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((float)r,-32);
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((float)r,-32);
                R = std * sqrtf(-2.0f*logf(1.0f-u1));
                //sm += R * cosf(M_2PI*u2); *Y++ = sm;
                //Idea to fix overflow
                sm += R * cosf(M_2PI*u2);
                if (isnan(sm) || !isfinite(sm)) { sm = *(Y-1) - R*cosf(M_2PI*u2); }
                *Y++ = sm;
            }
        }

        if (zmn)
        {
            sm = 0.0f;
            for (size_t n=N; n>0u; --n) { sm += *--Y; }
            sm /= (float)N;
            for (size_t n=N; n>0u; --n, ++Y) { *Y -= sm; }
        }
    }

    return 0;
}


int brown_d (double *Y, const size_t N, const double std, const int zmn)
{
    if (std<0.0) { fprintf(stderr, "error in brown_d: std must be nonnegative\n"); return 1; }

    if (N==0u) {}
    else if (std<DBL_EPSILON)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else
    {
        const double M_2PI = 2.0*M_PI;
        double u1, u2, R, sm = 0.0;
        uint32_t r, xorshifted, rot;
        uint64_t state = 0u;
        const uint64_t mul = 6364136223846793005u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
        struct timespec ts;

        //Init random num generator
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in brown_d: timespec_get.\n"); perror("timespec_get"); return 1; }
	    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

        //Generate
        if (std==1.0)
        {
            for (size_t n=N-1u; n>0u; n-=2u)
            {
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((double)r,-32);
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((double)r,-32);
                R = sqrt(-2.0*log(1.0-u1));
                //sm += R * cos(M_2PI*u2); *Y++ = sm;
                //sm += R * sin(M_2PI*u2); *Y++ = sm;
                //Idea to fix overflow
                sm += R * cos(M_2PI*u2);
                if (isnan(sm) || !isfinite(sm)) { sm = *(Y-1) - R*cos(M_2PI*u2); }
                *Y++ = sm;
                sm += R * sin(M_2PI*u2);
                if (isnan(sm) || !isfinite(sm)) { sm = *(Y-1) - R*sin(M_2PI*u2); }
                *Y++ = sm;
            }
            if (N%2u==1u)
            {
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((double)r,-32);
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((double)r,-32);
                R = sqrt(-2.0*log(1.0-u1));
                //sm += R * cos(M_2PI*u2); *Y++ = sm;
                //Idea to fix overflow
                sm += R * cos(M_2PI*u2);
                if (isnan(sm) || !isfinite(sm)) { sm = *(Y-1) - R*cos(M_2PI*u2); }
                *Y++ = sm;
            }
        }
        else
        {
            for (size_t n=N-1u; n>0u; n-=2u)
            {
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((double)r,-32);
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((double)r,-32);
                R = std * sqrt(-2.0*log(1.0-u1));
                //sm += R * cos(M_2PI*u2); *Y++ = sm;
                //sm += R * sin(M_2PI*u2); *Y++ = sm;
                //Idea to fix overflow
                sm += R * cos(M_2PI*u2);
                if (isnan(sm) || !isfinite(sm)) { sm = *(Y-1) - R*cos(M_2PI*u2); }
                *Y++ = sm;
                sm += R * sin(M_2PI*u2);
                if (isnan(sm) || !isfinite(sm)) { sm = *(Y-1) - R*sin(M_2PI*u2); }
                *Y++ = sm;
            }
            if (N%2u==1u)
            {
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((double)r,-32);
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((double)r,-32);
                R = std * sqrt(-2.0*log(1.0-u1));
                //sm += R * cos(M_2PI*u2); *Y++ = sm;
                //Idea to fix overflow
                sm += R * cos(M_2PI*u2);
                if (isnan(sm) || !isfinite(sm)) { sm = *(Y-1) - R*cos(M_2PI*u2); }
                *Y++ = sm;
            }
        }

        if (zmn)
        {
            sm = 0.0;
            for (size_t n=N; n>0u; --n) { sm += *--Y; }
            sm /= (double)N;
            for (size_t n=N; n>0u; --n, ++Y) { *Y -= sm; }
        }
    }

    return 0;
}


int brown_c (float *Y, const size_t N, const float std, const int zmn)
{
    if (std<0.0f) { fprintf(stderr, "error in brown_c: std must be nonnegative\n"); return 1; }

    if (N==0u) {}
    else if (std<FLT_EPSILON)
    {
        for (size_t n=2u*N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else
    {
        const float M_2PI = (float)(2.0*M_PI);
        float u1, u2, R, smr=0.0f, smi=0.0f;
        uint32_t r, xorshifted, rot;
        uint64_t state = 0u;
        const uint64_t mul = 6364136223846793005u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
        struct timespec ts;

        //Init random num generator
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in brown_c: timespec_get.\n"); perror("timespec_get"); return 1; }
	    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

        //Generate
        if (std==1.0f)
        {
            for (size_t n=N; n>0u; --n)
            {
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((float)r,-32);
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((float)r,-32);
                R = sqrtf(-2.0f*logf(1.0f-u1));
                //smr += R * cosf(M_2PI*u2); *Y++ = smr;
                //smi += R * sinf(M_2PI*u2); *Y++ = smi;
                //Idea to fix overflow
                smr += R * cosf(M_2PI*u2);
                if (isnan(smr) || !isfinite(smr)) { smr = *(Y-1) - R*cosf(M_2PI*u2); }
                *Y++ = smr;
                smi += R * sinf(M_2PI*u2);
                if (isnan(smi) || !isfinite(smi)) { smi = *(Y-1) - R*sinf(M_2PI*u2); }
                *Y++ = smi;
            }
        }
        else
        {
            for (size_t n=N; n>0u; --n)
            {
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((float)r,-32);
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((float)r,-32);
                R = std * sqrtf(-2.0f*logf(1.0f-u1));
                //smr += R * cosf(M_2PI*u2); *Y++ = smr;
                //smi += R * sinf(M_2PI*u2); *Y++ = smi;
                //Idea to fix overflow
                smr += R * cosf(M_2PI*u2);
                if (isnan(smr) || !isfinite(smr)) { smr = *(Y-1) - R*cosf(M_2PI*u2); }
                *Y++ = smr;
                smi += R * sinf(M_2PI*u2);
                if (isnan(smi) || !isfinite(smi)) { smi = *(Y-1) - R*sinf(M_2PI*u2); }
                *Y++ = smi;
            }
        }

        if (zmn)
        {
            smr = smi = 0.0f;
            for (size_t n=N; n>0u; --n) { smi += *--Y; smr += *--Y; }
            smr /= (float)N; smi /= (float)N;
            for (size_t n=N; n>0u; --n) { *Y++ -= smr; *Y++ -= smi; }
        }
    }

    return 0;
}


int brown_z (double *Y, const size_t N, const double std, const int zmn)
{
    if (std<0.0) { fprintf(stderr, "error in brown_z: std must be nonnegative\n"); return 1; }

    if (N==0u) {}
    else if (std<DBL_EPSILON)
    {
        for (size_t n=2u*N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else
    {
        const double M_2PI = 2.0*M_PI;
        double u1, u2, R, smr=0.0, smi=0.0;
        uint32_t r, xorshifted, rot;
        uint64_t state = 0u;
        const uint64_t mul = 6364136223846793005u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
        struct timespec ts;

        //Init random num generator
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in brown_z: timespec_get.\n"); perror("timespec_get"); return 1; }
	    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

        //Generate
        if (std==1.0)
        {
            for (size_t n=N; n>0u; --n)
            {
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((double)r,-32);
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((double)r,-32);
                R = sqrt(-2.0*log(1.0-u1));
                //smr += R * cos(M_2PI*u2); *Y++ = smr;
                //smi += R * sin(M_2PI*u2); *Y++ = smi;
                //Idea to fix overflow
                smr += R * cos(M_2PI*u2);
                if (isnan(smr) || !isfinite(smr)) { smr = *(Y-1) - R*cos(M_2PI*u2); }
                *Y++ = smr;
                smi += R * sin(M_2PI*u2);
                if (isnan(smi) || !isfinite(smi)) { smi = *(Y-1) - R*sin(M_2PI*u2); }
                *Y++ = smi;
            }
        }
        else
        {
            for (size_t n=N; n>0u; --n)
            {
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((double)r,-32);
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((double)r,-32);
                R = std * sqrt(-2.0*log(1.0-u1));
                //smr += R * cos(M_2PI*u2); *Y++ = smr;
                //smi += R * sin(M_2PI*u2); *Y++ = smi;
                //Idea to fix overflow
                smr += R * cos(M_2PI*u2);
                if (isnan(smr) || !isfinite(smr)) { smr = *(Y-1) - R*cos(M_2PI*u2); }
                *Y++ = smr;
                smi += R * sin(M_2PI*u2);
                if (isnan(smi) || !isfinite(smi)) { smi = *(Y-1) - R*sin(M_2PI*u2); }
                *Y++ = smi;
            }
        }

        if (zmn)
        {
            smr = smi = 0.0;
            for (size_t n=N; n>0u; --n) { smi += *--Y; smr += *--Y; }
            smr /= (double)N; smi /= (double)N;
            for (size_t n=N; n>0u; --n) { *Y++ -= smr; *Y++ -= smi; }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
