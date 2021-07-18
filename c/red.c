//Gets random floats from normal distributions with mean=0 and std,
//takes fft, multiplies by 1/f amplitude, and takes ifft.
//Although this requires an extra fft compared to Octave,
//it requires half as many random numbers, and allows control of std.

//This uses modified code from PCG randoms minimal C library,
//but remains stand-alone (no install of PCG required).
//I use the textbook cos and sin formulae to make Gaussian.

#include <stdio.h>
#include <stdint.h>
#include <float.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int red_s (float *Y, const size_t N, const float std, const char zmn);
int red_d (double *Y, const size_t N, const double std, const char zmn);
int red_c (float *Y, const size_t N, const float std, const char zmn);
int red_z (double *Y, const size_t N, const double std, const char zmn);


int red_s (float *Y, const size_t N, const float std, const char zmn)
{
    if (std<0.0f) { fprintf(stderr, "error in red_s: std must be nonnegative\n"); return 1; }

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
        const uint64_t mul = 6364136223846793005u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
        struct timespec ts;

        //Init random num generator
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in red_s: timespec_get.\n"); perror("timespec_get"); return 1; }
	    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

        //Init FFT
        float *X1, *Y1;
        X1 = (float *)fftwf_malloc(N*sizeof(float));
        Y1 = (float *)fftwf_malloc(2u*N*sizeof(float));
        fftwf_plan plan = fftwf_plan_dft_r2c_1d((int)N,X1,(fftwf_complex *)Y1,FFTW_ESTIMATE);
        if (!plan) { fprintf(stderr,"error in red_s: problem creating fftw plan"); return 1; }

        //Init IFFT
        float *Xi;
        Xi = (float *)fftwf_malloc(2u*N*sizeof(float));
        fftwf_plan iplan = fftwf_plan_dft_1d((int)N,(fftwf_complex *)Y1,(fftwf_complex *)Xi,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in red_s: problem creating fftw plan"); return 1; }

        //Generate white noise
        if (std==1.0f)
        {
            for (size_t n=0u; n<N-1u; n+=2u)
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
                *X1++ = R * cosf(M_2PI*u2);
                *X1++ = R * sinf(M_2PI*u2);
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
                *X1++ = R * cosf(M_2PI*u2);
            }
        }
        else
        {
            for (size_t n=0u; n<N-1u; n+=2u)
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
                *X1++ = R * cosf(M_2PI*u2);
                *X1++ = R * sinf(M_2PI*u2);
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
                *X1++ = R * cosf(M_2PI*u2);
            }
        }

        //FFT
        X1 -= N;
        fftwf_execute(plan);

        //Power law (1/f) characteristic
        //Can use b to set some baseline overall gain
        //Here, this is set so unity gain at fundamental freq
        const float b = 1.0f;
        float a;
        Y1 += 2u;
        for (size_t n=1u; n<N/2u+1u; ++n, ++Y1)
        {
            a = b/(float)(n);
            *Y1 *= a; *++Y1 *= a;
        }
        for (size_t n=N/2u+1u; n<N; ++n, ++Y1)
        {
            a = b/(float)(N-n);
            *Y1 *= a; *++Y1 *= a;
        }
        Y1 -= 2u*N;

        //IFFT
        fftwf_execute(iplan);
        for (size_t n=0u; n<N; ++n, Xi+=2u, ++Y) { *Y = *Xi / (float)(N); }
        Xi -= 2u*N;

        //Finish
        fftwf_destroy_plan(plan); fftwf_destroy_plan(iplan);
        fftwf_free(X1); fftwf_free(Y1); fftwf_free(Xi);
    }

    if (zmn)
    {
        float sm = 0.0f;
        for (size_t n=0u; n<N; ++n) { sm += *--Y; }
        sm /= (float)N;
        for (size_t n=0u; n<N; ++n, ++Y) { *Y -= sm; }
    }
    
    return 0;
}


int red_d (double *Y, const size_t N, const double std, const char zmn)
{
    if (std<0.0) { fprintf(stderr, "error in red_d: std must be nonnegative\n"); return 1; }

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
        const uint64_t mul = 6364136223846793005u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
        struct timespec ts;

        //Init random num generator
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in red_d: timespec_get.\n"); perror("timespec_get"); return 1; }
	    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

        //Init FFT
        double *X1, *Y1;
        X1 = (double *)fftw_malloc(N*sizeof(double));
        Y1 = (double *)fftw_malloc(2u*N*sizeof(double));
        fftw_plan plan = fftw_plan_dft_r2c_1d((int)N,X1,(fftw_complex *)Y1,FFTW_ESTIMATE);
        if (!plan) { fprintf(stderr,"error in red_d: problem creating fftw plan"); return 1; }

        //Init IFFT
        double *Xi;
        Xi = (double *)fftw_malloc(2u*N*sizeof(double));
        fftw_plan iplan = fftw_plan_dft_1d((int)N,(fftw_complex *)Y1,(fftw_complex *)Xi,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in red_d: problem creating fftw plan"); return 1; }

        //Generate white noise
        if (std==1.0)
        {
            for (size_t n=0u; n<N-1u; n+=2u)
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
                *X1++ = R * cos(M_2PI*u2);
                *X1++ = R * sin(M_2PI*u2);
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
                *X1++ = R * cos(M_2PI*u2);
            }
        }
        else
        {
            for (size_t n=0u; n<N-1u; n+=2u)
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
                *X1++ = R * cos(M_2PI*u2);
                *X1++ = R * sin(M_2PI*u2);
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
                *X1++ = R * cos(M_2PI*u2);
            }
        }

        //FFT
        X1 -= N;
        fftw_execute(plan);

        //Power law (1/f) characteristic
        //Can use b to set some baseline overall gain
        //Here, this is set so unity gain at fundamental freq
        const double b = 1.0;
        double a;
        Y1 += 2u;
        for (size_t n=1u; n<N/2u+1u; ++n, ++Y1)
        {
            a = b/(double)(n);
            *Y1 *= a; *++Y1 *= a;
        }
        for (size_t n=N/2u+1u; n<N; ++n, ++Y1)
        {
            a = b/(double)(N-n);
            *Y1 *= a; *++Y1 *= a;
        }
        Y1 -= 2u*N;

        //IFFT
        fftw_execute(iplan);
        for (size_t n=0u; n<N; ++n, Xi+=2u, ++Y) { *Y = *Xi / (double)(N); }
        Xi -= 2u*N;

        //Finish
        fftw_destroy_plan(plan); fftw_destroy_plan(iplan);
        fftw_free(X1); fftw_free(Y1); fftw_free(Xi);
    }

    if (zmn)
    {
        double sm = 0.0;
        for (size_t n=0u; n<N; ++n) { sm += *--Y; }
        sm /= (double)N;
        for (size_t n=0u; n<N; ++n, ++Y) { *Y -= sm; }
    }

    return 0;
}


int red_c (float *Y, const size_t N, const float std, const char zmn)
{
    if (std<0.0f) { fprintf(stderr, "error in red_c: std must be nonnegative\n"); return 1; }

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
        const uint64_t mul = 6364136223846793005u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
        struct timespec ts;

        //Init random num generator
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in red_c: timespec_get.\n"); perror("timespec_get"); return 1; }
	    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

        //Init FFT
        float *X1, *Y1;
        X1 = (float *)fftwf_malloc(2u*N*sizeof(float));
        Y1 = (float *)fftwf_malloc(2u*N*sizeof(float));
        fftwf_plan plan = fftwf_plan_dft_1d((int)N,(fftwf_complex *)X1,(fftwf_complex *)Y1,FFTW_FORWARD,FFTW_ESTIMATE);
        if (!plan) { fprintf(stderr,"error in red_c: problem creating fftw plan"); return 1; }

        //Init IFFT
        float *Xi;
        Xi = (float *)fftwf_malloc(2u*N*sizeof(float));
        fftwf_plan iplan = fftwf_plan_dft_1d((int)N,(fftwf_complex *)Y1,(fftwf_complex *)Xi,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in red_c: problem creating fftw plan"); return 1; }

        //Generate white noise
        if (std==1.0f)
        {
            for (size_t n=0u; n<N; ++n)
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
                *X1++ = R * cosf(M_2PI*u2);
                *X1++ = R * sinf(M_2PI*u2);
            }
        }
        else
        {
            for (size_t n=0u; n<N; ++n)
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
                *X1++ = R * cosf(M_2PI*u2);
                *X1++ = R * sinf(M_2PI*u2);
            }
        }

        //FFT
        X1 -= 2u*N;
        fftwf_execute(plan);

        //Power law (1/f) characteristic
        //Can use b to set some baseline overall gain
        //Here, this is set so unity gain at fundamental freq
        const float b = 1.0f;
        float a;
        Y1 += 2u;
        for (size_t n=1u; n<N/2u+1u; ++n, ++Y1)
        {
            a = b/(float)(n);
            *Y1 *= a; *++Y1 *= a;
        }
        for (size_t n=N/2u+1u; n<N; ++n, ++Y1)
        {
            a = b/(float)(N-n);
            *Y1 *= a; *++Y1 *= a;
        }
        Y1 -= 2u*N;

        //IFFT
        fftwf_execute(iplan);
        for (size_t n=0u; n<2u*N; ++n, ++Xi, ++Y) { *Y = *Xi / (float)(N); }
        Xi -= 2u*N;

        //Finish
        fftwf_destroy_plan(plan); fftwf_destroy_plan(iplan);
        fftwf_free(X1); fftwf_free(Y1); fftwf_free(Xi);
    }

    if (zmn)
    {
        float smr=0.0f, smi=0.0f;
        for (size_t n=0u; n<N; ++n) { smi += *--Y; smr += *--Y; }
        smr /= (float)N; smi /= (float)N;
        for (size_t n=0u; n<N; ++n) { *Y++ -= smr; *Y++ -= smi; }
    }

    return 0;
}


int red_z (double *Y, const size_t N, const double std, const char zmn)
{
    if (std<0.0) { fprintf(stderr, "error in red_z: std must be nonnegative\n"); return 1; }

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
        const uint64_t mul = 6364136223846793005u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
        struct timespec ts;

        //Init random num generator
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in red_z: timespec_get.\n"); perror("timespec_get"); return 1; }
	    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

        //Init FFT
        double *X1, *Y1;
        X1 = (double *)fftw_malloc(2u*N*sizeof(double));
        Y1 = (double *)fftw_malloc(2u*N*sizeof(double));
        fftw_plan plan = fftw_plan_dft_1d((int)N,(fftw_complex *)X1,(fftw_complex *)Y1,FFTW_FORWARD,FFTW_ESTIMATE);
        if (!plan) { fprintf(stderr,"error in red_z: problem creating fftw plan"); return 1; }

        //Init IFFT
        double *Xi;
        Xi = (double *)fftw_malloc(2u*N*sizeof(double));
        fftw_plan iplan = fftw_plan_dft_1d((int)N,(fftw_complex *)Y1,(fftw_complex *)Xi,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in red_z: problem creating fftw plan"); return 1; }

        //Generate white noise
        if (std==1.0)
        {
            for (size_t n=0u; n<N; ++n)
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
                *X1++ = R * cos(M_2PI*u2);
                *X1++ = R * sin(M_2PI*u2);
            }
        }
        else
        {
            for (size_t n=0u; n<N; ++n)
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
                *X1++ = R * cos(M_2PI*u2);
                *X1++ = R * sin(M_2PI*u2);
            }
        }

        //FFT
        X1 -= 2u*N;
        fftw_execute(plan);

        //Power law (1/f) characteristic
        //Can use b to set some baseline overall gain
        //Here, this is set so unity gain at fundamental freq
        const double b = 1.0;
        double a;
        Y1 += 2u;
        for (size_t n=1u; n<N/2u+1u; ++n, ++Y1)
        {
            a = b/(double)(n);
            *Y1 *= a; *++Y1 *= a;
        }
        for (size_t n=N/2u+1u; n<N; ++n, ++Y1)
        {
            a = b/(double)(N-n);
            *Y1 *= a; *++Y1 *= a;
        }
        Y1 -= 2u*N;

        //IFFT
        fftw_execute(iplan);
        for (size_t n=0u; n<2u*N; ++n, ++Xi, ++Y) { *Y = *Xi / (double)(N); }
        Xi -= 2u*N;

        //Finish
        fftw_destroy_plan(plan); fftw_destroy_plan(iplan);
        fftw_free(X1); fftw_free(Y1); fftw_free(Xi);
    }

    if (zmn)
    {
        double smr=0.0, smi = 0.0;
        for (size_t n=0u; n<N; ++n) { smi += *--Y; smr += *--Y; }
        smr /= (double)N; smi /= (double)N;
        for (size_t n=0u; n<N; ++n) { *Y++ -= smr; *Y++ -= smi; }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
