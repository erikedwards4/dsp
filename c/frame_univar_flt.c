//This takes a continuous univariate time series (vector X),
//and outputs a set of frames (matrix Y).

//There uses different conventions than frame_univar.c.

//The present framing function is for flexibility; specifically,
//this allows stp to be a floating-point value (which can easily arise in practice).
//The first center-sample (c0) is also a floating-point value for full flexibility.
//Since W is specified, the user can sample only part of X, or well beyond X, if desired.

//Conceptually, this works like a quick interpolation of X from 0:N-1 to c0:stp:c0+stp*W,
//where interpolated points are set to the nearest-neighbor, and extrapolated points to zero.

#include <stdio.h>
#include <float.h>
#include <math.h>
#include "codee_dsp.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int frame_univar_flt_s (float *Y, const float *X, const size_t N, const size_t L, const size_t W, const float c0, const float stp);
int frame_univar_flt_d (double *Y, const double *X, const size_t N, const size_t L, const size_t W, const double c0, const double stp);
int frame_univar_flt_c (float *Y, const float *X, const size_t N, const size_t L, const size_t W, const float c0, const float stp);
int frame_univar_flt_z (double *Y, const double *X, const size_t N, const size_t L, const size_t W, const double c0, const double stp);


int frame_univar_flt_s (float *Y, const float *X, const size_t N, const size_t L, const size_t W, const float c0, const float stp)
{
    if (W>N) { fprintf(stderr,"error in frame_univar_flt_s: W must be <= N (length X)\n"); return 1; }
    if (L>N) { fprintf(stderr,"error in frame_univar_flt_s: L must be <= N (length X)\n"); return 1; }
    if (L<1u) { fprintf(stderr,"error in frame_univar_flt_s: L must be positive\n"); return 1; }
    if (N<1u) { fprintf(stderr,"error in frame_univar_flt_s: N (length X) must be positive\n"); return 1; }
    if (c0>(float)(N-1u)) { fprintf(stderr,"error in frame_univar_flt_s: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (stp<FLT_EPSILON) { fprintf(stderr,"error in frame_univar_flt_s: stp (step size) must be positive\n"); return 1; }

    if (N==0u || W==0u) {}
    else
    {
        const size_t Lpre = L/2u;                   //nsamps before center samp
        const size_t Lpost = L-Lpre-1u;             //nsamps after center samp
        size_t w = 0u;                              //current frame
        float cc = c0;                              //current exact center-samp
        int cs = (int)roundf(c0);                   //current rounded center-samp
        int ss, es = cs + (int)Lpost;               //current rounded start-samp, end-samp
        int prev_cs;                                //previous rounded center-samp

        while (es<0 && w<W)
        {
            for (size_t l=L; l>0u; --l, ++Y) { *Y = 0.0f; }
            ++w; cc += stp; cs = (int)roundf(cc);
            es = cs + (int)Lpost;
        }
        ss = cs - (int)Lpre;
        prev_cs = cs;
        while (ss<0 && w<W)
        {
            for (size_t l=(size_t)(-ss); l>0u; --l, ++Y) { *Y = 0.0f; }
            for (size_t l=(size_t)(-ss); l<L; ++l, ++X, ++Y) { *Y = *X; }
            ++w; cc += stp; cs = (int)roundf(cc);
            X -= (int)L + ss;
            ss = cs - (int)Lpre; prev_cs = cs;
        }
        X += ss;
        es = cs + (int)Lpost;
        while (es<(int)N && w<W)
        {
            for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
            ++w; cc += stp; cs = (int)roundf(cc);
            X += cs - prev_cs - (int)L;
            es = cs + (int)Lpost; prev_cs = cs;
        }
        ss = cs - (int)Lpre;
        while (ss<(int)N && w<W)
        {
            for (size_t l=0u; l<N-(size_t)ss; ++l, ++X, ++Y) { *Y = *X; }
            for (size_t l=N-(size_t)ss; l<L; ++l, ++Y) { *Y = 0.0f; }
            ++w; cc += stp; cs = (int)roundf(cc);
            X += cs - prev_cs - (int)N + ss;
            ss = cs - (int)Lpre; prev_cs = cs;
        }
        while (w<W)
        {
            for (size_t l=L; l>0u; --l, ++Y) { *Y = 0.0f; }
            ++w;
        }
    }

    return 0;
}


int frame_univar_flt_d (double *Y, const double *X, const size_t N, const size_t L, const size_t W, const double c0, const double stp)
{
    if (W>N) { fprintf(stderr,"error in frame_univar_flt_d: W must be <= N (length X)\n"); return 1; }
    if (L>N) { fprintf(stderr,"error in frame_univar_flt_d: L must be <= N (length X)\n"); return 1; }
    if (L<1u) { fprintf(stderr,"error in frame_univar_flt_d: L must be positive\n"); return 1; }
    if (N<1u) { fprintf(stderr,"error in frame_univar_flt_d: N (length X) must be positive\n"); return 1; }
    if (c0>(double)(N-1u)) { fprintf(stderr,"error in frame_univar_flt_d: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (stp<DBL_EPSILON) { fprintf(stderr,"error in frame_univar_flt_d: stp (step size) must be positive\n"); return 1; }

    if (N==0u || W==0u) {}
    else
    {
        const size_t Lpre = L/2u;                   //nsamps before center samp
        const size_t Lpost = L-Lpre-1u;             //nsamps after center samp
        size_t w = 0u;                              //current frame
        double cc = c0;                             //current exact center-samp
        int cs = (int)round(c0);                    //current rounded center-samp
        int ss, es = cs + (int)Lpost;               //current rounded start-samp, end-samp
        int prev_cs;                                //previous rounded center-samp

        while (es<0 && w<W)
        {
            for (size_t l=L; l>0u; --l, ++Y) { *Y = 0.0; }
            ++w; cc += stp; cs = (int)round(cc);
            es = cs + (int)Lpost;
        }
        ss = cs - (int)Lpre;
        prev_cs = cs;
        while (ss<0 && w<W)
        {
            for (size_t l=(size_t)(-ss); l>0u; --l, ++Y) { *Y = 0.0; }
            for (size_t l=(size_t)(-ss); l<L; ++l, ++X, ++Y) { *Y = *X; }
            ++w; cc += stp; cs = (int)round(cc);
            X -= (int)L + ss;
            ss = cs - (int)Lpre; prev_cs = cs;
        }
        X += ss;
        es = cs + (int)Lpost;
        while (es<(int)N && w<W)
        {
            for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
            ++w; cc += stp; cs = (int)round(cc);
            X += cs - prev_cs - (int)L;
            es = cs + (int)Lpost; prev_cs = cs;
        }
        ss = cs - (int)Lpre;
        while (ss<(int)N && w<W)
        {
            for (size_t l=0u; l<N-(size_t)ss; ++l, ++X, ++Y) { *Y = *X; }
            for (size_t l=N-(size_t)ss; l<L; ++l, ++Y) { *Y = 0.0; }
            ++w; cc += stp; cs = (int)round(cc);
            X += cs - prev_cs - (int)N + ss;
            ss = cs - (int)Lpre; prev_cs = cs;
        }
        while (w<W)
        {
            for (size_t l=L; l>0u; --l, ++Y) { *Y = 0.0; }
            ++w;
        }
    }

    return 0;
}


int frame_univar_flt_c (float *Y, const float *X, const size_t N, const size_t L, const size_t W, const float c0, const float stp)
{
    if (W>N) { fprintf(stderr,"error in frame_univar_flt_c: W must be <= N (length X)\n"); return 1; }
    if (L>N) { fprintf(stderr,"error in frame_univar_flt_c: L must be <= N (length X)\n"); return 1; }
    if (L<1u) { fprintf(stderr,"error in frame_univar_flt_c: L must be positive\n"); return 1; }
    if (N<1u) { fprintf(stderr,"error in frame_univar_flt_c: N (length X) must be positive\n"); return 1; }
    if (c0>(float)(N-1u)) { fprintf(stderr,"error in frame_univar_flt_c: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (stp<FLT_EPSILON) { fprintf(stderr,"error in frame_univar_flt_c: stp (step size) must be positive\n"); return 1; }

    if (N==0u || W==0u) {}
    else
    {
        const size_t Lpre = L/2u;                   //nsamps before center samp
        const size_t Lpost = L-Lpre-1u;             //nsamps after center samp
        size_t w = 0u;                              //current frame
        float cc = c0;                              //current exact center-samp
        int cs = (int)roundf(c0);                   //current rounded center-samp
        int ss, es = cs + (int)Lpost;               //current rounded start-samp, end-samp
        int prev_cs;                                //previous rounded center-samp

        while (es<0 && w<W)
        {
            for (size_t l=L; l>0u; --l, ++Y) { *Y = 0.0f; *++Y = 0.0f; }
            ++w; cc += stp; cs = (int)roundf(cc);
            es = cs + (int)Lpost;
        }
        ss = cs - (int)Lpre;
        prev_cs = cs;
        while (ss<0 && w<W)
        {
            for (size_t l=(size_t)(-ss); l>0u; --l, ++Y) { *Y = 0.0f; *++Y = 0.0f; }
            for (size_t l=(size_t)(-ss); l<L; ++l, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            ++w; cc += stp; cs = (int)roundf(cc);
            X -= 2*((int)L+ss);
            ss = cs - (int)Lpre; prev_cs = cs;
        }
        X += 2*ss;
        es = cs + (int)Lpost;
        while (es<(int)N && w<W)
        {
            for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            ++w; cc += stp; cs = (int)roundf(cc);
            X += 2*(cs-prev_cs-(int)L);
            es = cs + (int)Lpost; prev_cs = cs;
        }
        ss = cs - (int)Lpre;
        while (ss<(int)N && w<W)
        {
            for (size_t l=0u; l<N-(size_t)ss; ++l, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            for (size_t l=N-(size_t)ss; l<L; ++l, ++Y) { *Y = 0.0f; *++Y = 0.0f; }
            ++w; cc += stp; cs = (int)roundf(cc);
            X += 2*(cs-prev_cs-(int)N+ss);
            ss = cs - (int)Lpre; prev_cs = cs;
        }
        while (w<W)
        {
            for (size_t l=L; l>0u; --l, ++Y) { *Y = 0.0f; *++Y = 0.0f; }
            ++w;
        }
    }

    return 0;
}


int frame_univar_flt_z (double *Y, const double *X, const size_t N, const size_t L, const size_t W, const double c0, const double stp)
{
    if (W>N) { fprintf(stderr,"error in frame_univar_flt_z: W must be <= N (length X)\n"); return 1; }
    if (L>N) { fprintf(stderr,"error in frame_univar_flt_z: L must be <= N (length X)\n"); return 1; }
    if (L<1u) { fprintf(stderr,"error in frame_univar_flt_z: L must be positive\n"); return 1; }
    if (N<1u) { fprintf(stderr,"error in frame_univar_flt_z: N (length X) must be positive\n"); return 1; }
    if (c0>(double)(N-1u)) { fprintf(stderr,"error in frame_univar_flt_z: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (stp<DBL_EPSILON) { fprintf(stderr,"error in frame_univar_flt_z: stp (step size) must be positive\n"); return 1; }

    if (N==0u || W==0u) {}
    else
    {
        const size_t Lpre = L/2u;                   //nsamps before center samp
        const size_t Lpost = L-Lpre-1u;             //nsamps after center samp
        size_t w = 0u;                              //current frame
        double cc = c0;                             //current exact center-samp
        int cs = (int)round(c0);                    //current rounded center-samp
        int ss, es = cs + (int)Lpost;               //current rounded start-samp, end-samp
        int prev_cs;                                //previous rounded center-samp

        while (es<0 && w<W)
        {
            for (size_t l=L; l>0u; --l, ++Y) { *Y = 0.0; *++Y = 0.0; }
            ++w; cc += stp; cs = (int)round(cc);
            es = cs + (int)Lpost;
        }
        ss = cs - (int)Lpre;
        prev_cs = cs;
        while (ss<0 && w<W)
        {
            for (size_t l=(size_t)(-ss); l>0u; --l, ++Y) { *Y = 0.0; *++Y = 0.0; }
            for (size_t l=(size_t)(-ss); l<L; ++l, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            ++w; cc += stp; cs = (int)round(cc);
            X -= 2*((int)L+ss);
            ss = cs - (int)Lpre; prev_cs = cs;
        }
        X += 2*ss;
        es = cs + (int)Lpost;
        while (es<(int)N && w<W)
        {
            for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            ++w; cc += stp; cs = (int)round(cc);
            X += 2*(cs-prev_cs-(int)L);
            es = cs + (int)Lpost; prev_cs = cs;
        }
        ss = cs - (int)Lpre;
        while (ss<(int)N && w<W)
        {
            for (size_t l=0u; l<N-(size_t)ss; ++l, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            for (size_t l=N-(size_t)ss; l<L; ++l, ++Y) { *Y = 0.0; *++Y = 0.0; }
            ++w; cc += stp; cs = (int)round(cc);
            X += 2*(cs-prev_cs-(int)N+ss);
            ss = cs - (int)Lpre; prev_cs = cs;
        }
        while (w<W)
        {
            for (size_t l=L; l>0u; --l, ++Y) { *Y = 0.0; *++Y = 0.0; }
            ++w;
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
