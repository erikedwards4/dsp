//This takes a continuous univariate time series (vector X),
//and outputs a set of frames (matrix Y).

//There is an established convention from HTK through Kaldi to Librosa,
//involving a setting "snip-edges", which usually has default "true".
//For those conventions and integer-only stp (frame-shift), see frame_univar.c.

//The present framing function is for flexibility; specifically,
//this allows stp to be a floating-point value (which can easily arise in practice).
//The first center-sample (c0) is also a floating-point value for full flexibility.
//Since W is still specified, the user can sample only part of X, or well beyond X, if desired.

//Conceptually, this is identical to a quick interpolation of X from 0:N-1 to c0:stp:c0+stp*W,
//where interpolated points are set to the nearest-neighbor, and extrapolated points to zero.
//(Thus, this could easily be changed to a linear interpolation.)

#include <stdio.h>
#include <float.h>
#include <math.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int frame_univar_float_s (float *Y, const float *X, const size_t N, const size_t W, const size_t L, const float c0, const float stp);
int frame_univar_float_d (double *Y, const double *X, const size_t N, const size_t W, const size_t L, const double c0, const double stp);
int frame_univar_float_c (float *Y, const float *X, const size_t N, const size_t W, const size_t L, const float c0, const float stp);
int frame_univar_float_z (double *Y, const double *X, const size_t N, const size_t W, const size_t L, const double c0, const double stp);


int frame_univar_float_s (float *Y, const float *X, const size_t N, const size_t W, const size_t L, const float c0, const float stp)
{
    if (L>=N) { fprintf(stderr,"error in frame_univar_float_s: L must be < N (length X)\n"); return 1; }
    if (L<1u) { fprintf(stderr,"error in frame_univar_float_s: L must be positive\n"); return 1; }
    if (N<1u) { fprintf(stderr,"error in frame_univar_float_s: N (length X) must be positive\n"); return 1; }
    if (c0>N-1u) { fprintf(stderr,"error in frame_univar_float_s: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (stp<FLT_EPSILON) { fprintf(stderr,"error in frame_univar_float_s: stp (step size) must be positive\n"); return 1; }

    if (N==0u || W==0u) {}
    else
    {
        const size_t Lpre = L/2u;                   //nsamps before center samp
        const size_t Lpost = L-L/2u-1u;             //nsamps after center samp
        size_t w = 0u;                              //current frame
        float ss, cs = c0, es=cs+(float)Lpost;      //current start-samp, center-samp, end-samp
        float cc = c0;                              //current exact center-samp
        int cs = (int)roundf(cc);                   //current rounded center-samp
        int ss, es = cs + (int)Lpost;               //current rounded start-samp, end-samp
        int prev_cs;                                //previous rounded center-samp

        while (es<0 && w<W)
        {
            for (size_t l=0u; l<L; ++l, ++Y) { *Y = 0.0f; }
            ++w; cc += stp; cs = (int)roundf(cc);
            es = cs + (int)Lpost;
        }
        ss = cs - (int)Lpre;
        prev_cs = cs;
        while (ss<0 && w<W)
        {
            for (size_t l=0u; l<(size_t)(-ss); ++l, ++Y) { *Y = 0.0f; }
            for (size_t l=(size_t)(-ss); l<L; ++l, ++X, ++Y) { *Y = *X; }
            X -= L - (size_t)(-ss);
            ++w; cc += stp; cs = (int)roundf(cc);
            X += cs - prev_cs; prev_cs = cs;
        }
        es = cs + (int)Lpost;
        while (es<(int)N && w<W)
        {
            for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = *X; }
            X -= L;
            es = cs + (int)Lpost;
            ++w; cc += stp; cs = (int)roundf(cc);
            X += cs - prev_cs; prev_cs = cs;
        }
        ss = cs - (int)Lpre;
        while (ss<(int)N && w<W)
        {
            for (size_t l=0u; l<N+Lpre-(size_t)cs; ++l, ++X, ++Y) { *Y = *X; }
            for (size_t l=N+Lpre-(size_t)cs; l<L; ++l, ++Y) { *Y = 0.0f; }
            X -= N + Lpre - (size_t)cs;
            ss = cs - (int)Lpre;
            ++w; cc += stp; cs = (int)roundf(cc);
            X += cs - prev_cs; prev_cs = cs;
        }
        while (w<W)
        {
            for (size_t l=0u; l<L; ++l, ++Y) { *Y = 0.0f; }
            ++w;
        }
    }
}



#ifdef __cplusplus
}
}
#endif
