//This gets ZCs as usual from X1, and then applies window X2 to get rate.
//
//No internal normalization by the sum of X2 is done here:
//X2 should sum to 1 if the output Y is interpreted as a moving average.

//I profiled a few ways to do the dot product between bool *Z and float *X2.
//1. for (size_t l=0u; l<Lw; ++l, ++Z, ++X2) { sm += (float)*Z * *X2; }
//2. for (size_t l=0u; l<Lw; ++l, ++Z, ++X2) { if (*Z) { sm += *X2; } }
//3. a while loop method.
//1. was ~10% faster than 2., which was ~10% faster than 3.

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int zcr_windowed_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const size_t dim, const float c0, const float stp, const int going);
int zcr_windowed_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const size_t dim, const double c0, const double stp, const int going);
int zcr_windowed_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const size_t dim, const float c0, const float stp, const int going);
int zcr_windowed_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const size_t dim, const double c0, const double stp, const int going);


int zcr_windowed_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const size_t dim, const float c0, const float stp, const int going)
{
    if (dim>3u) { fprintf(stderr,"error in zcr_windowed_s: dim must be in [0 3]\n"); return 1; }
    if (Lw<1u) { fprintf(stderr,"error in zcr_windowed_s: Lw (winlength) must be positive\n"); return 1; }
    if (stp<FLT_EPSILON) { fprintf(stderr,"error in zcr_windowed_s: stp (step size) must be positive\n"); return 1; }
    if (going!=0 && going!=1 && going!=-1) { fprintf(stderr,"error in zcr_windowed_s: going must be in {-1,0,1}\n"); return 1; }

    const size_t N = R*C*S*H;
    if (N==0u) {fprintf(stderr,"error in zcr_windowed_s: X1 is empty (size of at least one dim is 0)\n"); return 1; }

    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lx<1u) { fprintf(stderr,"error in zcr_windowed_s: Lx (length of vecs in X1) must be positive\n"); return 1; }
    if (Lw>Lx) { fprintf(stderr,"error in zcr_windowed_s: Lw (winlength) must be <= Lx (length of vecs in X1)\n"); return 1; }
    if (W>Lx) { fprintf(stderr,"error in zcr_windowed_s: W must be <= Lx (length of vecs in X1)\n"); return 1; }
    if (c0>(float)(Lx-1u)) { fprintf(stderr,"error in zcr_windowed_s: c0 (center samp of 1st frame) must be < Lx (length of vecs in X1)\n"); return 1; }

    if (W==0u) {}
    else
    {
        const size_t Lpre = Lw/2u;                  //nsamps before center samp
        const size_t Lpost = Lw-Lpre-1u;            //nsamps after center samp
        size_t w = 0u;                              //current frame
        float cc = c0;                              //current exact center-samp
        int cs = (int)roundf(c0);                   //current rounded center-samp
        int ss, es = cs + (int)Lpost;               //current rounded start-samp, end-samp
        int prev_cs;                                //previous rounded center-samp
        int *Z, s, sp;                              //ints for ZCR
        float sm;                                   //local sum for ZCR

        //Allocate
        if (!(Z=(int *)malloc(Lx*sizeof(int)))) { fprintf(stderr,"error in zcr_windowed_s: problem with malloc. "); perror("malloc"); return 1; }

        //Get block and vec sizes
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            //For each vec in X1
            for (size_t v=0u; v<V; ++v)
            {
                //Get Z
                if (going==0)
                {
                    sp = (*X1++<0.0f); *Z++ = 0;
                    for (size_t l=1u; l<Lx; ++l, ++X1, ++Z) { s = (*X1<0.0f); *Z = (s!=sp); sp = s; }
                }
                else if (going==1)
                {
                    sp = (*X1++>=0.0f); *Z++ = 0;
                    for (size_t l=1u; l<Lx; ++l, ++X1, ++Z) { s = (*X1>=0.0f); *Z = s*(s!=sp); sp = s; }
                }
                else if (going==-1)
                {
                    sp = (*X1++<0.0f); *Z++ = 0;
                    for (size_t l=1u; l<Lx; ++l, ++X1, ++Z) { s = (*X1<0.0f); *Z = s*(s!=sp); sp = s; }
                }
                Z -= Lx;

                //Initialize windowing
                w = 0u; cc = c0; cs = (int)round(c0); es = cs + (int)Lpost;

                //Windows before first samp
                while (es<0 && w<W)
                {
                    *Y++ = 0.0f;
                    ++w; cc += stp; cs = (int)roundf(cc);
                    es = cs + (int)Lpost;
                }
                ss = cs - (int)Lpre;
                prev_cs = cs;

                //Windows overlapping first samp
                while (ss<0 && w<W)
                {
                    sm = 0.0f; X2 -= ss;
                    for (size_t l=(size_t)(-ss); l<Lw; ++l, ++X2, ++Z) { sm += (float)(*Z) * *X2; }
                    *Y++ = sm;
                    ++w; cc += stp; cs = (int)roundf(cc);
                    Z -= (int)Lw + ss; X2 -= Lw;
                    ss = cs - (int)Lpre; prev_cs = cs;
                }
                Z += ss;
                es = cs + (int)Lpost;
                
                //Windows fully within sig
                while (es<(int)Lx && w<W)
                {
                    sm = 0.0f;
                    for (size_t l=0u; l<Lw; ++l, ++Z, ++X2) { sm += (float)*Z * *X2; }
                    *Y++ = sm;
                    ++w; cc += stp; cs = (int)roundf(cc);
                    Z += cs - prev_cs - (int)Lw; X2 -= Lw;
                    es = cs + (int)Lpost; prev_cs = cs;
                }
                ss = cs - (int)Lpre;

                //Windows overlapping last samp
                while (ss<(int)Lx && w<W)
                {
                    sm = 0.0f;
                    for (size_t l=0u; l<Lx-(size_t)ss; ++l, ++Z, ++X2) { sm += (float)(*Z) * *X2; }
                    *Y++ = sm;
                    X2 += ss - (int)Lx;
                    ++w; cc += stp; cs = (int)roundf(cc);
                    Z += cs - prev_cs - (int)Lx + ss;
                    ss = cs - (int)Lpre; prev_cs = cs;
                }
                Z -= cs - prev_cs + ss;

                //Windows after last samp
                while (w<W) { *Y++ = 0.0f; ++w; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X1+=B*(Lx-1u), Y+=B*(W-1u))
            {
                for (size_t b=0u; b<B; ++b, X1-=K*Lx-1u, Y-=K*W-1u)
                {
                    //Get Z
                    if (going==0)
                    {
                        sp = (*X1<0.0f); X1+=K; *Z++ = 0;
                        for (size_t l=1u; l<Lx; ++l, X1+=K, ++Z) { s = (*X1<0.0f); *Z = (s!=sp); sp = s; }
                    }
                    else if (going==1)
                    {
                        sp = (*X1>=0.0f); X1+=K; *Z++ = 0;
                        for (size_t l=1u; l<Lx; ++l, X1+=K, ++Z) { s = (*X1>=0.0f); *Z = s*(s!=sp); sp = s; }
                    }
                    else if (going==-1)
                    {
                        sp = (*X1<0.0f); X1+=K; *Z++ = 0;
                        for (size_t l=1u; l<Lx; ++l, X1+=K, ++Z) { s = (*X1<0.0f); *Z = s*(s!=sp); sp = s; }
                    }
                    Z -= Lx;

                    //Initialize windowing
                    w = 0u; cc = c0; cs = (int)round(c0); es = cs + (int)Lpost;

                    //Windows before first samp
                    while (es<0 && w<W)
                    {
                        *Y = 0.0f; Y += K;
                        ++w; cc += stp; cs = (int)roundf(cc);
                        es = cs + (int)Lpost;
                    }
                    ss = cs - (int)Lpre;
                    prev_cs = cs;

                    //Windows overlapping first samp
                    while (ss<0 && w<W)
                    {
                        sm = 0.0f; X2 -= ss;
                        for (size_t l=(size_t)(-ss); l<Lw; ++l, ++X2, ++Z) { sm += (float)(*Z) * *X2; }
                        *Y = sm; Y += K;
                        ++w; cc += stp; cs = (int)roundf(cc);
                        Z -= (int)Lw + ss; X2 -= Lw;
                        ss = cs - (int)Lpre; prev_cs = cs;
                    }
                    Z += ss;
                    es = cs + (int)Lpost;
                    
                    //Windows fully within sig
                    while (es<(int)Lx && w<W)
                    {
                        sm = 0.0f;
                        for (size_t l=0u; l<Lw; ++l, ++Z, ++X2) { sm += (float)*Z * *X2; }
                        *Y = sm; Y += K;
                        ++w; cc += stp; cs = (int)roundf(cc);
                        Z += cs - prev_cs - (int)Lw; X2 -= Lw;
                        es = cs + (int)Lpost; prev_cs = cs;
                    }
                    ss = cs - (int)Lpre;

                    //Windows overlapping last samp
                    while (ss<(int)Lx && w<W)
                    {
                        sm = 0.0f;
                        for (size_t l=0u; l<Lx-(size_t)ss; ++l, ++Z, ++X2) { sm += (float)(*Z) * *X2; }
                        *Y = sm; Y += K;
                        X2 += ss - (int)Lx;
                        ++w; cc += stp; cs = (int)roundf(cc);
                        Z += cs - prev_cs - (int)Lx + ss;
                        ss = cs - (int)Lpre; prev_cs = cs;
                    }
                    Z -= cs - prev_cs + ss;

                    //Windows after last samp
                    while (w<W) { *Y = 0.0f; Y += K; ++w; }
                }
            }
        }
        
        free(Z);
    }

    return 0;
}


int zcr_windowed_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const size_t dim, const double c0, const double stp, const int going)
{
    if (dim>3u) { fprintf(stderr,"error in zcr_windowed_d: dim must be in [0 3]\n"); return 1; }
    if (Lw<1u) { fprintf(stderr,"error in zcr_windowed_d: Lw (winlength) must be positive\n"); return 1; }
    if (stp<DBL_EPSILON) { fprintf(stderr,"error in zcr_windowed_d: stp (step size) must be positive\n"); return 1; }
    if (going!=0 && going!=1 && going!=-1) { fprintf(stderr,"error in zcr_windowed_d: going must be in {-1,0,1}\n"); return 1; }

    const size_t N = R*C*S*H;
    if (N==0u) {fprintf(stderr,"error in zcr_windowed_d: X1 is empty (size of at least one dim is 0)\n"); return 1; }

    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lx<1u) { fprintf(stderr,"error in zcr_windowed_d: Lx (length of vecs in X1) must be positive\n"); return 1; }
    if (Lw>Lx) { fprintf(stderr,"error in zcr_windowed_d: Lw (winlength) must be <= Lx (length of vecs in X1)\n"); return 1; }
    if (W>Lx) { fprintf(stderr,"error in zcr_windowed_d: W must be <= Lx (length of vecs in X1)\n"); return 1; }
    if (c0>(double)(Lx-1u)) { fprintf(stderr,"error in zcr_windowed_d: c0 (center samp of 1st frame) must be < Lx (length of vecs in X1)\n"); return 1; }

    if (W==0u) {}
    else
    {
        const size_t Lpre = Lw/2u;                  //nsamps before center samp
        const size_t Lpost = Lw-Lpre-1u;            //nsamps after center samp
        size_t w = 0u;                              //current frame
        double cc = c0;                             //current exact center-samp
        int cs = (int)round(c0);                    //current rounded center-samp
        int ss, es = cs + (int)Lpost;               //current rounded start-samp, end-samp
        int prev_cs;                                //previous rounded center-samp
        int *Z, s, sp;                              //ints for ZCR
        double sm;                                  //local sum for ZCR

        //Allocate
        if (!(Z=(int *)malloc(Lx*sizeof(int)))) { fprintf(stderr,"error in zcr_windowed_d: problem with malloc. "); perror("malloc"); return 1; }

        //Get block and vec sizes
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            //For each vec in X1
            for (size_t v=0u; v<V; ++v)
            {
                //Get Z
                if (going==0)
                {
                    sp = (*X1++<0.0); *Z++ = 0;
                    for (size_t l=1u; l<Lx; ++l, ++X1, ++Z) { s = (*X1<0.0); *Z = (s!=sp); sp = s; }
                }
                else if (going==1)
                {
                    sp = (*X1++>=0.0); *Z++ = 0;
                    for (size_t l=1u; l<Lx; ++l, ++X1, ++Z) { s = (*X1>=0.0); *Z = s*(s!=sp); sp = s; }
                }
                else if (going==-1)
                {
                    sp = (*X1++<0.0); *Z++ = 0;
                    for (size_t l=1u; l<Lx; ++l, ++X1, ++Z) { s = (*X1<0.0); *Z = s*(s!=sp); sp = s; }
                }
                Z -= Lx;

                //Initialize windowing
                w = 0u; cc = c0; cs = (int)round(c0); es = cs + (int)Lpost;

                //Windows before first samp
                while (es<0 && w<W)
                {
                    *Y++ = 0.0;
                    ++w; cc += stp; cs = (int)round(cc);
                    es = cs + (int)Lpost;
                }
                ss = cs - (int)Lpre;
                prev_cs = cs;

                //Windows overlapping first samp
                while (ss<0 && w<W)
                {
                    sm = 0.0; X2 -= ss;
                    for (size_t l=(size_t)(-ss); l<Lw; ++l, ++X2, ++Z) { sm += (double)(*Z) * *X2; }
                    *Y++ = sm;
                    ++w; cc += stp; cs = (int)round(cc);
                    Z -= (int)Lw + ss; X2 -= Lw;
                    ss = cs - (int)Lpre; prev_cs = cs;
                }
                Z += ss;
                es = cs + (int)Lpost;
                
                //Windows fully within sig
                while (es<(int)Lx && w<W)
                {
                    sm = 0.0;
                    for (size_t l=0u; l<Lw; ++l, ++Z, ++X2) { sm += (double)*Z * *X2; }
                    *Y++ = sm;
                    ++w; cc += stp; cs = (int)round(cc);
                    Z += cs - prev_cs - (int)Lw; X2 -= Lw;
                    es = cs + (int)Lpost; prev_cs = cs;
                }
                ss = cs - (int)Lpre;

                //Windows overlapping last samp
                while (ss<(int)Lx && w<W)
                {
                    sm = 0.0;
                    for (size_t l=0u; l<Lx-(size_t)ss; ++l, ++Z, ++X2) { sm += (double)(*Z) * *X2; }
                    *Y++ = sm;
                    X2 += ss - (int)Lx;
                    ++w; cc += stp; cs = (int)round(cc);
                    Z += cs - prev_cs - (int)Lx + ss;
                    ss = cs - (int)Lpre; prev_cs = cs;
                }
                Z -= cs - prev_cs + ss;

                //Windows after last samp
                while (w<W) { *Y++ = 0.0; ++w; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X1+=B*(Lx-1u), Y+=B*(W-1u))
            {
                for (size_t b=0u; b<B; ++b, X1-=K*Lx-1u, Y-=K*W-1u)
                {
                    //Get Z
                    if (going==0)
                    {
                        sp = (*X1<0.0); X1+=K; *Z++ = 0;
                        for (size_t l=1u; l<Lx; ++l, X1+=K, ++Z) { s = (*X1<0.0); *Z = (s!=sp); sp = s; }
                    }
                    else if (going==1)
                    {
                        sp = (*X1>=0.0); X1+=K; *Z++ = 0;
                        for (size_t l=1u; l<Lx; ++l, X1+=K, ++Z) { s = (*X1>=0.0); *Z = s*(s!=sp); sp = s; }
                    }
                    else if (going==-1)
                    {
                        sp = (*X1<0.0); X1+=K; *Z++ = 0;
                        for (size_t l=1u; l<Lx; ++l, X1+=K, ++Z) { s = (*X1<0.0); *Z = s*(s!=sp); sp = s; }
                    }
                    Z -= Lx;

                    //Initialize windowing
                    w = 0u; cc = c0; cs = (int)round(c0); es = cs + (int)Lpost;

                    //Windows before first samp
                    while (es<0 && w<W)
                    {
                        *Y = 0.0; Y += K;
                        ++w; cc += stp; cs = (int)round(cc);
                        es = cs + (int)Lpost;
                    }
                    ss = cs - (int)Lpre;
                    prev_cs = cs;

                    //Windows overlapping first samp
                    while (ss<0 && w<W)
                    {
                        sm = 0.0; X2 -= ss;
                        for (size_t l=(size_t)(-ss); l<Lw; ++l, ++X2, ++Z) { sm += (double)(*Z) * *X2; }
                        *Y = sm; Y += K;
                        ++w; cc += stp; cs = (int)round(cc);
                        Z -= (int)Lw + ss; X2 -= Lw;
                        ss = cs - (int)Lpre; prev_cs = cs;
                    }
                    Z += ss;
                    es = cs + (int)Lpost;
                    
                    //Windows fully within sig
                    while (es<(int)Lx && w<W)
                    {
                        sm = 0.0;
                        for (size_t l=0u; l<Lw; ++l, ++Z, ++X2) { sm += (double)*Z * *X2; }
                        *Y = sm; Y += K;
                        ++w; cc += stp; cs = (int)round(cc);
                        Z += cs - prev_cs - (int)Lw; X2 -= Lw;
                        es = cs + (int)Lpost; prev_cs = cs;
                    }
                    ss = cs - (int)Lpre;

                    //Windows overlapping last samp
                    while (ss<(int)Lx && w<W)
                    {
                        sm = 0.0;
                        for (size_t l=0u; l<Lx-(size_t)ss; ++l, ++Z, ++X2) { sm += (double)(*Z) * *X2; }
                        *Y = sm; Y += K;
                        X2 += ss - (int)Lx;
                        ++w; cc += stp; cs = (int)round(cc);
                        Z += cs - prev_cs - (int)Lx + ss;
                        ss = cs - (int)Lpre; prev_cs = cs;
                    }
                    Z -= cs - prev_cs + ss;

                    //Windows after last samp
                    while (w<W) { *Y = 0.0; Y += K; ++w; }
                }
            }
        }
        
        free(Z);
    }

    return 0;
}


int zcr_windowed_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const size_t dim, const float c0, const float stp, const int going)
{
    if (dim>3u) { fprintf(stderr,"error in zcr_windowed_c: dim must be in [0 3]\n"); return 1; }
    if (Lw<1u) { fprintf(stderr,"error in zcr_windowed_c: Lw (winlength) must be positive\n"); return 1; }
    if (stp<FLT_EPSILON) { fprintf(stderr,"error in zcr_windowed_c: stp (step size) must be positive\n"); return 1; }
    if (going!=0 && going!=1 && going!=-1) { fprintf(stderr,"error in zcr_windowed_c: going must be in {-1,0,1}\n"); return 1; }

    const size_t N = R*C*S*H;
    if (N==0u) {fprintf(stderr,"error in zcr_windowed_c: X1 is empty (size of at least one dim is 0)\n"); return 1; }

    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lx<1u) { fprintf(stderr,"error in zcr_windowed_c: Lx (length of vecs in X1) must be positive\n"); return 1; }
    if (Lw>Lx) { fprintf(stderr,"error in zcr_windowed_c: Lw (winlength) must be <= Lx (length of vecs in X1)\n"); return 1; }
    if (W>Lx) { fprintf(stderr,"error in zcr_windowed_c: W must be <= Lx (length of vecs in X1)\n"); return 1; }
    if (c0>(float)(Lx-1u)) { fprintf(stderr,"error in zcr_windowed_c: c0 (center samp of 1st frame) must be < Lx (length of vecs in X1)\n"); return 1; }

    if (W==0u) {}
    else
    {
        const size_t Lpre = Lw/2u;                  //nsamps before center samp
        const size_t Lpost = Lw-Lpre-1u;            //nsamps after center samp
        size_t w = 0u;                              //current frame
        float cc = c0;                              //current exact center-samp
        int cs = (int)roundf(c0);                   //current rounded center-samp
        int ss, es = cs + (int)Lpost;               //current rounded start-samp, end-samp
        int prev_cs;                                //previous rounded center-samp
        int *Z, s, sp;                              //ints for ZCR
        float sm;                                   //local sum for ZCR

        //Allocate
        if (!(Z=(int *)malloc(Lx*sizeof(int)))) { fprintf(stderr,"error in zcr_windowed_c: problem with malloc. "); perror("malloc"); return 1; }

        //Get block and vec sizes
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        //Comment this out to use real part of X1 for ZCs
        ++X1;

        for (size_t g=0u; g<G; ++g, X1+=2u*B*(Lx-1u), Y+=B*(W-1u))
        {
            for (size_t b=0u; b<B; ++b, X1-=2u*K*Lx-2u, Y-=K*W-1u)
            {
                //Get Z
                if (going==0)
                {
                    sp = (*X1<0.0f); X1+=2u*K; *Z++ = 0;
                    for (size_t l=1u; l<Lx; ++l, X1+=2u*K, ++Z) { s = (*X1<0.0f); *Z = (s!=sp); sp = s; }
                }
                else if (going==1)
                {
                    sp = (*X1>=0.0f); X1+=2u*K; *Z++ = 0;
                    for (size_t l=1u; l<Lx; ++l, X1+=2u*K, ++Z) { s = (*X1>=0.0f); *Z = s*(s!=sp); sp = s; }
                }
                else if (going==-1)
                {
                    sp = (*X1<0.0f); X1+=2u*K; *Z++ = 0;
                    for (size_t l=1u; l<Lx; ++l, X1+=2u*K, ++Z) { s = (*X1<0.0f); *Z = s*(s!=sp); sp = s; }
                }
                Z -= Lx;

                //Initialize windowing
                w = 0u; cc = c0; cs = (int)round(c0); es = cs + (int)Lpost;

                //Windows before first samp
                while (es<0 && w<W)
                {
                    *Y = 0.0f; Y += K;
                    ++w; cc += stp; cs = (int)roundf(cc);
                    es = cs + (int)Lpost;
                }
                ss = cs - (int)Lpre;
                prev_cs = cs;

                //Windows overlapping first samp
                while (ss<0 && w<W)
                {
                    sm = 0.0f; X2 -= ss;
                    for (size_t l=(size_t)(-ss); l<Lw; ++l, ++X2, ++Z) { sm += (float)(*Z) * *X2; }
                    *Y = sm; Y += K;
                    ++w; cc += stp; cs = (int)roundf(cc);
                    Z -= (int)Lw + ss; X2 -= Lw;
                    ss = cs - (int)Lpre; prev_cs = cs;
                }
                Z += ss;
                es = cs + (int)Lpost;
                
                //Windows fully within sig
                while (es<(int)Lx && w<W)
                {
                    sm = 0.0f;
                    for (size_t l=0u; l<Lw; ++l, ++Z, ++X2) { sm += (float)*Z * *X2; }
                    *Y = sm; Y += K;
                    ++w; cc += stp; cs = (int)roundf(cc);
                    Z += cs - prev_cs - (int)Lw; X2 -= Lw;
                    es = cs + (int)Lpost; prev_cs = cs;
                }
                ss = cs - (int)Lpre;

                //Windows overlapping last samp
                while (ss<(int)Lx && w<W)
                {
                    sm = 0.0f;
                    for (size_t l=0u; l<Lx-(size_t)ss; ++l, ++Z, ++X2) { sm += (float)(*Z) * *X2; }
                    *Y = sm; Y += K;
                    X2 += ss - (int)Lx;
                    ++w; cc += stp; cs = (int)roundf(cc);
                    Z += cs - prev_cs - (int)Lx + ss;
                    ss = cs - (int)Lpre; prev_cs = cs;
                }
                Z -= cs - prev_cs + ss;

                //Windows after last samp
                while (w<W) { *Y = 0.0f; Y += K; ++w; }
            }
        }
        
        free(Z);
    }

    return 0;
}


int zcr_windowed_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const size_t dim, const double c0, const double stp, const int going)
{
    if (dim>3u) { fprintf(stderr,"error in zcr_windowed_z: dim must be in [0 3]\n"); return 1; }
    if (Lw<1u) { fprintf(stderr,"error in zcr_windowed_z: Lw (winlength) must be positive\n"); return 1; }
    if (stp<DBL_EPSILON) { fprintf(stderr,"error in zcr_windowed_z: stp (step size) must be positive\n"); return 1; }
    if (going!=0 && going!=1 && going!=-1) { fprintf(stderr,"error in zcr_windowed_z: going must be in {-1,0,1}\n"); return 1; }

    const size_t N = R*C*S*H;
    if (N==0u) {fprintf(stderr,"error in zcr_windowed_z: X1 is empty (size of at least one dim is 0)\n"); return 1; }

    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lx<1u) { fprintf(stderr,"error in zcr_windowed_z: Lx (length of vecs in X1) must be positive\n"); return 1; }
    if (Lw>Lx) { fprintf(stderr,"error in zcr_windowed_z: Lw (winlength) must be <= Lx (length of vecs in X1)\n"); return 1; }
    if (W>Lx) { fprintf(stderr,"error in zcr_windowed_z: W must be <= Lx (length of vecs in X1)\n"); return 1; }
    if (c0>(double)(Lx-1u)) { fprintf(stderr,"error in zcr_windowed_z: c0 (center samp of 1st frame) must be < Lx (length of vecs in X1)\n"); return 1; }

    if (W==0u) {}
    else
    {
        const size_t Lpre = Lw/2u;                  //nsamps before center samp
        const size_t Lpost = Lw-Lpre-1u;            //nsamps after center samp
        size_t w = 0u;                              //current frame
        double cc = c0;                              //current exact center-samp
        int cs = (int)round(c0);                   //current rounded center-samp
        int ss, es = cs + (int)Lpost;               //current rounded start-samp, end-samp
        int prev_cs;                                //previous rounded center-samp
        int *Z, s, sp;                              //ints for ZCR
        double sm;                                   //local sum for ZCR

        //Allocate
        if (!(Z=(int *)malloc(Lx*sizeof(int)))) { fprintf(stderr,"error in zcr_windowed_z: problem with malloc. "); perror("malloc"); return 1; }

        //Get block and vec sizes
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        //Comment this out to use real part of X1 for ZCs
        ++X1;

        for (size_t g=0u; g<G; ++g, X1+=2u*B*(Lx-1u), Y+=B*(W-1u))
        {
            for (size_t b=0u; b<B; ++b, X1-=2u*K*Lx-2u, Y-=K*W-1u)
            {
                //Get Z
                if (going==0)
                {
                    sp = (*X1<0.0); X1+=2u*K; *Z++ = 0;
                    for (size_t l=1u; l<Lx; ++l, X1+=2u*K, ++Z) { s = (*X1<0.0); *Z = (s!=sp); sp = s; }
                }
                else if (going==1)
                {
                    sp = (*X1>=0.0); X1+=2u*K; *Z++ = 0;
                    for (size_t l=1u; l<Lx; ++l, X1+=2u*K, ++Z) { s = (*X1>=0.0); *Z = s*(s!=sp); sp = s; }
                }
                else if (going==-1)
                {
                    sp = (*X1<0.0); X1+=2u*K; *Z++ = 0;
                    for (size_t l=1u; l<Lx; ++l, X1+=2u*K, ++Z) { s = (*X1<0.0); *Z = s*(s!=sp); sp = s; }
                }
                Z -= Lx;

                //Initialize windowing
                w = 0u; cc = c0; cs = (int)round(c0); es = cs + (int)Lpost;

                //Windows before first samp
                while (es<0 && w<W)
                {
                    *Y = 0.0; Y += K;
                    ++w; cc += stp; cs = (int)round(cc);
                    es = cs + (int)Lpost;
                }
                ss = cs - (int)Lpre;
                prev_cs = cs;

                //Windows overlapping first samp
                while (ss<0 && w<W)
                {
                    sm = 0.0; X2 -= ss;
                    for (size_t l=(size_t)(-ss); l<Lw; ++l, ++X2, ++Z) { sm += (double)(*Z) * *X2; }
                    *Y = sm; Y += K;
                    ++w; cc += stp; cs = (int)round(cc);
                    Z -= (int)Lw + ss; X2 -= Lw;
                    ss = cs - (int)Lpre; prev_cs = cs;
                }
                Z += ss;
                es = cs + (int)Lpost;
                
                //Windows fully within sig
                while (es<(int)Lx && w<W)
                {
                    sm = 0.0;
                    for (size_t l=0u; l<Lw; ++l, ++Z, ++X2) { sm += (double)*Z * *X2; }
                    *Y = sm; Y += K;
                    ++w; cc += stp; cs = (int)round(cc);
                    Z += cs - prev_cs - (int)Lw; X2 -= Lw;
                    es = cs + (int)Lpost; prev_cs = cs;
                }
                ss = cs - (int)Lpre;

                //Windows overlapping last samp
                while (ss<(int)Lx && w<W)
                {
                    sm = 0.0;
                    for (size_t l=0u; l<Lx-(size_t)ss; ++l, ++Z, ++X2) { sm += (double)(*Z) * *X2; }
                    *Y = sm; Y += K;
                    X2 += ss - (int)Lx;
                    ++w; cc += stp; cs = (int)round(cc);
                    Z += cs - prev_cs - (int)Lx + ss;
                    ss = cs - (int)Lpre; prev_cs = cs;
                }
                Z -= cs - prev_cs + ss;

                //Windows after last samp
                while (w<W) { *Y = 0.0; Y += K; ++w; }
            }
        }
        
        free(Z);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
