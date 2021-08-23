//This gets mean crossings (MCs) as usual from X1,
//and then applies window X2 to get MC rate (MCR).
//
//No internal normalization by the sum of X2 is done here:
//X2 should sum to 1 if the output Y is interpreted as a moving average.

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int mcr_windowed_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const size_t dim, const float c0, const float stp, const int going);
int mcr_windowed_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const size_t dim, const double c0, const double stp, const int going);


int mcr_windowed_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const size_t dim, const float c0, const float stp, const int going)
{
    if (dim>3u) { fprintf(stderr,"error in mcr_windowed_s: dim must be in [0 3]\n"); return 1; }
    if (Lw<1u) { fprintf(stderr,"error in mcr_windowed_s: Lw (winlength) must be positive\n"); return 1; }
    if (stp<FLT_EPSILON) { fprintf(stderr,"error in mcr_windowed_s: stp (step size) must be positive\n"); return 1; }
    if (going!=0 && going!=1 && going!=-1) { fprintf(stderr,"error in mcr_windowed_s: going must be in {-1,0,1}\n"); return 1; }

    const size_t N = R*C*S*H;
    if (N==0u) {fprintf(stderr,"error in mcr_windowed_s: X1 is empty (size of at least one dim is 0)\n"); return 1; }

    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lx<1u) { fprintf(stderr,"error in mcr_windowed_s: Lx (length of vecs in X1) must be positive\n"); return 1; }
    if (Lw>Lx) { fprintf(stderr,"error in mcr_windowed_s: Lw (winlength) must be <= Lx (length of vecs in X1)\n"); return 1; }
    if (W>Lx) { fprintf(stderr,"error in mcr_windowed_s: W must be <= Lx (length of vecs in X1)\n"); return 1; }
    if (c0>(float)(Lx-1u)) { fprintf(stderr,"error in mcr_windowed_s: c0 (center samp of 1st frame) must be < Lx (length of vecs in X1)\n"); return 1; }

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
        int *Z, s, sp;                              //ints for MCR
        float sm, mn;                               //local sum, mean for MCR

        //Allocate
        if (!(Z=(int *)malloc(Lx*sizeof(int)))) { fprintf(stderr,"error in mcr_windowed_s: problem with malloc. "); perror("malloc"); return 1; }

        //Get block and vec sizes
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            //For each vec in X1
            for (size_t v=V; v>0u; --v)
            {
                //Mean
                mn = 0.0f;
                for (size_t l=Lx; l>0u; --l, ++X1) { mn += *X1; }
                mn /= (float)Lx; X1 -= Lx;

                //Get Z
                if (going==0)
                {
                    sp = (*X1++<mn); *Z++ = 0;
                    for (size_t l=Lx; l>1u; --l, ++X1, ++Z) { s = (*X1<mn); *Z = (s!=sp); sp = s; }
                }
                else if (going==1)
                {
                    sp = (*X1++>=mn); *Z++ = 0;
                    for (size_t l=Lx; l>1u; --l, ++X1, ++Z) { s = (*X1>=mn); *Z = s*(s!=sp); sp = s; }
                }
                else if (going==-1)
                {
                    sp = (*X1++<mn); *Z++ = 0;
                    for (size_t l=Lx; l>1u; --l, ++X1, ++Z) { s = (*X1<mn); *Z = s*(s!=sp); sp = s; }
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
                    for (size_t l=Lw-(size_t)(-ss); l>0u; --l, ++X2, ++Z) { sm += (float)(*Z) * *X2; }
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
                    for (size_t l=Lw; l>0u; --l, ++Z, ++X2) { sm += (float)*Z * *X2; }
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
                    for (size_t l=Lx-(size_t)ss; l>0u; --l, ++Z, ++X2) { sm += (float)(*Z) * *X2; }
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
            for (size_t g=G; g>0u; --g, X1+=B*(Lx-1u), Y+=B*(W-1u))
            {
                for (size_t b=B; b>0u; --b, X1-=K*Lx-1u, Y-=K*W-1u)
                {
                    //Mean
                    mn = 0.0f;
                    for (size_t l=Lx; l>0u; --l, X1+=K) { mn += *X1; }
                    mn /= (float)Lx; X1 -= K*Lx;

                    //Get Z
                    if (going==0)
                    {
                        sp = (*X1<mn); X1+=K; *Z++ = 0;
                        for (size_t l=Lx; l>1u; --l, X1+=K, ++Z) { s = (*X1<mn); *Z = (s!=sp); sp = s; }
                    }
                    else if (going==1)
                    {
                        sp = (*X1>=mn); X1+=K; *Z++ = 0;
                        for (size_t l=Lx; l>1u; --l, X1+=K, ++Z) { s = (*X1>=mn); *Z = s*(s!=sp); sp = s; }
                    }
                    else if (going==-1)
                    {
                        sp = (*X1<mn); X1+=K; *Z++ = 0;
                        for (size_t l=Lx; l>1u; --l, X1+=K, ++Z) { s = (*X1<mn); *Z = s*(s!=sp); sp = s; }
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
                        for (size_t l=Lw-(size_t)(-ss); l>0u; --l, ++X2, ++Z) { sm += (float)(*Z) * *X2; }
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
                        for (size_t l=Lw; l>0u; --l, ++Z, ++X2) { sm += (float)*Z * *X2; }
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
                        for (size_t l=Lx-(size_t)ss; l>0u; --l, ++Z, ++X2) { sm += (float)(*Z) * *X2; }
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


int mcr_windowed_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const size_t dim, const double c0, const double stp, const int going)
{
    if (dim>3u) { fprintf(stderr,"error in mcr_windowed_d: dim must be in [0 3]\n"); return 1; }
    if (Lw<1u) { fprintf(stderr,"error in mcr_windowed_d: Lw (winlength) must be positive\n"); return 1; }
    if (stp<DBL_EPSILON) { fprintf(stderr,"error in mcr_windowed_d: stp (step size) must be positive\n"); return 1; }
    if (going!=0 && going!=1 && going!=-1) { fprintf(stderr,"error in mcr_windowed_d: going must be in {-1,0,1}\n"); return 1; }

    const size_t N = R*C*S*H;
    if (N==0u) {fprintf(stderr,"error in mcr_windowed_d: X1 is empty (size of at least one dim is 0)\n"); return 1; }

    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lx<1u) { fprintf(stderr,"error in mcr_windowed_d: Lx (length of vecs in X1) must be positive\n"); return 1; }
    if (Lw>Lx) { fprintf(stderr,"error in mcr_windowed_d: Lw (winlength) must be <= Lx (length of vecs in X1)\n"); return 1; }
    if (W>Lx) { fprintf(stderr,"error in mcr_windowed_d: W must be <= Lx (length of vecs in X1)\n"); return 1; }
    if (c0>(double)(Lx-1u)) { fprintf(stderr,"error in mcr_windowed_d: c0 (center samp of 1st frame) must be < Lx (length of vecs in X1)\n"); return 1; }

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
        int *Z, s, sp;                              //ints for MCR
        double sm, mn;                              //local sum, mean for MCR

        //Allocate
        if (!(Z=(int *)malloc(Lx*sizeof(int)))) { fprintf(stderr,"error in mcr_windowed_d: problem with malloc. "); perror("malloc"); return 1; }

        //Get block and vec sizes
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            //For each vec in X1
            for (size_t v=V; v>0u; --v)
            {
                //Mean
                mn = 0.0;
                for (size_t l=Lx; l>0u; --l, ++X1) { mn += *X1; }
                mn /= (double)Lx; X1 -= Lx;

                //Get Z
                if (going==0)
                {
                    sp = (*X1++<mn); *Z++ = 0;
                    for (size_t l=Lx; l>1u; --l, ++X1, ++Z) { s = (*X1<mn); *Z = (s!=sp); sp = s; }
                }
                else if (going==1)
                {
                    sp = (*X1++>=mn); *Z++ = 0;
                    for (size_t l=Lx; l>1u; --l, ++X1, ++Z) { s = (*X1>=mn); *Z = s*(s!=sp); sp = s; }
                }
                else if (going==-1)
                {
                    sp = (*X1++<mn); *Z++ = 0;
                    for (size_t l=Lx; l>1u; --l, ++X1, ++Z) { s = (*X1<mn); *Z = s*(s!=sp); sp = s; }
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
                    for (size_t l=Lw-(size_t)(-ss); l>0u; --l, ++X2, ++Z) { sm += (double)(*Z) * *X2; }
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
                    for (size_t l=Lw; l>0u; --l, ++Z, ++X2) { sm += (double)*Z * *X2; }
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
                    for (size_t l=Lx-(size_t)ss; l>0u; --l, ++Z, ++X2) { sm += (double)(*Z) * *X2; }
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
            for (size_t g=G; g>0u; --g, X1+=B*(Lx-1u), Y+=B*(W-1u))
            {
                for (size_t b=B; b>0u; --b, X1-=K*Lx-1u, Y-=K*W-1u)
                {
                    //Mean
                    mn = 0.0;
                    for (size_t l=Lx; l>0u; --l, X1+=K) { mn += *X1; }
                    mn /= (double)Lx; X1 -= K*Lx;

                    //Get Z
                    if (going==0)
                    {
                        sp = (*X1<mn); X1+=K; *Z++ = 0;
                        for (size_t l=Lx; l>1u; --l, X1+=K, ++Z) { s = (*X1<mn); *Z = (s!=sp); sp = s; }
                    }
                    else if (going==1)
                    {
                        sp = (*X1>=mn); X1+=K; *Z++ = 0;
                        for (size_t l=Lx; l>1u; --l, X1+=K, ++Z) { s = (*X1>=mn); *Z = s*(s!=sp); sp = s; }
                    }
                    else if (going==-1)
                    {
                        sp = (*X1<mn); X1+=K; *Z++ = 0;
                        for (size_t l=Lx; l>1u; --l, X1+=K, ++Z) { s = (*X1<mn); *Z = s*(s!=sp); sp = s; }
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
                        for (size_t l=Lw-(size_t)(-ss); l>0u; --l, ++X2, ++Z) { sm += (double)(*Z) * *X2; }
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
                        for (size_t l=Lw; l>0u; --l, ++Z, ++X2) { sm += (double)*Z * *X2; }
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
                        for (size_t l=Lx-(size_t)ss; l>0u; --l, ++Z, ++X2) { sm += (double)(*Z) * *X2; }
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


#ifdef __cplusplus
}
}
#endif
