//1D convolution of each vector in X1 by X2 along dim.
//Each vector in X1 has length L1. X2 has length L2.

//FIR filtering is similar, except FIR is causal and conv is non-causal.
//The overlap is controlled here by a shape param 'full', 'same', 'valid'.

//Note that some "convolution" functions actually do cross-correlation.
//For actual cross-correlation (no flip of X2), see xcorr.

//Profile note: sm = fma(X1,X2,sm) was definitely slower than sm += X1*X2.

#include <stdio.h>
#include <string.h>
#include "codee_dsp.h"
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int conv_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in conv_s: N (total num elements in X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in conv_s: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in conv_s: L2 (length X2) must be positive\n"); return 1; }

    const int inc = 1 - (int)L2;       //fixed increment for X1 below
    size_t w=0u, W;                    //current frame and total frames
    int ss, es;                        //current start-samp, end-samp
    float sm;                          //intermediate sum

    //Set ss, es, W according to shape
    if (strncmp(shape,"full",4u)==0)
    {
        es = 0; ss = 1 - (int)L2;
        W = L1 + L2 - 1u;
    }
    else if (strncmp(shape,"same",4u)==0)
    {
        es = (int)L2/2; ss = es - (int)L2 + 1;
        W = L1;
    }
    else if (strncmp(shape,"valid",5u)==0)
    {
        ss = 0; es = (int)L2 - 1;
        W = (L1<L2) ? 0u : L1-L2+1u;
    }
    else
    {
        fprintf(stderr,"error in conv_s: shape string must be 'full', 'same' or 'valid'\n"); return 1;
    }

    if (W==0u) {}
    else if (L1==N)
    {
        //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

        //X2 overlaps first samp of X1
        while (ss<0 && es<(int)L1 && w<W)
        {
            sm = 0.0f; X1 += es;
            //for (int n=es; n>0; --n, --X1, ++X2) { sm = fmaf(*X1,*X2,sm); }
            for (int n=es; n>0; --n, --X1, ++X2) { sm += *X1 * *X2; }
            *Y++ = sm + *X1**X2;
            X2 -= es;
            ++ss; ++es; ++w;
        }
        X2 += es;

        if (L2>L1)  //X1 fully within X2
        {
            while (ss<0 && w<W)
            {
                sm = 0.0f;
                for (size_t l=L1; l>0u; --l, ++X1, --X2) { sm += *X1 * *X2; }
                *Y++ = sm;
                X1 -= L1; X2 += L1+1u;
                ++ss; ++w;
            }
            es = ss + (int)L2 - 1;
        }
        else        //X2 fully within X1
        {        
            while (es<(int)L1 && w<W)
            {
                sm = 0.0f;
                for (size_t l=L2; l>0u; --l, ++X1, --X2) { sm += *X1 * *X2; }
                *Y++ = sm;
                X1 += inc; X2 += L2;
                ++es; ++w;
            }
            ss = es - (int)L2 + 1;
        }

        //X2 overlaps end samp of X1
        while (ss<(int)L1 && w<W)
        {
            sm = 0.0f;
            for (int n=ss; n<(int)L1; ++n, ++X1, --X2) { sm += *X1 * *X2; }
            *Y++ = sm;
            X1 += 1 - (int)L1 + ss; X2 += (int)L1 - ss;
            ++ss; ++es; ++w;
        }

        //clock_gettime(CLOCK_REALTIME,&toc);
        //fprintf(stderr,"elapsed time = %.6f ms\n",(double)(toc.tv_sec-tic.tv_sec)*1e3+(double)(toc.tv_nsec-tic.tv_nsec)/1e6);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L1, G = V/B;

        for (size_t g=G; g>0u; --g, X1+=B*(L1-1u), Y+=B*(W-1u))
        {
            for (size_t b=B; b>0u; --b, ++X1, Y-=K*W-1u)
            {
                //X2 overlaps first samp of X1
                while (ss<0 && es<(int)L1 && w<W)
                {
                    sm = 0.0f; X1 += es*(int)K;
                    for (int n=es; n>0; --n, X1-=K, ++X2) { sm += *X1 * *X2; }
                    *Y = sm + *X1**X2; Y += K;
                    X2 -= es;
                    ++ss; ++es; ++w;
                }
                X2 += es;

                if (L2>L1)  //X1 fully within X2
                {
                    while (ss<0 && w<W)
                    {
                        sm = 0.0f;
                        for (size_t l=L1; l>0u; --l, X1+=K, --X2) { sm += *X1 * *X2; }
                        *Y = sm; Y += K;
                        X1 -= K*L1; X2 += L1+1u;
                        ++ss; ++w;
                    }
                    es = ss + (int)L2 - 1;
                }
                else        //X2 fully within X1
                {        
                    while (es<(int)L1 && w<W)
                    {
                        sm = 0.0f;
                        for (size_t l=L2; l>0u; --l, X1+=K, --X2) { sm += *X1 * *X2; }
                        *Y = sm; Y += K;
                        X1 += inc*(int)K; X2 += L2;
                        ++es; ++w;
                    }
                    ss = es - (int)L2 + 1;
                }

                //Frames overlapping end samp
                while (ss<(int)L1 && w<W)
                {
                    sm = 0.0f;
                    for (int n=ss; n<(int)L1; ++n, X1+=K, --X2) { sm += *X1 * *X2; }
                    *Y = sm; Y += K;
                    X1 += (int)K*(1-(int)L1+ss);
                    X2 += (int)L1 - ss;
                    ++ss; ++w;
                }
                X1 -= (int)K*ss; X2 -= L2 - 1u;
                ss -= (int)W; es = ss + (int)L2 - 1; w = 0u;
            }
        }
    }

    return 0;
}


int conv_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv_d: dim must be in [0 3]\n"); return 1; }
    
    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in conv_d: N (total num elements in X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in conv_d: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in conv_d: L2 (length X2) must be positive\n"); return 1; }

    const int inc = 1 - (int)L2;       //fixed increment for X1 below
    size_t w=0u, W;                    //current frame and total frames
    int ss, es;                        //current start-samp, end-samp
    double sm;                         //intermediate sum

    //Set ss, es, W according to shape
    if (strncmp(shape,"full",4u)==0)
    {
        es = 0; ss = 1 - (int)L2;
        W = L1 + L2 - 1u;
    }
    else if (strncmp(shape,"same",4u)==0)
    {
        es = (int)L2/2; ss = es - (int)L2 + 1;
        W = L1;
    }
    else if (strncmp(shape,"valid",5u)==0)
    {
        ss = 0; es = (int)L2 - 1;
        W = (L1<L2) ? 0u : L1-L2+1u;
    }
    else
    {
        fprintf(stderr,"error in conv_d: shape string must be 'full', 'same' or 'valid'\n"); return 1;
    }

    if (W==0u) {}
    else if (L1==N)
    {
        //X2 overlaps first samp of X1
        while (ss<0 && es<(int)L1 && w<W)
        {
            sm = 0.0; X1 += es;
            for (int n=es; n>0; --n, --X1, ++X2) { sm += *X1 * *X2; }
            *Y++ = sm + *X1**X2;
            X2 -= es;
            ++ss; ++es; ++w;
        }
        X2 += es;

        if (L2>L1)  //X1 fully within X2
        {
            while (ss<0 && w<W)
            {
                sm = 0.0;
                for (size_t l=L1; l>0u; --l, ++X1, --X2) { sm += *X1 * *X2; }
                *Y++ = sm;
                X1 -= L1; X2 += L1+1u;
                ++ss; ++w;
            }
            es = ss + (int)L2 - 1;
        }
        else        //X2 fully within X1
        {        
            while (es<(int)L1 && w<W)
            {
                sm = 0.0;
                for (size_t l=L2; l>0u; --l, ++X1, --X2) { sm += *X1 * *X2; }
                *Y++ = sm;
                X1 += inc; X2 += L2;
                ++es; ++w;
            }
            ss = es - (int)L2 + 1;
        }

        //X2 overlaps end samp of X1
        while (ss<(int)L1 && w<W)
        {
            sm = 0.0;
            for (int n=ss; n<(int)L1; ++n, ++X1, --X2) { sm += *X1 * *X2; }
            *Y++ = sm;
            X1 += 1 - (int)L1 + ss; X2 += (int)L1 - ss;
            ++ss; ++es; ++w;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L1, G = V/B;

        for (size_t g=G; g>0u; --g, X1+=B*(L1-1u), Y+=B*(W-1u))
        {
            for (size_t b=B; b>0u; --b, ++X1, Y-=K*W-1u)
            {
                //X2 overlaps first samp of X1
                while (ss<0 && es<(int)L1 && w<W)
                {
                    sm = 0.0; X1 += es*(int)K;
                    for (int n=es; n>0; --n, X1-=K, ++X2) { sm += *X1 * *X2; }
                    *Y = sm + *X1**X2; Y += K;
                    X2 -= es;
                    ++ss; ++es; ++w;
                }
                X2 += es;

                if (L2>L1)  //X1 fully within X2
                {
                    while (ss<0 && w<W)
                    {
                        sm = 0.0;
                        for (size_t l=L1; l>0u; --l, X1+=K, --X2) { sm += *X1 * *X2; }
                        *Y = sm; Y += K;
                        X1 -= K*L1; X2 += L1+1u;
                        ++ss; ++w;
                    }
                    es = ss + (int)L2 - 1;
                }
                else        //X2 fully within X1
                {        
                    while (es<(int)L1 && w<W)
                    {
                        sm = 0.0;
                        for (size_t l=L2; l>0u; --l, X1+=K, --X2) { sm += *X1 * *X2; }
                        *Y = sm; Y += K;
                        X1 += inc*(int)K; X2 += L2;
                        ++es; ++w;
                    }
                    ss = es - (int)L2 + 1;
                }

                //Frames overlapping end samp
                while (ss<(int)L1 && w<W)
                {
                    sm = 0.0;
                    for (int n=ss; n<(int)L1; ++n, X1+=K, --X2) { sm += *X1 * *X2; }
                    *Y = sm; Y += K;
                    X1 += (int)K*(1-(int)L1+ss);
                    X2 += (int)L1 - ss;
                    ++ss; ++w;
                }
                X1 -= (int)K*ss; X2 -= L2 - 1u;
                ss -= (int)W; es = ss + (int)L2 - 1; w = 0u;
            }
        }
    }

    return 0;
}


int conv_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in conv_c: N (total num elements in X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in conv_c: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in conv_c: L2 (length X2) must be positive\n"); return 1; }

    const int inc = 2 - 2*(int)L2;     //fixed increment for X1 below
    size_t w=0u, W;                    //current frame and total frames
    int ss, es;                        //current start-samp, end-samp
    float smr, smi;                    //intermediate sums

    //Set ss, es, W according to shape
    if (strncmp(shape,"full",4u)==0)
    {
        es = 0; ss = 1 - (int)L2;
        W = L1 + L2 - 1u;
    }
    else if (strncmp(shape,"same",4u)==0)
    {
        es = (int)L2/2; ss = es - (int)L2 + 1;
        W = L1;
    }
    else if (strncmp(shape,"valid",5u)==0)
    {
        ss = 0; es = (int)L2 - 1;
        W = (L1<L2) ? 0u : L1-L2+1u;
    }
    else
    {
        fprintf(stderr,"error in conv_c: shape string must be 'full', 'same' or 'valid'\n"); return 1;
    }

    if (W==0u) {}
    else if (L1==N)
    {
        //X2 overlaps first samp of X1
        while (ss<0 && es<(int)L1 && w<W)
        {
            smr = smi = 0.0f; X1 += 2*es;
            for (int n=es; n>0; --n, X1-=2, X2+=2)
            {
                smr += *X1**X2 - *(X1+1)**(X2+1);
                smi += *X1**(X2+1) + *(X1+1)**X2;
            }
            *Y++ = smr + *X1**X2 - *(X1+1)**(X2+1);
            *Y++ = smi + *X1**(X2+1) + *(X1+1)**X2;
            X2 -= 2*es;
            ++ss; ++es; ++w;
        }
        X2 += 2*es;

        if (L2>L1)  //X1 fully within X2
        {
            while (ss<0 && w<W)
            {
                smr = smi = 0.0f;
                for (size_t l=L1; l>0u; --l, X1+=2, X2-=2)
                {
                    smr += *X1**X2 - *(X1+1)**(X2+1);
                    smi += *X1**(X2+1) + *(X1+1)**X2;
                }
                *Y++ = smr; *Y++ = smi;
                X1 -= 2u*L1; X2 += 2u*L1+2u;
                ++ss; ++w;
            }
            es = ss + (int)L2 - 1;
        }
        else        //X2 fully within X1
        {        
            while (es<(int)L1 && w<W)
            {
                smr = smi = 0.0f;
                for (size_t l=L2; l>0u; --l, X1+=2, X2-=2)
                {
                    smr += *X1**X2 - *(X1+1)**(X2+1);
                    smi += *X1**(X2+1) + *(X1+1)**X2;
                }
                *Y++ = smr; *Y++ = smi;
                X1 += inc; X2 += 2u*L2;
                ++es; ++w;
            }
            ss = es - (int)L2 + 1;
        }

        //X2 overlaps end samp of X1
        while (ss<(int)L1 && w<W)
        {
            smr = smi = 0.0f;
            for (int n=ss; n<(int)L1; ++n, X1+=2, X2-=2)
            {
                smr += *X1**X2 - *(X1+1)**(X2+1);
                smi += *X1**(X2+1) + *(X1+1)**X2;
            }
            *Y++ = smr; *Y++ = smi;
            X1 += 2*(1-(int)L1+ss); X2 += 2*((int)L1-ss);
            ++ss; ++es; ++w;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L1, G = V/B;

        for (size_t g=G; g>0u; --g, X1+=2u*B*(L1-1u), Y+=2u*B*(W-1u))
        {
            for (size_t b=B; b>0u; --b, X1+=2, Y-=2u*K*W-2u)
            {
                //X2 overlaps first samp of X1
                while (ss<0 && es<(int)L1 && w<W)
                {
                    smr = smi = 0.0f; X1 += 2*es*(int)K;
                    for (int n=es; n>0; --n, X1-=2u*K, X2+=2)
                    {
                        smr += *X1**X2 - *(X1+1)**(X2+1);
                        smi += *X1**(X2+1) + *(X1+1)**X2;
                    }
                    *Y = smr + *X1**X2 - *(X1+1)**(X2+1);
                    *(Y+1) = smi + *X1**(X2+1) + *(X1+1)**X2;
                    Y += 2u*K;
                    X2 -= 2*es;
                    ++ss; ++es; ++w;
                }
                X2 += 2*es;

                if (L2>L1)  //X1 fully within X2
                {
                    while (ss<0 && w<W)
                    {
                        smr = smi = 0.0f;
                        for (size_t l=L1; l>0u; --l, X1+=2u*K, X2-=2)
                        {
                            smr += *X1**X2 - *(X1+1)**(X2+1);
                            smi += *X1**(X2+1) + *(X1+1)**X2;
                        }
                        *Y = smr; *(Y+1) = smi; Y += 2u*K;
                        X1 -= 2u*K*L1; X2 += 2u*L1 + 2u;
                        ++ss; ++w;
                    }
                    es = ss + (int)L2 - 1;
                }
                else        //X2 fully within X1
                {        
                    while (es<(int)L1 && w<W)
                    {
                        smr = smi = 0.0f;
                        for (size_t l=L2; l>0u; --l, X1+=2u*K, X2-=2)
                        {
                            smr += *X1**X2 - *(X1+1)**(X2+1);
                            smi += *X1**(X2+1) + *(X1+1)**X2;
                        }
                        *Y = smr; *(Y+1) = smi; Y += 2u*K;
                        X1 += inc*(int)K; X2 += 2u*L2;
                        ++es; ++w;
                    }
                    ss = es - (int)L2 + 1;
                }

                //Frames overlapping end samp
                while (ss<(int)L1 && w<W)
                {
                    smr = smi = 0.0f;
                    for (int n=ss; n<(int)L1; ++n, X1+=2u*K, X2-=2)
                    {
                        smr += *X1**X2 - *(X1+1)**(X2+1);
                        smi += *X1**(X2+1) + *(X1+1)**X2;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    X1 += 2*(int)K*(1-(int)L1+ss);
                    X2 += 2*((int)L1-ss);
                    ++ss; ++w;
                }
                X1 -= 2*(int)K*ss; X2 -= 2u*L2 - 2u;
                ss -= (int)W; es = ss + (int)L2 - 1; w = 0u;
            }
        }
    }

    return 0;
}


int conv_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in conv_z: N (total num elements in X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in conv_z: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in conv_z: L2 (length X2) must be positive\n"); return 1; }

    const int inc = 2 - 2*(int)L2;     //fixed increment for X1 below
    size_t w=0u, W;                    //current frame and total frames
    int ss, es;                        //current start-samp, end-samp
    double smr, smi;                   //intermediate sums

    //Set ss, es, W according to shape
    if (strncmp(shape,"full",4u)==0)
    {
        es = 0; ss = 1 - (int)L2;
        W = L1 + L2 - 1u;
    }
    else if (strncmp(shape,"same",4u)==0)
    {
        es = (int)L2/2; ss = es - (int)L2 + 1;
        W = L1;
    }
    else if (strncmp(shape,"valid",5u)==0)
    {
        ss = 0; es = (int)L2 - 1;
        W = (L1<L2) ? 0u : L1-L2+1u;
    }
    else
    {
        fprintf(stderr,"error in conv_z: shape string must be 'full', 'same' or 'valid'\n"); return 1;
    }

    if (W==0u) {}
    else if (L1==N)
    {
        //X2 overlaps first samp of X1
        while (ss<0 && es<(int)L1 && w<W)
        {
            smr = smi = 0.0; X1 += 2*es;
            for (int n=es; n>0; --n, X1-=2, X2+=2)
            {
                smr += *X1**X2 - *(X1+1)**(X2+1);
                smi += *X1**(X2+1) + *(X1+1)**X2;
            }
            *Y++ = smr + *X1**X2 - *(X1+1)**(X2+1);
            *Y++ = smi + *X1**(X2+1) + *(X1+1)**X2;
            X2 -= 2*es;
            ++ss; ++es; ++w;
        }
        X2 += 2*es;

        if (L2>L1)  //X1 fully within X2
        {
            while (ss<0 && w<W)
            {
                smr = smi = 0.0;
                for (size_t l=L1; l>0u; --l, X1+=2, X2-=2)
                {
                    smr += *X1**X2 - *(X1+1)**(X2+1);
                    smi += *X1**(X2+1) + *(X1+1)**X2;
                }
                *Y++ = smr; *Y++ = smi;
                X1 -= 2u*L1; X2 += 2u*L1+2u;
                ++ss; ++w;
            }
            es = ss + (int)L2 - 1;
        }
        else        //X2 fully within X1
        {        
            while (es<(int)L1 && w<W)
            {
                smr = smi = 0.0;
                for (size_t l=L2; l>0u; --l, X1+=2, X2-=2)
                {
                    smr += *X1**X2 - *(X1+1)**(X2+1);
                    smi += *X1**(X2+1) + *(X1+1)**X2;
                }
                *Y++ = smr; *Y++ = smi;
                X1 += inc; X2 += 2u*L2;
                ++es; ++w;
            }
            ss = es - (int)L2 + 1;
        }

        //X2 overlaps end samp of X1
        while (ss<(int)L1 && w<W)
        {
            smr = smi = 0.0;
            for (int n=ss; n<(int)L1; ++n, X1+=2, X2-=2)
            {
                smr += *X1**X2 - *(X1+1)**(X2+1);
                smi += *X1**(X2+1) + *(X1+1)**X2;
            }
            *Y++ = smr; *Y++ = smi;
            X1 += 2*(1-(int)L1+ss); X2 += 2*((int)L1-ss);
            ++ss; ++es; ++w;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L1, G = V/B;

        for (size_t g=G; g>0u; --g, X1+=2u*B*(L1-1u), Y+=2u*B*(W-1u))
        {
            for (size_t b=B; b>0u; --b, X1+=2, Y-=2u*K*W-2u)
            {
                //X2 overlaps first samp of X1
                while (ss<0 && es<(int)L1 && w<W)
                {
                    smr = smi = 0.0; X1 += 2*es*(int)K;
                    for (int n=es; n>0; --n, X1-=2u*K, X2+=2)
                    {
                        smr += *X1**X2 - *(X1+1)**(X2+1);
                        smi += *X1**(X2+1) + *(X1+1)**X2;
                    }
                    *Y = smr + *X1**X2 - *(X1+1)**(X2+1);
                    *(Y+1) = smi + *X1**(X2+1) + *(X1+1)**X2;
                    Y += 2u*K;
                    X2 -= 2*es;
                    ++ss; ++es; ++w;
                }
                X2 += 2*es;

                if (L2>L1)  //X1 fully within X2
                {
                    while (ss<0 && w<W)
                    {
                        smr = smi = 0.0;
                        for (size_t l=L1; l>0u; --l, X1+=2u*K, X2-=2)
                        {
                            smr += *X1**X2 - *(X1+1)**(X2+1);
                            smi += *X1**(X2+1) + *(X1+1)**X2;
                        }
                        *Y = smr; *(Y+1) = smi; Y += 2u*K;
                        X1 -= 2u*K*L1; X2 += 2u*L1 + 2u;
                        ++ss; ++w;
                    }
                    es = ss + (int)L2 - 1;
                }
                else        //X2 fully within X1
                {        
                    while (es<(int)L1 && w<W)
                    {
                        smr = smi = 0.0;
                        for (size_t l=L2; l>0u; --l, X1+=2u*K, X2-=2)
                        {
                            smr += *X1**X2 - *(X1+1)**(X2+1);
                            smi += *X1**(X2+1) + *(X1+1)**X2;
                        }
                        *Y = smr; *(Y+1) = smi; Y += 2u*K;
                        X1 += inc*(int)K; X2 += 2u*L2;
                        ++es; ++w;
                    }
                    ss = es - (int)L2 + 1;
                }

                //Frames overlapping end samp
                while (ss<(int)L1 && w<W)
                {
                    smr = smi = 0.0;
                    for (int n=ss; n<(int)L1; ++n, X1+=2u*K, X2-=2)
                    {
                        smr += *X1**X2 - *(X1+1)**(X2+1);
                        smi += *X1**(X2+1) + *(X1+1)**X2;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    X1 += 2*(int)K*(1-(int)L1+ss);
                    X2 += 2*((int)L1-ss);
                    ++ss; ++w;
                }
                X1 -= 2*(int)K*ss; X2 -= 2u*L2 - 2u;
                ss -= (int)W; es = ss + (int)L2 - 1; w = 0u;
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
