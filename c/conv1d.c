//1D convolution by X2 of each vector in X1 along dim.
//Each vector in X1 has length L1. X2 has length L2.

//FIR filtering is similar, except FIR is causal and conv is non-causal.
//Note that some "convolution" functions actually do cross-correlation.
//For actual cross-correlation (no flip of X2), see xcorr1d.

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int conv1d_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim);
int conv1d_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim);
int conv1d_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim);
int conv1d_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim);


int conv1d_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv1d_s: dim must be in [0 3]\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in conv1d_s: stp (step size) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in conv1d_s: dil (dilation) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in conv1d_s: N (total length of X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in conv1d_s: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if (c0>=(int)L1) { fprintf(stderr,"error in conv1d_s: c0 (center samp of 1st frame) must be < L1 (length of vecs in X1)\n"); return 1; }
    //if (W>N) { fprintf(stderr,"error in conv1d_s: W must be <= N (length X1)\n"); return 1; }

    if (W==0u) {}
    else if (L1==N)
    {
        const int Nb = (int)(dil*(L2-1u)) + 1;          //full length of X2 including dil
        const int inc = (int)stp - (int)(dil*L2);       //fixed increment for X1 below
        size_t w = 0u;                                  //current frame
        int ss, es = c0 + (int)(dil*(L2/2u));           //current start-samp, end-samp
        float sm;                                       //intermediate sum

        //Frames before first samp
        while (es<0 && w<W) { *Y++ = 0.0f; ++w; es += stp; }
        ss = es - Nb + 1;

        //Frames overlapping first samp
        while (ss<0 && es<(int)L1 && w<W)
        {
            sm = 0.0f; X1 += es;
            for (int n=es; n>=0; n-=(int)dil, X1-=dil, ++X2) { sm = fmaf(*X2,*X1,sm); }
            *Y++ = sm;
            X2 -= 1 + es/(int)dil; X1 += dil - (size_t)es%dil;
            ++w; ss += stp; es += stp;
        }
        X2 += L2 - 1u;

        //In case of L2>L1
        es = (int)L1 - 1;
        while (ss<0 && w<W)
        {
            sm = 0.0f; X1 += es; X2 -= L2 - L1;
            for (int n=es; n>=0; n-=(int)dil, X1-=dil, ++X2) { sm = fmaf(*X2,*X1,sm); }
            *Y++ = sm;
            X2 -= 1 + es/(int)dil; X1 += dil - (size_t)es%dil;
            ++w; ss += stp;
        }
        X2 += L2 - 1u; X1 += ss;

        //Frames fully within sig
        while (es<(int)L1 && w<W)
        {
            sm = 0.0f;
            for (size_t l=0u; l<L2; ++l, X1+=dil, --X2) { sm = fmaf(*X2,*X1,sm); }
            *Y++ = sm;
            X2 += L2; X1 += inc;
            ++w; es += stp;
        }
        ss = es - Nb + 1;

        //Frames overlapping end samp
        while (ss<(int)L1 && w<W)
        {
            size_t l = 0u; sm = 0.0f;
            for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=dil, --X2) { sm = fmaf(*X2,*X1,sm); }
            *Y++ = sm;
            X2 += l; X1 += (int)stp - (int)(dil*l);
            ++w; ss += stp;
        }

        //Frames after end samp
        while (w<W) { *Y++ = 0.0f; ++w; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L1, G = V/B;

        const int Nb = (int)(dil*(L2-1u));                  //full length of X2 including dil
        const int inc = (int)K*((int)stp-(int)(dil*L2));    //fixed increment for X1 below
        size_t w;                                           //current frame
        int ss, es;                                         //current start-samp, end-samp
        float sm;                                           //intermediate sum

        for (size_t g=0u; g<G; ++g, X1+=B*(L1-1u), Y+=B*(W-1u))
        {
            for (size_t b=0u; b<B; ++b, X1-=K*L1-1u, Y-=K*W-1u)
            {
                //Start frames
                w = 0u; es = c0 + (int)(dil*(L2/2u));

                //Frames before first samp
                while (es<0 && w<W) { *Y = 0.0f; Y+=K; ++w; es += stp; }
                ss = es - Nb + 1;

                //Frames overlapping first samp
                while (ss<0 && w<W)
                {
                    sm = 0.0f; X1 += (int)K*es;
                    for (int n=es; n>=0; n-=(int)dil, X1-=K*dil, ++X2) { sm = fmaf(*X2,*X1,sm); }
                    *Y = sm; Y += K;
                    X2 -= 1 + es/(int)dil; X1 += K*(dil-(size_t)es%dil);
                    ++w; ss += stp; es += stp;
                }
                X2 += L2 - 1u; X1 += ss*(int)K;

                //Frames fully within sig
                while (es<(int)L1 && w<W)
                {
                    sm = 0.0f;
                    for (size_t l=0u; l<L2; ++l, X1+=K*dil, --X2) { sm = fmaf(*X2,*X1,sm); }
                    *Y = sm; Y += K;
                    X2 += L2; X1 += inc;
                    ++w; es += stp;
                }
                ss = es - Nb + 1;

                //Frames overlapping end samp
                while (ss<(int)L1 && w<W)
                {
                    size_t l = 0u; sm = 0.0f;
                    for (int n=ss; n<(int)L1; n+=dil, X1+=K*dil, --X2, ++l) { sm = fmaf(*X2,*X1,sm); }
                    *Y = sm; Y += K;
                    X2 += l; X1 += (int)K*((int)stp-(int)(dil*l));
                    ++w; ss += stp;
                }
                X2 -= L2 - 1u; X1 -= K*((size_t)ss-L1);

                //Frames after end samp
                while (w<W) { *Y = 0.0f; Y+=K; ++w; }
            }
        }
    }

    return 0;
}


int conv1d_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv1d_d: dim must be in [0 3]\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in conv1d_d: stp (step size) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in conv1d_d: dil (dilation) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in conv1d_d: N (total length of X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in conv1d_d: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if (c0>=(int)L1) { fprintf(stderr,"error in conv1d_d: c0 (center samp of 1st frame) must be < L1 (length of vecs in X1)\n"); return 1; }

    if (W==0u) {}
    else if (L1==N)
    {
        const int Nb = (int)(dil*(L2-1u)) + 1;          //full length of X2 including dil
        const int inc = (int)stp - (int)(dil*L2);       //fixed increment for X1 below
        size_t w = 0u;                                  //current frame
        int ss, es = c0 + (int)(dil*(L2/2u));           //current start-samp, end-samp
        double sm;                                      //intermediate sum

        //Frames before first samp
        while (es<0 && w<W) { *Y++ = 0.0; ++w; es += stp; }
        ss = es - Nb + 1;

        //Frames overlapping first samp
        while (ss<0 && w<W)
        {
            sm = 0.0; X1 += es;
            for (int n=es; n>=0; n-=(int)dil, X1-=dil, ++X2) { sm = fma(*X2,*X1,sm); }
            *Y++ = sm;
            X2 -= 1 + es/(int)dil; X1 += dil - (size_t)es%dil;
            ++w; ss += stp; es += stp;
        }
        X2 += L2 - 1u; X1 += ss;

        //Frames fully within sig
        while (es<(int)L1 && w<W)
        {
            sm = 0.0;
            for (size_t l=0u; l<L2; ++l, X1+=dil, --X2) { sm = fma(*X2,*X1,sm); }
            *Y++ = sm;
            X2 += L2; X1 += inc;
            ++w; es += stp;
        }
        ss = es - Nb + 1;

        //Frames overlapping end samp
        while (ss<(int)L1 && w<W)
        {
            size_t l = 0u; sm = 0.0;
            for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=dil, --X2) { sm = fma(*X2,*X1,sm); }
            *Y++ = sm;
            X2 += l; X1 += (int)stp - (int)(dil*l);
            ++w; ss += stp;
        }

        //Frames after end samp
        while (w<W) { *Y++ = 0.0; ++w; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L1, G = V/B;

        const int Nb = (int)(dil*(L2-1u));                  //full length of X2 including dil
        const int inc = (int)K*((int)stp-(int)(dil*L2));    //fixed increment for X1 below
        size_t w;                                           //current frame
        int ss, es;                                         //current start-samp, end-samp
        double sm;                                          //intermediate sum

        for (size_t g=0u; g<G; ++g, X1+=B*(L1-1u), Y+=B*(W-1u))
        {
            for (size_t b=0u; b<B; ++b, X1-=K*L1-1u, Y-=K*W-1u)
            {
                //Start frames
                w = 0u; es = c0 + (int)(dil*(L2/2u));

                //Frames before first samp
                while (es<0 && w<W) { *Y = 0.0; Y+=K; ++w; es += stp; }
                ss = es - Nb + 1;

                //Frames overlapping first samp
                while (ss<0 && w<W)
                {
                    sm = 0.0; X1 += (int)K*es;
                    for (int n=es; n>=0; n-=(int)dil, X1-=K*dil, ++X2) { sm = fma(*X2,*X1,sm); }
                    *Y = sm; Y += K;
                    X2 -= 1 + es/(int)dil; X1 += K*(dil-(size_t)es%dil);
                    ++w; ss += stp; es += stp;
                }
                X2 += L2 - 1u; X1 += ss*(int)K;

                //Frames fully within sig
                while (es<(int)L1 && w<W)
                {
                    sm = 0.0;
                    for (size_t l=0u; l<L2; ++l, X1+=K*dil, --X2) { sm = fma(*X2,*X1,sm); }
                    *Y = sm; Y += K;
                    X2 += L2; X1 += inc;
                    ++w; es += stp;
                }
                ss = es - Nb + 1;

                //Frames overlapping end samp
                while (ss<(int)L1 && w<W)
                {
                    size_t l = 0u; sm = 0.0;
                    for (int n=ss; n<(int)L1; n+=dil, X1+=K*dil, --X2, ++l) { sm = fma(*X2,*X1,sm); }
                    *Y = sm; Y += K;
                    X2 += l; X1 += (int)K*((int)stp-(int)(dil*l));
                    ++w; ss += stp;
                }
                X2 -= L2 - 1u; X1 -= K*((size_t)ss-L1);

                //Frames after end samp
                while (w<W) { *Y = 0.0; Y+=K; ++w; }
            }
        }
    }

    return 0;
}


int conv1d_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv1d_c: dim must be in [0 3]\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in conv1d_c: stp (step size) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in conv1d_c: dil (dilation) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in conv1d_c: N (total length of X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in conv1d_c: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if (c0>=(int)L1) { fprintf(stderr,"error in conv1d_c: c0 (center samp of 1st frame) must be < L1 (length of vecs in X1)\n"); return 1; }

    if (W==0u) {}
    else if (L1==N)
    {
        const int Nb = (int)(dil*(L2-1u)) + 1;          //full length of X2 including dil
        const int inc = 2*((int)stp-(int)(dil*L2));     //fixed increment for X1 below
        size_t w = 0u;                                  //current frame
        int ss, es = c0 + (int)(dil*(L2/2u));           //current start-samp, end-samp
        float smr, smi;                                 //intermediate sums

        //Frames before first samp
        while (es<0 && w<W) { *Y++ = 0.0f; *Y++ = 0.0f; ++w; es += stp; }
        ss = es - Nb + 1;

        //Frames overlapping first samp
        while (ss<0 && w<W)
        {
            smr = smi = 0.0f; X1 += 2*es;
            for (int n=es; n>=0; n-=(int)dil, X1-=2u*dil, X2+=2)
            {
                smr += *X2**X1 - *(X2+1)**(X1+1);
                smi += *X2**(X1+1) + *(X2+1)**X1;
            }
            *Y++ = smr; *Y++ = smi;
            X2 -= 2*(1+es/(int)dil);
            X1 += 2u*(dil-(size_t)es%dil);
            ++w; ss += stp; es += stp;
        }
        X2 += 2u*L2 - 2u; X1 += 2*ss;

        //Frames fully within sig
        while (es<(int)L1 && w<W)
        {
            smr = smi = 0.0f;
            for (size_t l=0u; l<L2; ++l, X1+=2u*dil, X2-=2)
            {
                smr += *X2**X1 - *(X2+1)**(X1+1);
                smi += *X2**(X1+1) + *(X2+1)**X1;
            }
            *Y++ = smr; *Y++ = smi;
            X2 += 2u*L2; X1 += inc;
            ++w; es += stp;
        }
        ss = es - Nb + 1;

        //Frames overlapping end samp
        while (ss<(int)L1 && w<W)
        {
            size_t l = 0u; smr = smi = 0.0f;
            for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=2u*dil, X2-=2)
            {
                smr += *X2**X1 - *(X2+1)**(X1+1);
                smi += *X2**(X1+1) + *(X2+1)**X1;
            }
            *Y++ = smr; *Y++ = smi;
            X2 += 2u*l; X1 += 2*((int)stp-(int)(dil*l));
            ++w; ss += stp;
        }

        //Frames after end samp
        while (w<W) { *Y++ = 0.0f; *Y++ = 0.0f; ++w; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L1, G = V/B;

        const int Nb = (int)(dil*(L2-1u));                  //full length of X2 including dil
        const int inc = 2*(int)K*((int)stp-(int)(dil*L2));  //fixed increment for X1 below
        size_t w;                                           //current frame
        int ss, es;                                         //current start-samp, end-samp
        float smr, smi;                                     //intermediate sums

        for (size_t g=0u; g<G; ++g, X1+=2u*B*(L1-1u), Y+=2u*B*(W-1u))
        {
            for (size_t b=0u; b<B; ++b, X1-=2u*K*L1-2u, Y-=2u*K*W-2u)
            {
                //Start frames
                w = 0u; es = c0 + (int)(dil*(L2/2u));

                //Frames before first samp
                while (es<0 && w<W) { *Y = *(Y+1) = 0.0f; Y+=2u*K; ++w; es += stp; }
                ss = es - Nb + 1;

                //Frames overlapping first samp
                while (ss<0 && w<W)
                {
                    smr = smi = 0.0f; X1 += 2*(int)K*es;
                    for (int n=es; n>=0; n-=(int)dil, X1-=2u*K*dil, X2+=2)
                    {
                        smr += *X2**X1 - *(X2+1)**(X1+1);
                        smi += *X2**(X1+1) + *(X2+1)**X1;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    X2 -= 2*(1+es/(int)dil);
                    X1 += 2u*K*(dil-(size_t)es%dil);
                    ++w; ss += stp; es += stp;
                }
                X2 += 2u*L2 - 2u; X1 += 2*ss*(int)K;

                //Frames fully within sig
                while (es<(int)L1 && w<W)
                {
                    smr = smi = 0.0f;
                    for (size_t l=0u; l<L2; ++l, X1+=2u*K*dil, X2-=2)
                    {
                        smr += *X2**X1 - *(X2+1)**(X1+1);
                        smi += *X2**(X1+1) + *(X2+1)**X1;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    X2 += 2u*L2; X1 += inc;
                    ++w; es += stp;
                }
                ss = es - Nb + 1;

                //Frames overlapping end samp
                while (ss<(int)L1 && w<W)
                {
                    size_t l = 0u; smr = smi = 0.0f;
                    for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=2u*K*dil, X2-=2)
                    {
                        smr += *X2**X1 - *(X2+1)**(X1+1);
                        smi += *X2**(X1+1) + *(X2+1)**X1;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    X2 += 2u*l; X1 += 2*(int)K*((int)stp-(int)(dil*l));
                    ++w; ss += stp;
                }
                X2 -= 2u*L2 - 2u; X1 -= 2u*K*((size_t)ss-L1);

                //Frames after end samp
                while (w<W) { *Y = *(Y+1) = 0.0f; Y+=2u*K; ++w; }
            }
        }
    }

    return 0;
}


int conv1d_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv1d_z: dim must be in [0 3]\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in conv1d_z: stp (step size) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in conv1d_z: dil (dilation) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in conv1d_z: N (total length of X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in conv1d_z: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if (c0>=(int)L1) { fprintf(stderr,"error in conv1d_z: c0 (center samp of 1st frame) must be < L1 (length of vecs in X1)\n"); return 1; }

    if (W==0u) {}
    else if (L1==N)
    {
        const int Nb = (int)(dil*(L2-1u)) + 1;          //full length of X2 including dil
        const int inc = 2*((int)stp-(int)(dil*L2));     //fixed increment for X1 below
        size_t w = 0u;                                  //current frame
        int ss, es = c0 + (int)(dil*(L2/2u));           //current start-samp, end-samp
        double smr, smi;                                //intermediate sums

        //Frames before first samp
        while (es<0 && w<W) { *Y++ = 0.0; *Y++ = 0.0; ++w; es += stp; }
        ss = es - Nb + 1;

        //Frames overlapping first samp
        while (ss<0 && w<W)
        {
            smr = smi = 0.0; X1 += 2*es;
            for (int n=es; n>=0; n-=(int)dil, X1-=2u*dil, X2+=2)
            {
                smr += *X2**X1 - *(X2+1)**(X1+1);
                smi += *X2**(X1+1) + *(X2+1)**X1;
            }
            *Y++ = smr; *Y++ = smi;
            X2 -= 2*(1+es/(int)dil);
            X1 += 2u*(dil-(size_t)es%dil);
            ++w; ss += stp; es += stp;
        }
        X2 += 2u*L2 - 2u; X1 += 2*ss;

        //Frames fully within sig
        while (es<(int)L1 && w<W)
        {
            smr = smi = 0.0;
            for (size_t l=0u; l<L2; ++l, X1+=2u*dil, X2-=2)
            {
                smr += *X2**X1 - *(X2+1)**(X1+1);
                smi += *X2**(X1+1) + *(X2+1)**X1;
            }
            *Y++ = smr; *Y++ = smi;
            X2 += 2u*L2; X1 += inc;
            ++w; es += stp;
        }
        ss = es - Nb + 1;

        //Frames overlapping end samp
        while (ss<(int)L1 && w<W)
        {
            size_t l = 0u; smr = smi = 0.0;
            for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=2u*dil, X2-=2)
            {
                smr += *X2**X1 - *(X2+1)**(X1+1);
                smi += *X2**(X1+1) + *(X2+1)**X1;
            }
            *Y++ = smr; *Y++ = smi;
            X2 += 2u*l; X1 += 2*((int)stp-(int)(dil*l));
            ++w; ss += stp;
        }

        //Frames after end samp
        while (w<W) { *Y++ = 0.0; *Y++ = 0.0; ++w; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L1, G = V/B;

        const int Nb = (int)(dil*(L2-1u));                  //full length of X2 including dil
        const int inc = 2*(int)K*((int)stp-(int)(dil*L2));  //fixed increment for X1 below
        size_t w;                                           //current frame
        int ss, es;                                         //current start-samp, end-samp
        double smr, smi;                                    //intermediate sums

        for (size_t g=0u; g<G; ++g, X1+=2u*B*(L1-1u), Y+=2u*B*(W-1u))
        {
            for (size_t b=0u; b<B; ++b, X1-=2u*K*L1-2u, Y-=2u*K*W-2u)
            {
                //Start frames
                w = 0u; es = c0 + (int)(dil*(L2/2u));

                //Frames before first samp
                while (es<0 && w<W) { *Y = *(Y+1) = 0.0; Y+=2u*K; ++w; es += stp; }
                ss = es - Nb + 1;

                //Frames overlapping first samp
                while (ss<0 && w<W)
                {
                    smr = smi = 0.0; X1 += 2*(int)K*es;
                    for (int n=es; n>=0; n-=(int)dil, X1-=2u*K*dil, X2+=2)
                    {
                        smr += *X2**X1 - *(X2+1)**(X1+1);
                        smi += *X2**(X1+1) + *(X2+1)**X1;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    X2 -= 2*(1+es/(int)dil);
                    X1 += 2u*K*(dil-(size_t)es%dil);
                    ++w; ss += stp; es += stp;
                }
                X2 += 2u*L2 - 2u; X1 += 2*ss*(int)K;

                //Frames fully within sig
                while (es<(int)L1 && w<W)
                {
                    smr = smi = 0.0;
                    for (size_t l=0u; l<L2; ++l, X1+=2u*K*dil, X2-=2)
                    {
                        smr += *X2**X1 - *(X2+1)**(X1+1);
                        smi += *X2**(X1+1) + *(X2+1)**X1;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    X2 += 2u*L2; X1 += inc;
                    ++w; es += stp;
                }
                ss = es - Nb + 1;

                //Frames overlapping end samp
                while (ss<(int)L1 && w<W)
                {
                    size_t l = 0u; smr = smi = 0.0;
                    for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=2u*K*dil, X2-=2)
                    {
                        smr += *X2**X1 - *(X2+1)**(X1+1);
                        smi += *X2**(X1+1) + *(X2+1)**X1;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    X2 += 2u*l; X1 += 2*(int)K*((int)stp-(int)(dil*l));
                    ++w; ss += stp;
                }
                X2 -= 2u*L2 - 2u; X1 -= 2u*K*((size_t)ss-L1);

                //Frames after end samp
                while (w<W) { *Y = *(Y+1) = 0.0; Y+=2u*K; ++w; }
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
