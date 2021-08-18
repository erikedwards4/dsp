//1D cross-correlation by B of each vector in X along dim.
//Each vector in X has length Lx. B has length Lb.

//This is identical to conv1d, except B is NOT flipped.
//Note that some "convolution" functions out there actually do cross-correlation.
//Note that this does NOT attempt to emmulate Octave's xcorr (but is similar in principle).

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int xcorr1d_s (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lb, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim);
int xcorr1d_d (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lb, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim);
int xcorr1d_c (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lb, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim);
int xcorr1d_z (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lb, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim);


int xcorr1d_s (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lb, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in xcorr1d_s: dim must be in [0 3]\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in xcorr1d_s: stp (step size) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in xcorr1d_s: dil (dilation) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in xcorr1d_s: N (total length of X) must be positive\n"); return 1; }
    if (Lx<1u) { fprintf(stderr,"error in xcorr1d_s: Lx (length of vecs in X) must be positive\n"); return 1; }
    if (c0>=(int)Lx) { fprintf(stderr,"error in xcorr1d_s: c0 (center samp of 1st frame) must be < Lx (length of vecs in X)\n"); return 1; }

    if (W==0u) {}
    else if (Lx==N)
    {
        const int Nb = (int)(dil*(Lb-1u)) + 1;          //full length of B including dil
        const int inc = (int)stp - (int)(dil*Lb);       //fixed increment for X below
        size_t w = 0u;                                  //current frame
        int ss, es = c0 + (int)(dil*(Lb/2u));           //current start-samp, end-samp
        float sm;                                       //intermediate sum

        //Don't flip B
        B += Lb - 1u;

        //Frames before first samp
        while (es<0 && w<W) { *Y++ = 0.0f; ++w; es += stp; }
        ss = es - Nb + 1;

        //Frames overlapping first samp
        while (ss<0 && es<(int)Lx && w<W)
        {
            sm = 0.0f; X += es;
            for (int n=es; n>=0; n-=(int)dil, X-=dil, --B) { sm = fmaf(*B,*X,sm); }
            *Y++ = sm;
            B += 1 + es/(int)dil; X += dil - (size_t)es%dil;
            ++w; ss += stp; es += stp;
        }
        B -= Lb - 1u;

        //In case Lb>Lx
        while (ss<0 && w<W)
        {
            int l = ss + (int)dil*((1-ss)/(int)dil);
            B += (l-ss)/(int)dil; X += l;
            sm = 0.0f;
            while (l<(int)Lx) { sm = fmaf(*B,*X,sm); X+=dil; l+=(int)dil; ++B; }
            *Y++ = sm;
            B -= (l-ss)/(int)dil; X -= l;
            ++w; ss += stp;
        }
        es = ss + Nb - 1;

        //Frames fully within sig
        while (es<(int)Lx && w<W)
        {
            sm = 0.0f;
            for (size_t l=0u; l<Lb; ++l, X+=dil, ++B) { sm = fmaf(*B,*X,sm); }
            *Y++ = sm;
            B -= Lb; X += inc;
            ++w; es += stp;
        }
        ss = es - Nb + 1;

        //Frames overlapping end samp
        while (ss<(int)Lx && w<W)
        {
            size_t l = 0u; sm = 0.0f;
            for (int n=ss; n<(int)Lx; n+=dil, ++l, X+=dil, ++B) { sm = fmaf(*B,*X,sm); }
            *Y++ = sm;
            B -= l; X += (int)stp - (int)(dil*l);
            ++w; ss += stp;
        }

        //Frames after end samp
        while (w<W) { *Y++ = 0.0f; ++w; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/BS;

        const int Nb = (int)(dil*(Lb-1u)) + 1;              //full length of B including dil
        const int inc = (int)K*((int)stp-(int)(dil*Lb));    //fixed increment for X below
        size_t w;                                           //current frame
        int ss, es;                                         //current start-samp, end-samp
        float sm;                                           //intermediate sum

        //Don't flip B
        B += Lb - 1u;

        for (size_t g=0u; g<G; ++g, X+=BS*(Lx-1u), Y+=BS*(W-1u))
        {
            for (size_t b=0u; b<BS; ++b, X-=K*Lx-1u, Y-=K*W-1u)
            {
                //Start frames
                w = 0u; es = c0 + (int)(dil*(Lb/2u));

                //Frames before first samp
                while (es<0 && w<W) { *Y = 0.0f; Y+=K; ++w; es += stp; }
                ss = es - Nb + 1;

                //Frames overlapping first samp
                while (ss<0 && es<(int)Lx && w<W)
                {
                    sm = 0.0f; X += (int)K*es;
                    for (int n=es; n>=0; n-=(int)dil, X-=K*dil, --B) { sm = fmaf(*B,*X,sm); }
                    *Y = sm; Y += K;
                    B += 1 + es/(int)dil; X += K*(dil-(size_t)es%dil);
                    ++w; ss += stp; es += stp;
                }
                B -= Lb - 1u;

                //In case Lb>Lx
                while (ss<0 && w<W)
                {
                    int l = ss + (int)dil*((1-ss)/(int)dil);
                    B += (l-ss)/(int)dil; X += l*(int)K;
                    sm = 0.0f;
                    while (l<(int)Lx) { sm = fmaf(*B,*X,sm); X+=dil*K; l+=(int)dil; ++B; }
                    *Y = sm; Y += K;
                    B -= (l-ss)/(int)dil; X -= l*(int)K;
                    ++w; ss += stp;
                }
                es = ss + Nb - 1;

                //Frames fully within sig
                while (es<(int)Lx && w<W)
                {
                    sm = 0.0f;
                    for (size_t l=0u; l<Lb; ++l, X+=K*dil, ++B) { sm = fmaf(*B,*X,sm); }
                    *Y = sm; Y += K;
                    B -= Lb; X += inc;
                    ++w; es += stp;
                }
                ss = es - Nb + 1;

                //Frames overlapping end samp
                while (ss<(int)Lx && w<W)
                {
                    size_t l = 0u; sm = 0.0f;
                    for (int n=ss; n<(int)Lx; n+=dil, X+=K*dil, ++B, ++l) { sm = fmaf(*B,*X,sm); }
                    *Y = sm; Y += K;
                    B -= l; X += (int)K*((int)stp-(int)(dil*l));
                    ++w; ss += stp;
                }
                B += Lb - 1u; X -= K*((size_t)ss-Lx);

                //Frames after end samp
                while (w<W) { *Y = 0.0f; Y+=K; ++w; }
            }
        }
    }

    return 0;
}


int xcorr1d_d (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lb, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in xcorr1d_d: dim must be in [0 3]\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in xcorr1d_d: stp (step size) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in xcorr1d_d: dil (dilation) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in xcorr1d_d: N (total length of X) must be positive\n"); return 1; }
    if (Lx<1u) { fprintf(stderr,"error in xcorr1d_d: Lx (length of vecs in X) must be positive\n"); return 1; }
    if (c0>=(int)Lx) { fprintf(stderr,"error in xcorr1d_d: c0 (start samp of 1st frame) must be < Lx (length of vecs in X)\n"); return 1; }

    if (W==0u) {}
    else if (Lx==N)
    {
        const int Nb = (int)(dil*(Lb-1u)) + 1;      //full length of B including dil
        const int inc = (int)stp - (int)(dil*Lb);   //fixed increment for X below
        size_t w = 0u;                              //current frame
        int ss, es = c0 + Nb/2;                     //current start-samp, end-samp
        double sm;                                  //intermediate sum

        //Don't flip B
        B += Lb - 1u;

        //Frames before first samp
        while (es<0 && w<W) { *Y++ = 0.0; es += stp; ++w; }
        ss = es - Nb + 1;

        //Frames overlapping first samp
        while (ss<0 && es<(int)Lx && w<W)
        {
            sm = 0.0; X += es;
            for (int n=es; n>=0; n-=(int)dil, X-=dil, --B) { sm = fma(*B,*X,sm); }
            *Y++ = sm;
            B += 1 + es/(int)dil; X += dil - (size_t)es%dil;
            ss += stp; es += stp; ++w;
        }
        B -= Lb - 1u;

        //In case Lb>Lx
        while (ss<0 && w<W)
        {
            int l = ss + (int)dil*((1-ss)/(int)dil);
            B += (l-ss)/(int)dil; X += l;
            sm = 0.0;
            while (l<(int)Lx) { sm = fma(*B,*X,sm); X+=dil; l+=(int)dil; ++B; }
            *Y++ = sm;
            B -= (l-ss)/(int)dil; X -= l;
            ss += stp; ++w;
        }
        es = ss + Nb - 1;

        //Frames fully within sig
        while (es<(int)Lx && w<W)
        {
            sm = 0.0;
            for (size_t l=0u; l<Lb; ++l, X+=dil, ++B) { sm = fma(*B,*X,sm); }
            *Y++ = sm;
            B -= Lb; X += inc;
            es += stp; ++w;
        }
        ss = es - Nb + 1;

        //Frames overlapping end samp
        while (ss<(int)Lx && w<W)
        {
            size_t l = 0u; sm = 0.0;
            for (int n=ss; n<(int)Lx; n+=dil, ++l, X+=dil, ++B) { sm = fma(*B,*X,sm); }
            *Y++ = sm;
            B -= l; X += (int)stp - (int)(dil*l);
            ss += stp; ++w;
        }

        //Frames after end samp
        while (w<W) { *Y++ = 0.0; ++w; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/BS;
        fprintf(stderr,"K=%lu, BS=%lu, V=%lu, G=%lu \n",K,BS,V,G);

        const int Nb = (int)(dil*(Lb-1u)) + 1;              //full length of B including dil
        const int inc = (int)K*((int)stp-(int)(dil*Lb));    //fixed increment for X below
        size_t w;                                           //current frame
        int ss, es = c0 + Nb/2;                             //current start-samp, end-samp
        double sm;                                          //intermediate sum

        //Don't flip B
        B += Lb - 1u; int cnt = 0;

        for (size_t g=0u; g<G; ++g, X+=BS*(Lx-1u), Y+=BS*(W-1u))
        {
            for (size_t bs=0u; bs<BS; ++bs, X-=K*Lx-1u, Y-=K*W-1u)
            {
                //Start frames
                ss = c0; es = ss + Nb - 1; w = 0u;

                //Frames before first samp
                while (es<0 && w<W) { *Y = 0.0; Y+=K; es+=stp; ++w; }
                ss = es - Nb + 1;

                //Frames overlapping first samp
                while (ss<0 && es<(int)Lx && w<W)
                {
                    sm = 0.0; X += (int)K*es;
                    for (int n=es; n>=0; n-=(int)dil, X-=K*dil, --B) { sm = fma(*B,*X,sm); }
                    *Y = sm; Y += K;
                    B += 1 + es/(int)dil; X += K*(dil-(size_t)es%dil);
                    ss += stp; es += stp; ++w;
                }
                B -= Lb - 1u; fprintf(stderr,"Y=%g\n",(double)(*(Y-K)));

                //In case Lb>Lx
                while (ss<0 && w<W)
                {
                    // int lb = 0, nb = 0, lx; sm = 0.0;
                    // while (lb<Lb)
                    // {
                    //     if (nb>-ss && nb<(int)Lx)
                    //     {
                    //         sm = fma(*B,*X,sm);
                    //     }
                    //     ++lb; nb+=(int)dil;
                    // }
                    int b = (1-ss)/(int)dil, l = ss + (int)dil*b;
                    //int l = ss + (int)dil*((1-ss)/(int)dil);
                    B += b; X += l*(int)K;
                    sm = 0.0;
                    while (l<(int)Lx) { sm = fma(*B,*X,sm); X+=dil*K; l+=(int)dil; ++B; ++b; }
                    *Y = sm; Y += K;
                    B -= b; X -= l*(int)K;
                    // B -= (l-ss)/(int)dil; X -= l*(int)K;
                    ss += stp; ++w;
                }
                es = ss + Nb - 1;

                //Frames fully within sig
                while (es<(int)Lx && w<W)
                {
                    sm = 0.0;
                    for (size_t l=0u; l<Lb; ++l, X+=K*dil, ++B) { sm = fma(*B,*X,sm); }
                    *Y = sm; Y += K;
                    B -= Lb; X += inc;
                    es += stp; ++w;
                }
                ss = es - Nb + 1;

                //Frames overlapping end samp
                while (ss<(int)Lx && w<W)
                {
                    size_t l = 0u; sm = 0.0;
                    for (int n=ss; n<(int)Lx; n+=dil, X+=K*dil, ++B, ++l) { sm = fma(*B,*X,sm); }
                    *Y = sm; Y += K;
                    B -= l; X += (int)K*((int)stp-(int)(dil*l));
                    ss += stp; ++w;
                }
                B += Lb - 1u; X -= K*dil*((size_t)ss-Lx);

                //Frames after end samp
                while (w<W) { *Y = 0.0; Y+=K; ++w; }
            }
        }
    }

    return 0;
}


int xcorr1d_c (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lb, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in xcorr1d_c: dim must be in [0 3]\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in xcorr1d_c: stp (step size) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in xcorr1d_c: dil (dilation) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in xcorr1d_c: N (total length of X) must be positive\n"); return 1; }
    if (Lx<1u) { fprintf(stderr,"error in xcorr1d_c: Lx (length of vecs in X) must be positive\n"); return 1; }
    if (c0>=(int)Lx) { fprintf(stderr,"error in xcorr1d_c: c0 (center samp of 1st frame) must be < Lx (length of vecs in X)\n"); return 1; }

    if (W==0u) {}
    else if (Lx==N)
    {
        const int Nb = (int)(dil*(Lb-1u)) + 1;          //full length of B including dil
        const int inc = 2*((int)stp-(int)(dil*Lb));     //fixed increment for X below
        size_t w = 0u;                                  //current frame
        int ss, es = c0 + (int)(dil*(Lb/2u));           //current start-samp, end-samp
        float smr, smi;                                 //intermediate sums

        //Don't flip B
        B += 2u*Lb - 2u;

        //Frames before first samp
        while (es<0 && w<W) { *Y++ = 0.0f; *Y++ = 0.0f; ++w; es += stp; }
        ss = es - Nb + 1;

        //Frames overlapping first samp
        while (ss<0 && w<W)
        {
            smr = smi = 0.0f; X += 2*es;
            for (int n=es; n>=0; n-=(int)dil, X-=2u*dil, B-=2)
            {
                smr += *B**X - *(B+1)**(X+1);
                smi += *B**(X+1) + *(B+1)**X;
            }
            *Y++ = smr; *Y++ = smi;
            B += 2*(1+es/(int)dil);
            X += 2u*(dil-(size_t)es%dil);
            ++w; ss += stp; es += stp;
        }
        B -= 2u*Lb - 2u; X += 2*ss;

        //Frames fully within sig
        while (es<(int)Lx && w<W)
        {
            smr = smi = 0.0f;
            for (size_t l=0u; l<Lb; ++l, X+=2u*dil, B+=2)
            {
                smr += *B**X - *(B+1)**(X+1);
                smi += *B**(X+1) + *(B+1)**X;
            }
            *Y++ = smr; *Y++ = smi;
            B -= 2u*Lb; X += inc;
            ++w; es += stp;
        }
        ss = es - Nb + 1;

        //Frames overlapping end samp
        while (ss<(int)Lx && w<W)
        {
            size_t l = 0u; smr = smi = 0.0f;
            for (int n=ss; n<(int)Lx; n+=dil, ++l, X+=2u*dil, B+=2)
            {
                smr += *B**X - *(B+1)**(X+1);
                smi += *B**(X+1) + *(B+1)**X;
            }
            *Y++ = smr; *Y++ = smi;
            B -= 2u*l; X += 2*((int)stp-(int)(dil*l));
            ++w; ss += stp;
        }

        //Frames after end samp
        while (w<W) { *Y++ = 0.0f; *Y++ = 0.0f; ++w; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/BS;

        const int Nb = (int)(dil*(Lb-1u)) + 1;              //full length of B including dil
        const int inc = 2*(int)K*((int)stp-(int)(dil*Lb));  //fixed increment for X below
        size_t w;                                           //current frame
        int ss, es;                                         //current start-samp, end-samp
        float smr, smi;                                     //intermediate sums

        //Don't flip B
        B += 2u*Lb - 2u;

        for (size_t g=0u; g<G; ++g, X+=2u*BS*(Lx-1u), Y+=2u*BS*(W-1u))
        {
            for (size_t b=0u; b<BS; ++b, X-=2u*K*Lx-2u, Y-=2u*K*W-2u)
            {
                //Start frames
                w = 0u; es = c0 + (int)(dil*(Lb/2u));

                //Frames before first samp
                while (es<0 && w<W) { *Y = *(Y+1) = 0.0f; Y+=2u*K; ++w; es += stp; }
                ss = es - Nb + 1;

                //Frames overlapping first samp
                while (ss<0 && w<W)
                {
                    smr = smi = 0.0f; X += 2*(int)K*es;
                    for (int n=es; n>=0; n-=(int)dil, X-=2u*K*dil, B-=2)
                    {
                        smr += *B**X - *(B+1)**(X+1);
                        smi += *B**(X+1) + *(B+1)**X;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    B += 2*(1+es/(int)dil);
                    X += 2u*K*(dil-(size_t)es%dil);
                    ++w; ss += stp; es += stp;
                }
                B -= 2u*Lb - 2u; X += 2*ss*(int)K;

                //Frames fully within sig
                while (es<(int)Lx && w<W)
                {
                    smr = smi = 0.0f;
                    for (size_t l=0u; l<Lb; ++l, X+=2u*K*dil, B+=2)
                    {
                        smr += *B**X - *(B+1)**(X+1);
                        smi += *B**(X+1) + *(B+1)**X;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    B -= 2u*Lb; X += inc;
                    ++w; es += stp;
                }
                ss = es - Nb + 1;

                //Frames overlapping end samp
                while (ss<(int)Lx && w<W)
                {
                    size_t l = 0u; smr = smi = 0.0f;
                    for (int n=ss; n<(int)Lx; n+=dil, ++l, X+=2u*K*dil, B+=2)
                    {
                        smr += *B**X - *(B+1)**(X+1);
                        smi += *B**(X+1) + *(B+1)**X;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    B -= 2u*l; X += 2*(int)K*((int)stp-(int)(dil*l));
                    ++w; ss += stp;
                }
                B += 2u*Lb - 2u; X -= 2u*K*((size_t)ss-Lx);

                //Frames after end samp
                while (w<W) { *Y = *(Y+1) = 0.0f; Y+=2u*K; ++w; }
            }
        }
    }

    return 0;
}


int xcorr1d_z (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lb, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in xcorr1d_z: dim must be in [0 3]\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in xcorr1d_z: stp (step size) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in xcorr1d_z: dil (dilation) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in xcorr1d_z: N (total length of X) must be positive\n"); return 1; }
    if (Lx<1u) { fprintf(stderr,"error in xcorr1d_z: Lx (length of vecs in X) must be positive\n"); return 1; }
    if (c0>=(int)Lx) { fprintf(stderr,"error in xcorr1d_z: c0 (center samp of 1st frame) must be < Lx (length of vecs in X)\n"); return 1; }

    if (W==0u) {}
    else if (Lx==N)
    {
        const int Nb = (int)(dil*(Lb-1u)) + 1;          //full length of B including dil
        const int inc = 2*((int)stp-(int)(dil*Lb));     //fixed increment for X below
        size_t w = 0u;                                  //current frame
        int ss, es = c0 + (int)(dil*(Lb/2u));           //current start-samp, end-samp
        double smr, smi;                                //intermediate sums

        //Don't flip B
        B += 2u*Lb - 2u;

        //Frames before first samp
        while (es<0 && w<W) { *Y++ = 0.0; *Y++ = 0.0; ++w; es += stp; }
        ss = es - Nb + 1;

        //Frames overlapping first samp
        while (ss<0 && w<W)
        {
            smr = smi = 0.0; X += 2*es;
            for (int n=es; n>=0; n-=(int)dil, X-=2u*dil, B-=2)
            {
                smr += *B**X - *(B+1)**(X+1);
                smi += *B**(X+1) + *(B+1)**X;
            }
            *Y++ = smr; *Y++ = smi;
            B += 2*(1+es/(int)dil);
            X += 2u*(dil-(size_t)es%dil);
            ++w; ss += stp; es += stp;
        }
        B -= 2u*Lb - 2u; X += 2*ss;

        //Frames fully within sig
        while (es<(int)Lx && w<W)
        {
            smr = smi = 0.0;
            for (size_t l=0u; l<Lb; ++l, X+=2u*dil, B+=2)
            {
                smr += *B**X - *(B+1)**(X+1);
                smi += *B**(X+1) + *(B+1)**X;
            }
            *Y++ = smr; *Y++ = smi;
            B -= 2u*Lb; X += inc;
            ++w; es += stp;
        }
        ss = es - Nb + 1;

        //Frames overlapping end samp
        while (ss<(int)Lx && w<W)
        {
            size_t l = 0u; smr = smi = 0.0;
            for (int n=ss; n<(int)Lx; n+=dil, ++l, X+=2u*dil, B+=2)
            {
                smr += *B**X - *(B+1)**(X+1);
                smi += *B**(X+1) + *(B+1)**X;
            }
            *Y++ = smr; *Y++ = smi;
            B -= 2u*l; X += 2*((int)stp-(int)(dil*l));
            ++w; ss += stp;
        }

        //Frames after end samp
        while (w<W) { *Y++ = 0.0; *Y++ = 0.0; ++w; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/BS;

        const int Nb = (int)(dil*(Lb-1u)) + 1;              //full length of B including dil
        const int inc = 2*(int)K*((int)stp-(int)(dil*Lb));  //fixed increment for X below
        size_t w;                                           //current frame
        int ss, es;                                         //current start-samp, end-samp
        double smr, smi;                                    //intermediate sums

        //Don't flip B
        B += 2u*Lb - 2u;

        for (size_t g=0u; g<G; ++g, X+=2u*BS*(Lx-1u), Y+=2u*BS*(W-1u))
        {
            for (size_t b=0u; b<BS; ++b, X-=2u*K*Lx-2u, Y-=2u*K*W-2u)
            {
                //Start frames
                w = 0u; es = c0 + (int)(dil*(Lb/2u));

                //Frames before first samp
                while (es<0 && w<W) { *Y = *(Y+1) = 0.0; Y+=2u*K; ++w; es += stp; }
                ss = es - Nb + 1;

                //Frames overlapping first samp
                while (ss<0 && w<W)
                {
                    smr = smi = 0.0; X += 2*(int)K*es;
                    for (int n=es; n>=0; n-=(int)dil, X-=2u*K*dil, B-=2)
                    {
                        smr += *B**X - *(B+1)**(X+1);
                        smi += *B**(X+1) + *(B+1)**X;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    B += 2*(1+es/(int)dil);
                    X += 2u*K*(dil-(size_t)es%dil);
                    ++w; ss += stp; es += stp;
                }
                B -= 2u*Lb - 2u; X += 2*ss*(int)K;

                //Frames fully within sig
                while (es<(int)Lx && w<W)
                {
                    smr = smi = 0.0;
                    for (size_t l=0u; l<Lb; ++l, X+=2u*K*dil, B+=2)
                    {
                        smr += *B**X - *(B+1)**(X+1);
                        smi += *B**(X+1) + *(B+1)**X;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    B -= 2u*Lb; X += inc;
                    ++w; es += stp;
                }
                ss = es - Nb + 1;

                //Frames overlapping end samp
                while (ss<(int)Lx && w<W)
                {
                    size_t l = 0u; smr = smi = 0.0;
                    for (int n=ss; n<(int)Lx; n+=dil, ++l, X+=2u*K*dil, B+=2)
                    {
                        smr += *B**X - *(B+1)**(X+1);
                        smi += *B**(X+1) + *(B+1)**X;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    B -= 2u*l; X += 2*(int)K*((int)stp-(int)(dil*l));
                    ++w; ss += stp;
                }
                B += 2u*Lb - 2u; X -= 2u*K*((size_t)ss-Lx);

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
