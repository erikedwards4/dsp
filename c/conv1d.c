//1D convolution by B of each vector in X along dim.
//Each vector in X has length Lx. B has length Lb.

//FIR filtering is similar, except FIR is causal and conv is non-causal.
//Note that some "convolution" functions actually do cross-correlation.
//For actual cross-correlation (no flip of B), see xcorr1d.

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int conv1d_s (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lb, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim);
int conv1d_d (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lb, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim);
int conv1d_c (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lb, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim);
int conv1d_z (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lb, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim);


int conv1d_s (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lb, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv1d_s: dim must be in [0 3]\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in conv1d_s: stp (step size) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in conv1d_s: dil (dilation) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in conv1d_s: N (total length of X) must be positive\n"); return 1; }
    if (Lx<1u) { fprintf(stderr,"error in conv1d_s: Lx (length of vecs in X) must be positive\n"); return 1; }
    if (c0>=(int)Lx) { fprintf(stderr,"error in conv1d_s: c0 (center samp of 1st frame) must be < Lx (length of vecs in X)\n"); return 1; }
    //if (W>N) { fprintf(stderr,"error in conv1d_s: W must be <= N (length X)\n"); return 1; }

    if (W==0u) {}
    else if (Lx==N)
    {
        const int Nb = (int)(dil*(Lb-1u)) + 1;          //full length of B including dil
        const int inc = (int)stp - (int)(dil*Lb);       //fixed increment for X below
        size_t w = 0u;                                  //current frame
        int ss, es = c0 + (int)(dil*(Lb/2u));           //current start-samp, end-samp
        float sm;                                       //intermediate sum

        //Frames before first samp
        while (es<0 && w<W) { *Y++ = 0.0f; ++w; es += stp; }
        ss = es - Nb + 1;

        //Frames overlapping first samp
        while (ss<0 && es<(int)Lx && w<W)
        {
            sm = 0.0f; X += es;
            for (int n=es; n>=0; n-=(int)dil, X-=dil, ++B) { sm = fmaf(*B,*X,sm); }
            *Y++ = sm;
            B -= 1 + es/(int)dil; X += dil - (size_t)es%dil;
            ++w; ss += stp; es += stp;
        }
        B += Lb - 1u;

        //In case of Lb>Lx
        es = (int)Lx - 1;
        while (ss<0 && w<W)
        {
            sm = 0.0f; X += es; B -= Lb - Lx;
            for (int n=es; n>=0; n-=(int)dil, X-=dil, ++B) { sm = fmaf(*B,*X,sm); }
            *Y++ = sm;
            B -= 1 + es/(int)dil; X += dil - (size_t)es%dil;
            ++w; ss += stp;
        }
        B += Lb - 1u; X += ss;

        //Frames fully within sig
        while (es<(int)Lx && w<W)
        {
            sm = 0.0f;
            for (size_t l=0u; l<Lb; ++l, X+=dil, --B) { sm = fmaf(*B,*X,sm); }
            *Y++ = sm;
            B += Lb; X += inc;
            ++w; es += stp;
        }
        ss = es - Nb + 1;

        //Frames overlapping end samp
        while (ss<(int)Lx && w<W)
        {
            size_t l = 0u; sm = 0.0f;
            for (int n=ss; n<(int)Lx; n+=dil, ++l, X+=dil, --B) { sm = fmaf(*B,*X,sm); }
            *Y++ = sm;
            B += l; X += (int)stp - (int)(dil*l);
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

        const int Nb = (int)(dil*(Lb-1u));                  //full length of B including dil
        const int inc = (int)K*((int)stp-(int)(dil*Lb));    //fixed increment for X below
        size_t w;                                           //current frame
        int ss, es;                                         //current start-samp, end-samp
        float sm;                                           //intermediate sum

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
                while (ss<0 && w<W)
                {
                    sm = 0.0f; X += (int)K*es;
                    for (int n=es; n>=0; n-=(int)dil, X-=K*dil, ++B) { sm = fmaf(*B,*X,sm); }
                    *Y = sm; Y += K;
                    B -= 1 + es/(int)dil; X += K*(dil-(size_t)es%dil);
                    ++w; ss += stp; es += stp;
                }
                B += Lb - 1u; X += ss*(int)K;

                //Frames fully within sig
                while (es<(int)Lx && w<W)
                {
                    sm = 0.0f;
                    for (size_t l=0u; l<Lb; ++l, X+=K*dil, --B) { sm = fmaf(*B,*X,sm); }
                    *Y = sm; Y += K;
                    B += Lb; X += inc;
                    ++w; es += stp;
                }
                ss = es - Nb + 1;

                //Frames overlapping end samp
                while (ss<(int)Lx && w<W)
                {
                    size_t l = 0u; sm = 0.0f;
                    for (int n=ss; n<(int)Lx; n+=dil, X+=K*dil, --B, ++l) { sm = fmaf(*B,*X,sm); }
                    *Y = sm; Y += K;
                    B += l; X += (int)K*((int)stp-(int)(dil*l));
                    ++w; ss += stp;
                }
                B -= Lb - 1u; X -= K*((size_t)ss-Lx);

                //Frames after end samp
                while (w<W) { *Y = 0.0f; Y+=K; ++w; }
            }
        }
    }

    return 0;
}


int conv1d_d (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lb, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv1d_d: dim must be in [0 3]\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in conv1d_d: stp (step size) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in conv1d_d: dil (dilation) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in conv1d_d: N (total length of X) must be positive\n"); return 1; }
    if (Lx<1u) { fprintf(stderr,"error in conv1d_d: Lx (length of vecs in X) must be positive\n"); return 1; }
    if (c0>=(int)Lx) { fprintf(stderr,"error in conv1d_d: c0 (center samp of 1st frame) must be < Lx (length of vecs in X)\n"); return 1; }

    if (W==0u) {}
    else if (Lx==N)
    {
        const int Nb = (int)(dil*(Lb-1u)) + 1;          //full length of B including dil
        const int inc = (int)stp - (int)(dil*Lb);       //fixed increment for X below
        size_t w = 0u;                                  //current frame
        int ss, es = c0 + (int)(dil*(Lb/2u));           //current start-samp, end-samp
        double sm;                                      //intermediate sum

        //Frames before first samp
        while (es<0 && w<W) { *Y++ = 0.0; ++w; es += stp; }
        ss = es - Nb + 1;

        //Frames overlapping first samp
        while (ss<0 && w<W)
        {
            sm = 0.0; X += es;
            for (int n=es; n>=0; n-=(int)dil, X-=dil, ++B) { sm = fma(*B,*X,sm); }
            *Y++ = sm;
            B -= 1 + es/(int)dil; X += dil - (size_t)es%dil;
            ++w; ss += stp; es += stp;
        }
        B += Lb - 1u; X += ss;

        //Frames fully within sig
        while (es<(int)Lx && w<W)
        {
            sm = 0.0;
            for (size_t l=0u; l<Lb; ++l, X+=dil, --B) { sm = fma(*B,*X,sm); }
            *Y++ = sm;
            B += Lb; X += inc;
            ++w; es += stp;
        }
        ss = es - Nb + 1;

        //Frames overlapping end samp
        while (ss<(int)Lx && w<W)
        {
            size_t l = 0u; sm = 0.0;
            for (int n=ss; n<(int)Lx; n+=dil, ++l, X+=dil, --B) { sm = fma(*B,*X,sm); }
            *Y++ = sm;
            B += l; X += (int)stp - (int)(dil*l);
            ++w; ss += stp;
        }

        //Frames after end samp
        while (w<W) { *Y++ = 0.0; ++w; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/BS;

        const int Nb = (int)(dil*(Lb-1u));                  //full length of B including dil
        const int inc = (int)K*((int)stp-(int)(dil*Lb));    //fixed increment for X below
        size_t w;                                           //current frame
        int ss, es;                                         //current start-samp, end-samp
        double sm;                                          //intermediate sum

        for (size_t g=0u; g<G; ++g, X+=BS*(Lx-1u), Y+=BS*(W-1u))
        {
            for (size_t b=0u; b<BS; ++b, X-=K*Lx-1u, Y-=K*W-1u)
            {
                //Start frames
                w = 0u; es = c0 + (int)(dil*(Lb/2u));

                //Frames before first samp
                while (es<0 && w<W) { *Y = 0.0; Y+=K; ++w; es += stp; }
                ss = es - Nb + 1;

                //Frames overlapping first samp
                while (ss<0 && w<W)
                {
                    sm = 0.0; X += (int)K*es;
                    for (int n=es; n>=0; n-=(int)dil, X-=K*dil, ++B) { sm = fma(*B,*X,sm); }
                    *Y = sm; Y += K;
                    B -= 1 + es/(int)dil; X += K*(dil-(size_t)es%dil);
                    ++w; ss += stp; es += stp;
                }
                B += Lb - 1u; X += ss*(int)K;

                //Frames fully within sig
                while (es<(int)Lx && w<W)
                {
                    sm = 0.0;
                    for (size_t l=0u; l<Lb; ++l, X+=K*dil, --B) { sm = fma(*B,*X,sm); }
                    *Y = sm; Y += K;
                    B += Lb; X += inc;
                    ++w; es += stp;
                }
                ss = es - Nb + 1;

                //Frames overlapping end samp
                while (ss<(int)Lx && w<W)
                {
                    size_t l = 0u; sm = 0.0;
                    for (int n=ss; n<(int)Lx; n+=dil, X+=K*dil, --B, ++l) { sm = fma(*B,*X,sm); }
                    *Y = sm; Y += K;
                    B += l; X += (int)K*((int)stp-(int)(dil*l));
                    ++w; ss += stp;
                }
                B -= Lb - 1u; X -= K*((size_t)ss-Lx);

                //Frames after end samp
                while (w<W) { *Y = 0.0; Y+=K; ++w; }
            }
        }
    }

    return 0;
}


int conv1d_c (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lb, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv1d_c: dim must be in [0 3]\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in conv1d_c: stp (step size) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in conv1d_c: dil (dilation) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in conv1d_c: N (total length of X) must be positive\n"); return 1; }
    if (Lx<1u) { fprintf(stderr,"error in conv1d_c: Lx (length of vecs in X) must be positive\n"); return 1; }
    if (c0>=(int)Lx) { fprintf(stderr,"error in conv1d_c: c0 (center samp of 1st frame) must be < Lx (length of vecs in X)\n"); return 1; }

    if (W==0u) {}
    else if (Lx==N)
    {
        const int Nb = (int)(dil*(Lb-1u)) + 1;          //full length of B including dil
        const int inc = 2*((int)stp-(int)(dil*Lb));     //fixed increment for X below
        size_t w = 0u;                                  //current frame
        int ss, es = c0 + (int)(dil*(Lb/2u));           //current start-samp, end-samp
        float smr, smi;                                 //intermediate sums

        //Frames before first samp
        while (es<0 && w<W) { *Y++ = 0.0f; *Y++ = 0.0f; ++w; es += stp; }
        ss = es - Nb + 1;

        //Frames overlapping first samp
        while (ss<0 && w<W)
        {
            smr = smi = 0.0f; X += 2*es;
            for (int n=es; n>=0; n-=(int)dil, X-=2u*dil, B+=2)
            {
                smr += *B**X - *(B+1)**(X+1);
                smi += *B**(X+1) + *(B+1)**X;
            }
            *Y++ = smr; *Y++ = smi;
            B -= 2*(1+es/(int)dil);
            X += 2u*(dil-(size_t)es%dil);
            ++w; ss += stp; es += stp;
        }
        B += 2u*Lb - 2u; X += 2*ss;

        //Frames fully within sig
        while (es<(int)Lx && w<W)
        {
            smr = smi = 0.0f;
            for (size_t l=0u; l<Lb; ++l, X+=2u*dil, B-=2)
            {
                smr += *B**X - *(B+1)**(X+1);
                smi += *B**(X+1) + *(B+1)**X;
            }
            *Y++ = smr; *Y++ = smi;
            B += 2u*Lb; X += inc;
            ++w; es += stp;
        }
        ss = es - Nb + 1;

        //Frames overlapping end samp
        while (ss<(int)Lx && w<W)
        {
            size_t l = 0u; smr = smi = 0.0f;
            for (int n=ss; n<(int)Lx; n+=dil, ++l, X+=2u*dil, B-=2)
            {
                smr += *B**X - *(B+1)**(X+1);
                smi += *B**(X+1) + *(B+1)**X;
            }
            *Y++ = smr; *Y++ = smi;
            B += 2u*l; X += 2*((int)stp-(int)(dil*l));
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

        const int Nb = (int)(dil*(Lb-1u));                  //full length of B including dil
        const int inc = 2*(int)K*((int)stp-(int)(dil*Lb));  //fixed increment for X below
        size_t w;                                           //current frame
        int ss, es;                                         //current start-samp, end-samp
        float smr, smi;                                     //intermediate sums

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
                    for (int n=es; n>=0; n-=(int)dil, X-=2u*K*dil, B+=2)
                    {
                        smr += *B**X - *(B+1)**(X+1);
                        smi += *B**(X+1) + *(B+1)**X;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    B -= 2*(1+es/(int)dil);
                    X += 2u*K*(dil-(size_t)es%dil);
                    ++w; ss += stp; es += stp;
                }
                B += 2u*Lb - 2u; X += 2*ss*(int)K;

                //Frames fully within sig
                while (es<(int)Lx && w<W)
                {
                    smr = smi = 0.0f;
                    for (size_t l=0u; l<Lb; ++l, X+=2u*K*dil, B-=2)
                    {
                        smr += *B**X - *(B+1)**(X+1);
                        smi += *B**(X+1) + *(B+1)**X;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    B += 2u*Lb; X += inc;
                    ++w; es += stp;
                }
                ss = es - Nb + 1;

                //Frames overlapping end samp
                while (ss<(int)Lx && w<W)
                {
                    size_t l = 0u; smr = smi = 0.0f;
                    for (int n=ss; n<(int)Lx; n+=dil, ++l, X+=2u*K*dil, B-=2)
                    {
                        smr += *B**X - *(B+1)**(X+1);
                        smi += *B**(X+1) + *(B+1)**X;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    B += 2u*l; X += 2*(int)K*((int)stp-(int)(dil*l));
                    ++w; ss += stp;
                }
                B -= 2u*Lb - 2u; X -= 2u*K*((size_t)ss-Lx);

                //Frames after end samp
                while (w<W) { *Y = *(Y+1) = 0.0f; Y+=2u*K; ++w; }
            }
        }
    }

    return 0;
}


int conv1d_z (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lb, const size_t W, const int c0, const size_t stp, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv1d_z: dim must be in [0 3]\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in conv1d_z: stp (step size) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in conv1d_z: dil (dilation) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in conv1d_z: N (total length of X) must be positive\n"); return 1; }
    if (Lx<1u) { fprintf(stderr,"error in conv1d_z: Lx (length of vecs in X) must be positive\n"); return 1; }
    if (c0>=(int)Lx) { fprintf(stderr,"error in conv1d_z: c0 (center samp of 1st frame) must be < Lx (length of vecs in X)\n"); return 1; }

    if (W==0u) {}
    else if (Lx==N)
    {
        const int Nb = (int)(dil*(Lb-1u)) + 1;          //full length of B including dil
        const int inc = 2*((int)stp-(int)(dil*Lb));     //fixed increment for X below
        size_t w = 0u;                                  //current frame
        int ss, es = c0 + (int)(dil*(Lb/2u));           //current start-samp, end-samp
        double smr, smi;                                //intermediate sums

        //Frames before first samp
        while (es<0 && w<W) { *Y++ = 0.0; *Y++ = 0.0; ++w; es += stp; }
        ss = es - Nb + 1;

        //Frames overlapping first samp
        while (ss<0 && w<W)
        {
            smr = smi = 0.0; X += 2*es;
            for (int n=es; n>=0; n-=(int)dil, X-=2u*dil, B+=2)
            {
                smr += *B**X - *(B+1)**(X+1);
                smi += *B**(X+1) + *(B+1)**X;
            }
            *Y++ = smr; *Y++ = smi;
            B -= 2*(1+es/(int)dil);
            X += 2u*(dil-(size_t)es%dil);
            ++w; ss += stp; es += stp;
        }
        B += 2u*Lb - 2u; X += 2*ss;

        //Frames fully within sig
        while (es<(int)Lx && w<W)
        {
            smr = smi = 0.0;
            for (size_t l=0u; l<Lb; ++l, X+=2u*dil, B-=2)
            {
                smr += *B**X - *(B+1)**(X+1);
                smi += *B**(X+1) + *(B+1)**X;
            }
            *Y++ = smr; *Y++ = smi;
            B += 2u*Lb; X += inc;
            ++w; es += stp;
        }
        ss = es - Nb + 1;

        //Frames overlapping end samp
        while (ss<(int)Lx && w<W)
        {
            size_t l = 0u; smr = smi = 0.0;
            for (int n=ss; n<(int)Lx; n+=dil, ++l, X+=2u*dil, B-=2)
            {
                smr += *B**X - *(B+1)**(X+1);
                smi += *B**(X+1) + *(B+1)**X;
            }
            *Y++ = smr; *Y++ = smi;
            B += 2u*l; X += 2*((int)stp-(int)(dil*l));
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

        const int Nb = (int)(dil*(Lb-1u));                  //full length of B including dil
        const int inc = 2*(int)K*((int)stp-(int)(dil*Lb));  //fixed increment for X below
        size_t w;                                           //current frame
        int ss, es;                                         //current start-samp, end-samp
        double smr, smi;                                    //intermediate sums

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
                    for (int n=es; n>=0; n-=(int)dil, X-=2u*K*dil, B+=2)
                    {
                        smr += *B**X - *(B+1)**(X+1);
                        smi += *B**(X+1) + *(B+1)**X;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    B -= 2*(1+es/(int)dil);
                    X += 2u*K*(dil-(size_t)es%dil);
                    ++w; ss += stp; es += stp;
                }
                B += 2u*Lb - 2u; X += 2*ss*(int)K;

                //Frames fully within sig
                while (es<(int)Lx && w<W)
                {
                    smr = smi = 0.0;
                    for (size_t l=0u; l<Lb; ++l, X+=2u*K*dil, B-=2)
                    {
                        smr += *B**X - *(B+1)**(X+1);
                        smi += *B**(X+1) + *(B+1)**X;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    B += 2u*Lb; X += inc;
                    ++w; es += stp;
                }
                ss = es - Nb + 1;

                //Frames overlapping end samp
                while (ss<(int)Lx && w<W)
                {
                    size_t l = 0u; smr = smi = 0.0;
                    for (int n=ss; n<(int)Lx; n+=dil, ++l, X+=2u*K*dil, B-=2)
                    {
                        smr += *B**X - *(B+1)**(X+1);
                        smi += *B**(X+1) + *(B+1)**X;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    B += 2u*l; X += 2*(int)K*((int)stp-(int)(dil*l));
                    ++w; ss += stp;
                }
                B -= 2u*Lb - 2u; X -= 2u*K*((size_t)ss-Lx);

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
