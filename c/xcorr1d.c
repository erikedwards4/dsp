//1D cross-correlation of each vector in X1 by X2.
//Each vector in X1 has length L1. X2 has length L2 (kernel_size).

//This is identical to conv1d, except X2 is NOT flipped.
//Note that some "convolution" functions out there actually do cross-correlation.
//Note that this does NOT attempt to emmulate Octave's xcorr (but is similar in principle).

//This emulates PyTorch Conv1d, but limited to C_out=C_in, padding_mode='zeros', groups=1, bias=false.
//That is, X2 is the convolving kernel with length kernel_size.
//Each vector in X1 has length L_in, and each vector in Y has length L_out, set by:
//L_out = floor[1 + (L_in + 2*pad - dil*(L2-1) - 1)/stride].

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int xcorr1d_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t pad, const size_t str, const size_t dil, const size_t dim);
int xcorr1d_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t pad, const size_t str, const size_t dil, const size_t dim);
int xcorr1d_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t pad, const size_t str, const size_t dil, const size_t dim);
int xcorr1d_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t pad, const size_t str, const size_t dil, const size_t dim);


int xcorr1d_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t pad, const size_t str, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in xcorr1d_s: dim must be in [0 3]\n"); return 1; }
    if (str<1u) { fprintf(stderr,"error in xcorr1d_s: str (stride) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in xcorr1d_s: dil (dilation) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in xcorr1d_s: L2 (length of X2) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in xcorr1d_s: N (total length of X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in xcorr1d_s: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if (L1+2u*pad<=dil*(L2-1u)) { fprintf(stderr,"error in xcorr1d_s: L1+2*pad must be > dil*(L2-1)\n"); return 1; }

    const int N1 = (int)(L1+2u*pad);            //full length of vecs in X1 including padding
    const int N2 = (int)(dil*(L2-1u)) + 1;      //full length of vec in X2 including dilation
    const int inc = (int)str - (int)(dil*L2);   //fixed increment for X1 below
    size_t w=0u, W;                             //current frame and total frames
    int ss, es;                                 //current start-samp, end-samp
    float sm;                                   //intermediate sum

    //Set W (L_out)
    W = 1u + (size_t)(N1-N2)/str;

    //Set ss, es (this currently matches xcorr and conv, i.e. Octave-like conventions)
    ss = -(int)pad; es = ss + N2 - 1;
    //fprintf(stderr,"L1=%lu, L2=%lu, N1=%d, N2=%d, W=%lu\n",L1,L2,N1,N2,W);

    //Don't flip X2
    X2 += L2 - 1u;

    if (W==0u) {}
    else if (L1==N)
    {
        //X2 before first samp of X1
        while (es<0 && w<W) { *Y++ = 0.0; es+=str; ++w; }
        ss = es - N2 + 1;
int cnt=0;
        //X2 overlaps first samp of X1
        while (ss<0 && es<(int)L1 && w<W)
        {
            sm = 0.0f; X1 += es; cnt += es;
            for (int n=es; n>=0; n-=dil, X1-=dil, cnt-=dil, --X2) { sm += *X1 * *X2; }
            *Y++ = sm;
            X1 += dil - (size_t)es%dil;
            cnt += dil - (size_t)es%dil;
            X2 += 1 + es/(int)dil;
            ss+=str; es+=str; ++w;
        }
        X2 -= L2 - 1u;
        fprintf(stderr,"ss=%d, cnt=%d\n",ss,cnt);

        if (N2>(int)L1) //X1 fully within X2
        {
            while (ss<=0 && w<W)
            {
                sm = 0.0f;
                for (int n=ss%(int)dil; n<es; n+=dil) { sm += X1[n] * X2[(n-ss)/(int)dil]; }
                *Y++ = sm;
                ss+=str; es+=str; ++w;
            }
            fprintf(stderr,"ss=%d, cnt=%d\n",ss,cnt);
            X1+=ss;
        }
        else            //X2 fully within X1
        {        
            while (es<(int)L1 && w<W)
            {
                sm = 0.0f;
                for (size_t l=0u; l<L2; ++l, X1+=dil, cnt+=dil, ++X2) { sm += *X1 * *X2; }
                *Y++ = sm;
                X1 += inc; cnt += inc; X2 -= L2;
                es+=str; ++w;
            }
            fprintf(stderr,"ss=%d, cnt=%d\n",ss,cnt);
            ss = es - N2 + 1;
        }

        //X2 overlaps end samp of X1
        while (ss<(int)L1 && w<W)
        {
            size_t l = 0; sm = 0.0f;
            for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=dil, ++X2) { sm += *X1 * *X2; }
            *Y++ = sm;
            X1 += (int)str - (int)(dil*l); X2 -= l;
            ss+=str; es+=str; ++w;
        }

        //X2 past end samp of X1
        while (w<W) { *Y++ = 0.0f; ++w; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L1, G = V/B;
        if (V) { return (int)G; }
    }

    //     for (size_t g=0u; g<G; ++g, X1+=B*(L1-1u), Y+=B*(W-1u))
    //     {
    //         for (size_t b=0u; b<B; ++b, ++X1, Y-=K*W-1u)
    //         {
    //             //X2 overlaps first samp of X1
    //             while (ss<0 && es<(int)L1 && w<W)
    //             {
    //                 sm = 0.0f; X1 += es*(int)K;
    //                 for (int n=0; n<es; ++n, X1-=K, --X2) { sm += *X1 * *X2; }
    //                 *Y = sm + *X1**X2; Y += K;
    //                 X2 += es;
    //                 ++ss; ++es; ++w;
    //             }
    //             X2 -= es;

    //             if (L2>L1)  //X1 fully within X2
    //             {
    //                 while (ss<=0 && w<W)
    //                 {
    //                     sm = 0.0f;
    //                     for (size_t l=0u; l<L1; ++l, X1+=K, ++X2) { sm += *X1 * *X2; }
    //                     *Y = sm; Y += K;
    //                     X1 -= K*L1; X2 -= L1+1u;
    //                     ++ss; ++w;
    //                 }
    //                 es = ss + (int)L2 - 1;
    //                 X1 += K; ++X2;
    //             }
    //             else        //X2 fully within X1
    //             {        
    //                 while (es<(int)L1 && w<W)
    //                 {
    //                     sm = 0.0f;
    //                     for (size_t l=0u; l<L2; ++l, X1+=K, ++X2) { sm += *X1 * *X2; }
    //                     *Y = sm; Y += K;
    //                     X1 += inc*(int)K; X2 -= L2;
    //                     ++es; ++w;
    //                 }
    //                 ss = es - (int)L2 + 1;
    //             }

    //             //Frames overlapping end samp
    //             while (ss<(int)L1 && w<W)
    //             {
    //                 sm = 0.0f;
    //                 for (int n=ss; n<(int)L1; ++n, X1+=K, ++X2) { sm += *X1 * *X2; }
    //                 *Y = sm; Y += K;
    //                 X1 += (int)K*(1-(int)L1+ss);
    //                 X2 -= (int)L1 - ss;
    //                 ++ss; ++w;
    //             }
    //             X1 -= (int)K*ss; X2 += L2 - 1u;
    //             ss -= (int)W; es = ss + (int)L2 - 1; w = 0;
    //         }
    //     }
    // }

    // if (W==0u) {}
    // else if (L1==N)
    // {
    //     const int Nb = (int)(dil*(L2-1u)) + 1;          //full length of X2 including dil
    //     const int inc = (int)str - (int)(dil*L2);       //fixed increment for X1 below
    //     size_t w = 0u;                                  //current frame
    //     int ss, es = c0 + (int)(dil*(L2/2u));           //current start-samp, end-samp
    //     float sm;                                       //intermediate sum

    //     //Don't flip X2
    //     X2 += L2 - 1u;

    //     //Frames before first samp
    //     while (es<0 && w<W) { *Y++ = 0.0f; ++w; es += str; }
    //     ss = es - Nb + 1;

    //     //Frames overlapping first samp
    //     while (ss<0 && es<(int)L1 && w<W)
    //     {
    //         sm = 0.0f; X1 += es;
    //         for (int n=es; n>=0; n-=(int)dil, X1-=dil, --X2) { sm = fmaf(*X2,*X1,sm); }
    //         *Y++ = sm;
    //         X2 += 1 + es/(int)dil; X1 += dil - (size_t)es%dil;
    //         ++w; ss += str; es += str;
    //     }
    //     X2 -= L2 - 1u;

    //     //In case L2>L1
    //     while (ss<0 && w<W)
    //     {
    //         int l = ss + (int)dil*((1-ss)/(int)dil);
    //         X2 += (l-ss)/(int)dil; X1 += l;
    //         sm = 0.0f;
    //         while (l<(int)L1) { sm = fmaf(*X2,*X1,sm); X1+=dil; l+=(int)dil; ++X2; }
    //         *Y++ = sm;
    //         X2 -= (l-ss)/(int)dil; X1 -= l;
    //         ++w; ss += str;
    //     }
    //     es = ss + Nb - 1;

    //     //Frames fully within sig
    //     while (es<(int)L1 && w<W)
    //     {
    //         sm = 0.0f;
    //         for (size_t l=0u; l<L2; ++l, X1+=dil, ++X2) { sm = fmaf(*X2,*X1,sm); }
    //         *Y++ = sm;
    //         X2 -= L2; X1 += inc;
    //         ++w; es += str;
    //     }
    //     ss = es - Nb + 1;

    //     //Frames overlapping end samp
    //     while (ss<(int)L1 && w<W)
    //     {
    //         size_t l = 0u; sm = 0.0f;
    //         for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=dil, ++X2) { sm = fmaf(*X2,*X1,sm); }
    //         *Y++ = sm;
    //         X2 -= l; X1 += (int)str - (int)(dil*l);
    //         ++w; ss += str;
    //     }

    //     //Frames after end samp
    //     while (w<W) { *Y++ = 0.0f; ++w; }
    // }
    // else
    // {
    //     const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
    //     const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
    //     const size_t V = N/L1, G = V/B;

    //     const int Nb = (int)(dil*(L2-1u)) + 1;              //full length of X2 including dil
    //     const int inc = (int)K*((int)str-(int)(dil*L2));    //fixed increment for X1 below
    //     size_t w;                                           //current frame
    //     int ss, es;                                         //current start-samp, end-samp
    //     float sm;                                           //intermediate sum

    //     //Don't flip X2
    //     X2 += L2 - 1u;

    //     for (size_t g=0u; g<G; ++g, X1+=B*(L1-1u), Y+=B*(W-1u))
    //     {
    //         for (size_t b=0u; b<B; ++b, X1-=K*L1-1u, Y-=K*W-1u)
    //         {
    //             //Start frames
    //             w = 0u; es = c0 + (int)(dil*(L2/2u));

    //             //Frames before first samp
    //             while (es<0 && w<W) { *Y = 0.0f; Y+=K; ++w; es += str; }
    //             ss = es - Nb + 1;

    //             //Frames overlapping first samp
    //             while (ss<0 && es<(int)L1 && w<W)
    //             {
    //                 sm = 0.0f; X1 += (int)K*es;
    //                 for (int n=es; n>=0; n-=(int)dil, X1-=K*dil, --X2) { sm = fmaf(*X2,*X1,sm); }
    //                 *Y = sm; Y += K;
    //                 X2 += 1 + es/(int)dil; X1 += K*(dil-(size_t)es%dil);
    //                 ++w; ss += str; es += str;
    //             }
    //             X2 -= L2 - 1u;

    //             //In case L2>L1
    //             while (ss<0 && w<W)
    //             {
    //                 int l = ss + (int)dil*((1-ss)/(int)dil);
    //                 X2 += (l-ss)/(int)dil; X1 += l*(int)K;
    //                 sm = 0.0f;
    //                 while (l<(int)L1) { sm = fmaf(*X2,*X1,sm); X1+=dil*K; l+=(int)dil; ++X2; }
    //                 *Y = sm; Y += K;
    //                 X2 -= (l-ss)/(int)dil; X1 -= l*(int)K;
    //                 ++w; ss += str;
    //             }
    //             es = ss + Nb - 1;

    //             //Frames fully within sig
    //             while (es<(int)L1 && w<W)
    //             {
    //                 sm = 0.0f;
    //                 for (size_t l=0u; l<L2; ++l, X1+=K*dil, ++X2) { sm = fmaf(*X2,*X1,sm); }
    //                 *Y = sm; Y += K;
    //                 X2 -= L2; X1 += inc;
    //                 ++w; es += str;
    //             }
    //             ss = es - Nb + 1;

    //             //Frames overlapping end samp
    //             while (ss<(int)L1 && w<W)
    //             {
    //                 size_t l = 0u; sm = 0.0f;
    //                 for (int n=ss; n<(int)L1; n+=dil, X1+=K*dil, ++X2, ++l) { sm = fmaf(*X2,*X1,sm); }
    //                 *Y = sm; Y += K;
    //                 X2 -= l; X1 += (int)K*((int)str-(int)(dil*l));
    //                 ++w; ss += str;
    //             }
    //             X2 += L2 - 1u; X1 -= K*((size_t)ss-L1);

    //             //Frames after end samp
    //             while (w<W) { *Y = 0.0f; Y+=K; ++w; }
    //         }
    //     }
    // }

    return 0;
}


// int xcorr1d_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t pad, const size_t str, const size_t dil, const size_t dim)
// {
//     if (dim>3u) { fprintf(stderr,"error in xcorr1d_d: dim must be in [0 3]\n"); return 1; }
//     if (str<1u) { fprintf(stderr,"error in xcorr1d_d: str (stride) must be positive\n"); return 1; }
//     if (dil<1u) { fprintf(stderr,"error in xcorr1d_d: dil (dilation) must be positive\n"); return 1; }

//     const size_t N = R*C*S*H;
//     const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
//     if (N<1u) { fprintf(stderr,"error in xcorr1d_d: N (total length of X1) must be positive\n"); return 1; }
//     if (L1<1u) { fprintf(stderr,"error in xcorr1d_d: L1 (length of vecs in X1) must be positive\n"); return 1; }
//     if (c0>=(int)L1) { fprintf(stderr,"error in xcorr1d_d: c0 (start samp of 1st frame) must be < L1 (length of vecs in X1)\n"); return 1; }

//     if (W==0u) {}
//     else if (L1==N)
//     {
//         const int Nb = (int)(dil*(L2-1u)) + 1;      //full length of X2 including dil
//         const int inc = (int)str - (int)(dil*L2);   //fixed increment for X1 below
//         size_t w = 0u;                              //current frame
//         int ss, es = c0 + Nb/2;                     //current start-samp, end-samp
//         double sm;                                  //intermediate sum

//         //Don't flip X2
//         X2 += L2 - 1u;

//         //Frames before first samp
//         while (es<0 && w<W) { *Y++ = 0.0; es += str; ++w; }
//         ss = es - Nb + 1;

//         //Frames overlapping first samp
//         while (ss<0 && es<(int)L1 && w<W)
//         {
//             sm = 0.0; X1 += es;
//             for (int n=es; n>=0; n-=(int)dil, X1-=dil, --X2) { sm = fma(*X2,*X1,sm); }
//             *Y++ = sm;
//             X2 += 1 + es/(int)dil; X1 += dil - (size_t)es%dil;
//             ss += str; es += str; ++w;
//         }
//         X2 -= L2 - 1u;

//         //In case L2>L1
//         while (ss<0 && w<W)
//         {
//             int l = ss + (int)dil*((1-ss)/(int)dil);
//             X2 += (l-ss)/(int)dil; X1 += l;
//             sm = 0.0;
//             while (l<(int)L1) { sm = fma(*X2,*X1,sm); X1+=dil; l+=(int)dil; ++X2; }
//             *Y++ = sm;
//             X2 -= (l-ss)/(int)dil; X1 -= l;
//             ss += str; ++w;
//         }
//         es = ss + Nb - 1;

//         //Frames fully within sig
//         while (es<(int)L1 && w<W)
//         {
//             sm = 0.0;
//             for (size_t l=0u; l<L2; ++l, X1+=dil, ++X2) { sm = fma(*X2,*X1,sm); }
//             *Y++ = sm;
//             X2 -= L2; X1 += inc;
//             es += str; ++w;
//         }
//         ss = es - Nb + 1;

//         //Frames overlapping end samp
//         while (ss<(int)L1 && w<W)
//         {
//             size_t l = 0u; sm = 0.0;
//             for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=dil, ++X2) { sm = fma(*X2,*X1,sm); }
//             *Y++ = sm;
//             X2 -= l; X1 += (int)str - (int)(dil*l);
//             ss += str; ++w;
//         }

//         //Frames after end samp
//         while (w<W) { *Y++ = 0.0; ++w; }
//     }
//     else
//     {
//         const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
//         const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
//         const size_t V = N/L1, G = V/B;
//         fprintf(stderr,"K=%lu, B=%lu, V=%lu, G=%lu \n",K,B,V,G);

//         const int Nb = (int)(dil*(L2-1u)) + 1;              //full length of X2 including dil
//         const int inc = (int)K*((int)str-(int)(dil*L2));    //fixed increment for X1 below
//         size_t w;                                           //current frame
//         int ss, es = c0 + Nb/2;                             //current start-samp, end-samp
//         double sm;                                          //intermediate sum

//         //Don't flip X2
//         X2 += L2 - 1u; int cnt = 0;

//         for (size_t g=0u; g<G; ++g, X1+=B*(L1-1u), Y+=B*(W-1u))
//         {
//             for (size_t bs=0u; bs<B; ++bs, X1-=K*L1-1u, Y-=K*W-1u)
//             {
//                 //Start frames
//                 ss = c0; es = ss + Nb - 1; w = 0u;

//                 //Frames before first samp
//                 while (es<0 && w<W) { *Y = 0.0; Y+=K; es+=str; ++w; }
//                 ss = es - Nb + 1;

//                 //Frames overlapping first samp
//                 while (ss<0 && es<(int)L1 && w<W)
//                 {
//                     sm = 0.0; X1 += (int)K*es;
//                     for (int n=es; n>=0; n-=(int)dil, X1-=K*dil, --X2) { sm = fma(*X2,*X1,sm); }
//                     *Y = sm; Y += K;
//                     X2 += 1 + es/(int)dil; X1 += K*(dil-(size_t)es%dil);
//                     ss += str; es += str; ++w;
//                 }
//                 X2 -= L2 - 1u; fprintf(stderr,"Y=%g\n",(double)(*(Y-K)));

//                 //In case L2>L1
//                 while (ss<0 && w<W)
//                 {
//                     // int lb = 0, nb = 0, lx; sm = 0.0;
//                     // while (lb<L2)
//                     // {
//                     //     if (nb>-ss && nb<(int)L1)
//                     //     {
//                     //         sm = fma(*X2,*X1,sm);
//                     //     }
//                     //     ++lb; nb+=(int)dil;
//                     // }
//                     int b = (1-ss)/(int)dil, l = ss + (int)dil*b;
//                     //int l = ss + (int)dil*((1-ss)/(int)dil);
//                     X2 += b; X1 += l*(int)K;
//                     sm = 0.0;
//                     while (l<(int)L1) { sm = fma(*X2,*X1,sm); X1+=dil*K; l+=(int)dil; ++X2; ++b; }
//                     *Y = sm; Y += K;
//                     X2 -= b; X1 -= l*(int)K;
//                     // X2 -= (l-ss)/(int)dil; X1 -= l*(int)K;
//                     ss += str; ++w;
//                 }
//                 es = ss + Nb - 1;

//                 //Frames fully within sig
//                 while (es<(int)L1 && w<W)
//                 {
//                     sm = 0.0;
//                     for (size_t l=0u; l<L2; ++l, X1+=K*dil, ++X2) { sm = fma(*X2,*X1,sm); }
//                     *Y = sm; Y += K;
//                     X2 -= L2; X1 += inc;
//                     es += str; ++w;
//                 }
//                 ss = es - Nb + 1;

//                 //Frames overlapping end samp
//                 while (ss<(int)L1 && w<W)
//                 {
//                     size_t l = 0u; sm = 0.0;
//                     for (int n=ss; n<(int)L1; n+=dil, X1+=K*dil, ++X2, ++l) { sm = fma(*X2,*X1,sm); }
//                     *Y = sm; Y += K;
//                     X2 -= l; X1 += (int)K*((int)str-(int)(dil*l));
//                     ss += str; ++w;
//                 }
//                 X2 += L2 - 1u; X1 -= K*dil*((size_t)ss-L1);

//                 //Frames after end samp
//                 while (w<W) { *Y = 0.0; Y+=K; ++w; }
//             }
//         }
//     }

//     return 0;
// }


// int xcorr1d_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t pad, const size_t str, const size_t dil, const size_t dim)
// {
//     if (dim>3u) { fprintf(stderr,"error in xcorr1d_c: dim must be in [0 3]\n"); return 1; }
//     if (str<1u) { fprintf(stderr,"error in xcorr1d_c: str (stride) must be positive\n"); return 1; }
//     if (dil<1u) { fprintf(stderr,"error in xcorr1d_c: dil (dilation) must be positive\n"); return 1; }

//     const size_t N = R*C*S*H;
//     const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
//     if (N<1u) { fprintf(stderr,"error in xcorr1d_c: N (total length of X1) must be positive\n"); return 1; }
//     if (L1<1u) { fprintf(stderr,"error in xcorr1d_c: L1 (length of vecs in X1) must be positive\n"); return 1; }
//     if (c0>=(int)L1) { fprintf(stderr,"error in xcorr1d_c: c0 (center samp of 1st frame) must be < L1 (length of vecs in X1)\n"); return 1; }

//     if (W==0u) {}
//     else if (L1==N)
//     {
//         const int Nb = (int)(dil*(L2-1u)) + 1;          //full length of X2 including dil
//         const int inc = 2*((int)str-(int)(dil*L2));     //fixed increment for X1 below
//         size_t w = 0u;                                  //current frame
//         int ss, es = c0 + (int)(dil*(L2/2u));           //current start-samp, end-samp
//         float smr, smi;                                 //intermediate sums

//         //Don't flip X2
//         X2 += 2u*L2 - 2u;

//         //Frames before first samp
//         while (es<0 && w<W) { *Y++ = 0.0f; *Y++ = 0.0f; ++w; es += str; }
//         ss = es - Nb + 1;

//         //Frames overlapping first samp
//         while (ss<0 && w<W)
//         {
//             smr = smi = 0.0f; X1 += 2*es;
//             for (int n=es; n>=0; n-=(int)dil, X1-=2u*dil, X2-=2)
//             {
//                 smr += *X2**X1 - *(X2+1)**(X1+1);
//                 smi += *X2**(X1+1) + *(X2+1)**X1;
//             }
//             *Y++ = smr; *Y++ = smi;
//             X2 += 2*(1+es/(int)dil);
//             X1 += 2u*(dil-(size_t)es%dil);
//             ++w; ss += str; es += str;
//         }
//         X2 -= 2u*L2 - 2u; X1 += 2*ss;

//         //Frames fully within sig
//         while (es<(int)L1 && w<W)
//         {
//             smr = smi = 0.0f;
//             for (size_t l=0u; l<L2; ++l, X1+=2u*dil, X2+=2)
//             {
//                 smr += *X2**X1 - *(X2+1)**(X1+1);
//                 smi += *X2**(X1+1) + *(X2+1)**X1;
//             }
//             *Y++ = smr; *Y++ = smi;
//             X2 -= 2u*L2; X1 += inc;
//             ++w; es += str;
//         }
//         ss = es - Nb + 1;

//         //Frames overlapping end samp
//         while (ss<(int)L1 && w<W)
//         {
//             size_t l = 0u; smr = smi = 0.0f;
//             for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=2u*dil, X2+=2)
//             {
//                 smr += *X2**X1 - *(X2+1)**(X1+1);
//                 smi += *X2**(X1+1) + *(X2+1)**X1;
//             }
//             *Y++ = smr; *Y++ = smi;
//             X2 -= 2u*l; X1 += 2*((int)str-(int)(dil*l));
//             ++w; ss += str;
//         }

//         //Frames after end samp
//         while (w<W) { *Y++ = 0.0f; *Y++ = 0.0f; ++w; }
//     }
//     else
//     {
//         const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
//         const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
//         const size_t V = N/L1, G = V/B;

//         const int Nb = (int)(dil*(L2-1u)) + 1;              //full length of X2 including dil
//         const int inc = 2*(int)K*((int)str-(int)(dil*L2));  //fixed increment for X1 below
//         size_t w;                                           //current frame
//         int ss, es;                                         //current start-samp, end-samp
//         float smr, smi;                                     //intermediate sums

//         //Don't flip X2
//         X2 += 2u*L2 - 2u;

//         for (size_t g=0u; g<G; ++g, X1+=2u*B*(L1-1u), Y+=2u*B*(W-1u))
//         {
//             for (size_t b=0u; b<B; ++b, X1-=2u*K*L1-2u, Y-=2u*K*W-2u)
//             {
//                 //Start frames
//                 w = 0u; es = c0 + (int)(dil*(L2/2u));

//                 //Frames before first samp
//                 while (es<0 && w<W) { *Y = *(Y+1) = 0.0f; Y+=2u*K; ++w; es += str; }
//                 ss = es - Nb + 1;

//                 //Frames overlapping first samp
//                 while (ss<0 && w<W)
//                 {
//                     smr = smi = 0.0f; X1 += 2*(int)K*es;
//                     for (int n=es; n>=0; n-=(int)dil, X1-=2u*K*dil, X2-=2)
//                     {
//                         smr += *X2**X1 - *(X2+1)**(X1+1);
//                         smi += *X2**(X1+1) + *(X2+1)**X1;
//                     }
//                     *Y = smr; *(Y+1) = smi; Y += 2u*K;
//                     X2 += 2*(1+es/(int)dil);
//                     X1 += 2u*K*(dil-(size_t)es%dil);
//                     ++w; ss += str; es += str;
//                 }
//                 X2 -= 2u*L2 - 2u; X1 += 2*ss*(int)K;

//                 //Frames fully within sig
//                 while (es<(int)L1 && w<W)
//                 {
//                     smr = smi = 0.0f;
//                     for (size_t l=0u; l<L2; ++l, X1+=2u*K*dil, X2+=2)
//                     {
//                         smr += *X2**X1 - *(X2+1)**(X1+1);
//                         smi += *X2**(X1+1) + *(X2+1)**X1;
//                     }
//                     *Y = smr; *(Y+1) = smi; Y += 2u*K;
//                     X2 -= 2u*L2; X1 += inc;
//                     ++w; es += str;
//                 }
//                 ss = es - Nb + 1;

//                 //Frames overlapping end samp
//                 while (ss<(int)L1 && w<W)
//                 {
//                     size_t l = 0u; smr = smi = 0.0f;
//                     for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=2u*K*dil, X2+=2)
//                     {
//                         smr += *X2**X1 - *(X2+1)**(X1+1);
//                         smi += *X2**(X1+1) + *(X2+1)**X1;
//                     }
//                     *Y = smr; *(Y+1) = smi; Y += 2u*K;
//                     X2 -= 2u*l; X1 += 2*(int)K*((int)str-(int)(dil*l));
//                     ++w; ss += str;
//                 }
//                 X2 += 2u*L2 - 2u; X1 -= 2u*K*((size_t)ss-L1);

//                 //Frames after end samp
//                 while (w<W) { *Y = *(Y+1) = 0.0f; Y+=2u*K; ++w; }
//             }
//         }
//     }

//     return 0;
// }


// int xcorr1d_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t pad, const size_t str, const size_t dil, const size_t dim)
// {
//     if (dim>3u) { fprintf(stderr,"error in xcorr1d_z: dim must be in [0 3]\n"); return 1; }
//     if (str<1u) { fprintf(stderr,"error in xcorr1d_z: str (stride) must be positive\n"); return 1; }
//     if (dil<1u) { fprintf(stderr,"error in xcorr1d_z: dil (dilation) must be positive\n"); return 1; }

//     const size_t N = R*C*S*H;
//     const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
//     if (N<1u) { fprintf(stderr,"error in xcorr1d_z: N (total length of X1) must be positive\n"); return 1; }
//     if (L1<1u) { fprintf(stderr,"error in xcorr1d_z: L1 (length of vecs in X1) must be positive\n"); return 1; }
//     if (c0>=(int)L1) { fprintf(stderr,"error in xcorr1d_z: c0 (center samp of 1st frame) must be < L1 (length of vecs in X1)\n"); return 1; }

//     if (W==0u) {}
//     else if (L1==N)
//     {
//         const int Nb = (int)(dil*(L2-1u)) + 1;          //full length of X2 including dil
//         const int inc = 2*((int)str-(int)(dil*L2));     //fixed increment for X1 below
//         size_t w = 0u;                                  //current frame
//         int ss, es = c0 + (int)(dil*(L2/2u));           //current start-samp, end-samp
//         double smr, smi;                                //intermediate sums

//         //Don't flip X2
//         X2 += 2u*L2 - 2u;

//         //Frames before first samp
//         while (es<0 && w<W) { *Y++ = 0.0; *Y++ = 0.0; ++w; es += str; }
//         ss = es - Nb + 1;

//         //Frames overlapping first samp
//         while (ss<0 && w<W)
//         {
//             smr = smi = 0.0; X1 += 2*es;
//             for (int n=es; n>=0; n-=(int)dil, X1-=2u*dil, X2-=2)
//             {
//                 smr += *X2**X1 - *(X2+1)**(X1+1);
//                 smi += *X2**(X1+1) + *(X2+1)**X1;
//             }
//             *Y++ = smr; *Y++ = smi;
//             X2 += 2*(1+es/(int)dil);
//             X1 += 2u*(dil-(size_t)es%dil);
//             ++w; ss += str; es += str;
//         }
//         X2 -= 2u*L2 - 2u; X1 += 2*ss;

//         //Frames fully within sig
//         while (es<(int)L1 && w<W)
//         {
//             smr = smi = 0.0;
//             for (size_t l=0u; l<L2; ++l, X1+=2u*dil, X2+=2)
//             {
//                 smr += *X2**X1 - *(X2+1)**(X1+1);
//                 smi += *X2**(X1+1) + *(X2+1)**X1;
//             }
//             *Y++ = smr; *Y++ = smi;
//             X2 -= 2u*L2; X1 += inc;
//             ++w; es += str;
//         }
//         ss = es - Nb + 1;

//         //Frames overlapping end samp
//         while (ss<(int)L1 && w<W)
//         {
//             size_t l = 0u; smr = smi = 0.0;
//             for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=2u*dil, X2+=2)
//             {
//                 smr += *X2**X1 - *(X2+1)**(X1+1);
//                 smi += *X2**(X1+1) + *(X2+1)**X1;
//             }
//             *Y++ = smr; *Y++ = smi;
//             X2 -= 2u*l; X1 += 2*((int)str-(int)(dil*l));
//             ++w; ss += str;
//         }

//         //Frames after end samp
//         while (w<W) { *Y++ = 0.0; *Y++ = 0.0; ++w; }
//     }
//     else
//     {
//         const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
//         const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
//         const size_t V = N/L1, G = V/B;

//         const int Nb = (int)(dil*(L2-1u)) + 1;              //full length of X2 including dil
//         const int inc = 2*(int)K*((int)str-(int)(dil*L2));  //fixed increment for X1 below
//         size_t w;                                           //current frame
//         int ss, es;                                         //current start-samp, end-samp
//         double smr, smi;                                    //intermediate sums

//         //Don't flip X2
//         X2 += 2u*L2 - 2u;

//         for (size_t g=0u; g<G; ++g, X1+=2u*B*(L1-1u), Y+=2u*B*(W-1u))
//         {
//             for (size_t b=0u; b<B; ++b, X1-=2u*K*L1-2u, Y-=2u*K*W-2u)
//             {
//                 //Start frames
//                 w = 0u; es = c0 + (int)(dil*(L2/2u));

//                 //Frames before first samp
//                 while (es<0 && w<W) { *Y = *(Y+1) = 0.0; Y+=2u*K; ++w; es += str; }
//                 ss = es - Nb + 1;

//                 //Frames overlapping first samp
//                 while (ss<0 && w<W)
//                 {
//                     smr = smi = 0.0; X1 += 2*(int)K*es;
//                     for (int n=es; n>=0; n-=(int)dil, X1-=2u*K*dil, X2-=2)
//                     {
//                         smr += *X2**X1 - *(X2+1)**(X1+1);
//                         smi += *X2**(X1+1) + *(X2+1)**X1;
//                     }
//                     *Y = smr; *(Y+1) = smi; Y += 2u*K;
//                     X2 += 2*(1+es/(int)dil);
//                     X1 += 2u*K*(dil-(size_t)es%dil);
//                     ++w; ss += str; es += str;
//                 }
//                 X2 -= 2u*L2 - 2u; X1 += 2*ss*(int)K;

//                 //Frames fully within sig
//                 while (es<(int)L1 && w<W)
//                 {
//                     smr = smi = 0.0;
//                     for (size_t l=0u; l<L2; ++l, X1+=2u*K*dil, X2+=2)
//                     {
//                         smr += *X2**X1 - *(X2+1)**(X1+1);
//                         smi += *X2**(X1+1) + *(X2+1)**X1;
//                     }
//                     *Y = smr; *(Y+1) = smi; Y += 2u*K;
//                     X2 -= 2u*L2; X1 += inc;
//                     ++w; es += str;
//                 }
//                 ss = es - Nb + 1;

//                 //Frames overlapping end samp
//                 while (ss<(int)L1 && w<W)
//                 {
//                     size_t l = 0u; smr = smi = 0.0;
//                     for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=2u*K*dil, X2+=2)
//                     {
//                         smr += *X2**X1 - *(X2+1)**(X1+1);
//                         smi += *X2**(X1+1) + *(X2+1)**X1;
//                     }
//                     *Y = smr; *(Y+1) = smi; Y += 2u*K;
//                     X2 -= 2u*l; X1 += 2*(int)K*((int)str-(int)(dil*l));
//                     ++w; ss += str;
//                 }
//                 X2 += 2u*L2 - 2u; X1 -= 2u*K*((size_t)ss-L1);

//                 //Frames after end samp
//                 while (w<W) { *Y = *(Y+1) = 0.0; Y+=2u*K; ++w; }
//             }
//         }
//     }

//     return 0;
// }


#ifdef __cplusplus
}
}
#endif
