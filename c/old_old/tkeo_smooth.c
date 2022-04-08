//Teager-Kaiser Energy Operator (TKEO) according to de Matos [2018],
//which is computed using outputs of smooth_diff and smooth_diffdiff:
//y[n] = x'[n]*x'[n] - x[n]*x''[n]

//For complex input, the output is real, and is equal to the sum
//of the TKEOs of the real and imaginary parts, as proven
//by Hamila et al. [1999], as cited in de Matos [2018].

//See de Matos MC. 2018. Seismic attributes from the complex Teager-Kaiser energy.
//and its excellent reference:
//Holoborodko P. 2008. Smooth noise robust differentiators. www.holoborodko.com.

//For real X, the exact method of de Matos [2018] would be:
//analytic_sig X | tkeo_smooth > Y\n";
//where Y is the "complex TK energy" (but is real-valued).
//(The input to tkeo_smooth will be complex, so sum rule above applies.)

//This paper also gives the complex TKV (variational Teager-Kaiser) method,
//which is a simple extension of this, but would require looking at a few of his references.
//Basically, a further Hilbert transform of the (real-valued) "complex TK energy" gives the TKV.

//I decided to use signal extrapolation by signal reversal for the edge samps.
//This is exactly like window_univar if stp=1;
//but for stp>1, it keeps the first center samp always at 0.
//Thus, the output vecs have length Ly = 1 + (Lx-1)/stp;
//This is the same as xcorr1d (which emulates PyTorch Conv1d) with pad = dil = 0:
//L_out = floor[1 + (L_in + 2*pad - dil*(L2-1) - 1)/stride].

#include <stdio.h>
#include <stdlib.h>
#include "xcorr.c"
#include "../../math/c/square.c"
#include "../../math/c/times.c"
#include "../../math/c/minus.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int tkeo_smooth_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t stp);
int tkeo_smooth_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t stp);
int tkeo_smooth_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t stp);
int tkeo_smooth_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t stp);


int tkeo_smooth_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t stp)
{
    if (dim>3u) { fprintf(stderr,"error in tkeo_smooth_s: dim must be in [0 3]\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in tkeo_smooth_s: stp must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = 1u + (Lx-1u)/stp;
    if (Lx<11u) { fprintf(stderr,"error in tkeo_smooth_s: Lx (length of vecs in X) must be >=11 (length of smooth_diff)\n"); return 1; }

    //Init 1st-order differentiator (n=2, N=11 of smooth_diff)
    const float D[11] = {-1.0f/512.0f,-8.0f/512.0f,-27.0f/512.0f,-48.0f/512.0f,-42.0f/512.0,0.0f,42.0f/512.0f,48.0f/512.0f,27.0/512.0f,8.0f/512.0f,1.0f/512.0f};

    //Init 2nd-order differentiator (n=3, N=5 of smooth_diffdiff)
    const float DD[9] = {-7.0f/192.0f,12.0f/192.0f,52.0f/192.0f,-12.0f/192.0f,-90.0f/192.0f,-12.0f/192.0f,52.0f/192.0f,12.0f/192.0f,-7.0f/192.0f};

    //Init window_univar
    int ss = -11/2, n, prev_n = 0;          //start-samp, current/prev samps

    //Init Xd for intermediate output of diff and diffdiff
    float *Xd  = (float *)calloc(Ly,sizeof(float));
    if (!Xd) { fprintf(stderr,"error in tkeo_smooth_s: problem with calloc. "); perror("calloc"); return 1; }

    if (N==0u) {}
    else if (Lx==N)
    {
        //Smooth diff of X
        //xcorr_s (Xd,X,D,R,C,S,H,iscolmajor,11u,"same",dim); //if stp==1
        //xcorr1d_s (Xd,X,D,R,C,S,H,iscolmajor,11u,0u,stp,0u,dim);
        xcorr1_s (Xd,X,D,R,C,S,H,iscolmajor,11u,5u,stp,Ly,dim);

        // int cnt = 0;
        // for (size_t w=Ly; w>0u; --w, ss+=(int)stp, ++Xd)
        // {
        //     if (ss<0 || ss>(int)Lx-11)
        //     {
        //         for (int s=ss; s<ss+11; ++s)
        //         {
        //             n = s; //This ensures extrapolation by signal reversal to any length
        //             while (n<0 || n>=(int)Lx) { n = (n<0) ? -n-1 : (n<(int)Lx) ? n : 2*(int)Lx-1-n; }
        //             X += n - prev_n; cnt += n - prev_n;
        //             *Xd += *X * D[s-ss];
        //             prev_n = n;
        //         }
        //     }
        //     else
        //     {
        //         X += ss - prev_n; cnt += ss - prev_n;
        //         for (size_t l=0u; l<11u; ++l, ++X, ++cnt) { *Xd += *X * D[l]; }
        //         X -= 11u - stp; prev_n = ss + (int)stp; cnt -= 11u - stp;
        //     }
        // }
        // fprintf(stderr,"cnt=%d\n",cnt);
        // X -= Ly - 5u;

        //Get Xd^2 (1st term of smooth TKEO)
        //square_s (Y,Xd,N);
        Y += Ly; Xd += Ly;
        for (size_t w=Ly; w>0u; --w) { --Xd; --Y; *Y = *Xd * *Xd; }

        //Smooth diffdiff of X
        //xcorr_s (Xd,X,DD,R,C,S,H,iscolmajor,9u,"same",dim); //if stp==1
        //xcorr1d_s (Xd,X,DD,R,C,S,H,iscolmajor,9u,0u,stp,0u,dim);
        xcorr1_s (Xd,X,DD,R,C,S,H,iscolmajor,9u,4u,stp,Ly,dim);
        // ss = -9/2; prev_n = 0;
        // for (size_t w=Ly; w>0u; --w, ss+=(int)stp, ++Xd)
        // {
        //     *Xd = 0.0f;
        //     if (ss<0 || ss>(int)Lx-9)
        //     {
        //         for (int s=ss; s<ss+9; ++s)
        //         {
        //             n = s; //This ensures extrapolation by signal reversal to any length
        //             while (n<0 || n>=(int)Lx) { n = (n<0) ? -n-1 : (n<(int)Lx) ? n : 2*(int)Lx-1-n; }
        //             X += n - prev_n;
        //             *Xd += *X * DD[s-ss];
        //             prev_n = n;
        //         }
        //     }
        //     else
        //     {
        //         X += ss - prev_n;
        //         for (size_t l=0u; l<9u; ++l, ++X) { *Xd += *X * DD[l]; }
        //         X -= 9u - stp; prev_n = ss + (int)stp;
        //     }
        // }
        // X += 4u;

        //Get X*Xdd (2nd term of smooth TKEO)
        //times_inplace_s(Xd,X,L,1,1,1,L,1,1,1,0);

        //Get smooth TKEO
        //minus_inplace_s(Y,Xd,L,1,1,1,L,1,1,1,0);
        Y += Ly;
        for (size_t w=Ly; w>0u; --w) { *--Y -= *--X * *--Xd; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                //Smooth diff of X
                ss = -11/2; prev_n = 0;
                for (size_t l=Lx; l>0u; --l, ++ss, ++Xd)
                {
                    *Xd = 0.0f;
                    if (ss<0 || ss>(int)Lx-11)
                    {
                        for (int s=ss; s<ss+11; ++s)
                        {
                            n = s; //This ensures extrapolation by signal reversal to any length
                            while (n<0 || n>=(int)Lx) { n = (n<0) ? -n-1 : (n<(int)Lx) ? n : 2*(int)Lx-1-n; }
                            X += n - prev_n;
                            *Xd += *X * D[s-ss];
                            prev_n = n;
                        }
                    }
                    else
                    {
                        X += ss - prev_n;
                        for (size_t l=0u; l<11u; ++l, ++X) { *Xd += *X * D[l]; }
                        X -= 11u - stp; prev_n = ss + 1;
                    }
                }
                X -= Ly - 5u;

                //Get Xd^2 (1st term of smooth TKEO)
                Y += Ly;
                for (size_t w=Ly; w>0u; --w) { --Xd; --Y; *Y = *Xd * *Xd; }

                //Smooth diffdiff of X
                ss = -9/2; prev_n = 0;
                for (size_t w=Ly; w>0u; --w, ++ss, ++Xd)
                {
                    *Xd = 0.0f;
                    if (ss<0 || ss>(int)Lx-9)
                    {
                        for (int s=ss; s<ss+9; ++s)
                        {
                            n = s; //This ensures extrapolation by signal reversal to any length
                            while (n<0 || n>=(int)Lx) { n = (n<0) ? -n-1 : (n<(int)Lx) ? n : 2*(int)Lx-1-n; }
                            X += n - prev_n;
                            *Xd += *X * DD[s-ss];
                            prev_n = n;
                        }
                    }
                    else
                    {
                        X += ss - prev_n;
                        for (size_t l=0u; l<9u; ++l, ++X) { *Xd += *X * DD[l]; }
                        X -= 9u - stp; prev_n = ss + 1;
                    }
                }
                X -= Ly - 4u; Xd -= Ly;

                //Get X*Xdd (2nd term of smooth TKEO) and subtrac to get smooth TKEO
                for (size_t w=Ly; w>0u; --w, ++X, ++Xd, ++Y) { *Y -= *X * *Xd; }
                Xd -= Ly;
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(Ly-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*Lx-1u, Y-=K*Ly-1u)
                {
                    //Smooth diff of X
                    ss = -11/2; prev_n = 0;
                    for (size_t w=Ly; w>0u; --w, ++ss, ++Xd)
                    {
                        *Xd = 0.0f;
                        if (ss<0 || ss>(int)Lx-11)
                        {
                            for (int s=ss; s<ss+11; ++s)
                            {
                                n = s; //This ensures extrapolation by signal reversal to any length
                                while (n<0 || n>=(int)Lx) { n = (n<0) ? -n-1 : (n<(int)Lx) ? n : 2*(int)Lx-1-n; }
                                X += (int)K * (n-prev_n);
                                *Xd += *X * D[s-ss];
                                prev_n = n;
                            }
                        }
                        else
                        {
                            X += (int)K * (ss-prev_n);
                            for (size_t l=0u; l<11u; ++l, X+=K) { *Xd += *X * D[l]; }
                            X -= K*(11u-stp); prev_n = ss + 1;
                        }
                    }
                    X -= K * (Ly-5u);

                    //Get Xd^2 (1st term of smooth TKEO)
                    Y += K*Ly;
                    for (size_t w=Ly; w>0u; --w) { --Xd; Y-=K; *Y = *Xd * *Xd; }

                    //Smooth diffdiff of X
                    ss = -9/2; prev_n = 0;
                    for (size_t w=Ly; w>0u; --w, ++ss, ++Xd)
                    {
                        *Xd = 0.0f;
                        if (ss<0 || ss>(int)Lx-9)
                        {
                            for (int s=ss; s<ss+9; ++s)
                            {
                                n = s; //This ensures extrapolation by signal reversal to any length
                                while (n<0 || n>=(int)Lx) { n = (n<0) ? -n-1 : (n<(int)Lx) ? n : 2*(int)Lx-1-n; }
                                X += (int)K * (n-prev_n);
                                *Xd += *X * DD[s-ss];
                                prev_n = n;
                            }
                        }
                        else
                        {
                            X += (int)K * (ss-prev_n);
                            for (size_t l=0u; l<9u; ++l, X+=K) { *Xd += *X * DD[l]; }
                            X -= K*(9u-stp); prev_n = ss + 1;
                        }
                    }
                    X -= K*(Ly-4u); Xd -= Ly;

                    //Get X*Xdd (2nd term of smooth TKEO) and subtrac to get smooth TKEO
                    for (size_t w=Ly; w>0u; --w, X+=K, ++Xd, Y+=K) { *Y -= *X * *Xd; }
                    Xd -= Ly;
                }
            }
        }
    }
    free(Xd);

    return 0;
}


int tkeo_smooth_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t stp)
{
    if (dim>3u) { fprintf(stderr,"error in tkeo_smooth_d: dim must be in [0 3]\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in tkeo_smooth_d: stp must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = 1u + (Lx-1u)/stp;
    if (Lx<11u) { fprintf(stderr,"error in tkeo_smooth_d: L (length of vecs in X) must be >=11 (length of smooth_diff)\n"); return 1; }

    //Init 1st-order differentiator (n=2, N=11 of smooth_diff)
    const double D[11] = {-1.0/512.0,-8.0/512.0,-27.0/512.0,-48.0/512.0,-42.0/512.0,0.0,42.0/512.0,48.0/512.0,27.0/512.0,8.0/512.0,1.0/512.0};

    //Init 2nd-order differentiator (n=3, N=5 of smooth_diffdiff)
    const double DD[9] = {-7.0/192.0,12.0/192.0,52.0/192.0,-12.0/192.0,-90.0/192.0,-12.0/192.0,52.0/192.0,12.0/192.0,-7.0/192.0};

    //Init window_univar
    int ss = -11/2, n, prev_n = 0;          //start-samp, current/prev samps

    //Init Xd for intermediate output of diff and diffdiff
    double *Xd  = (double *)calloc(Ly,sizeof(double));
    if (!Xd) { fprintf(stderr,"error in tkeo_smooth_d: problem with calloc. "); perror("calloc"); return 1; }

    if (N==0u) {}
    // else if (Lx==N)
    // {
    //     //Smooth diff of X
    //     for (size_t l=Lx; l>0u; --l, ++ss, ++Xd)
    //     {
    //         if (ss<0 || ss>(int)Lx-11)
    //         {
    //             for (int s=ss; s<ss+11; ++s)
    //             {
    //                 n = s; //This ensures extrapolation by signal reversal to any length
    //                 while (n<0 || n>=(int)Lx) { n = (n<0) ? -n-1 : (n<(int)Lx) ? n : 2*(int)Lx-1-n; }
    //                 X += n - prev_n;
    //                 *Xd += *X * D[s-ss];
    //                 prev_n = n;
    //             }
    //         }
    //         else
    //         {
    //             X += ss - prev_n;
    //             for (size_t l=0u; l<11u; ++l, ++X) { *Xd += *X * D[l]; }
    //             X -= 10; prev_n = ss + 1;
    //         }
    //     }
    //     X -= L - 5u;

    //     //Get Xd^2 (1st term of smooth TKEO)
    //     Y += L;
    //     for (size_t l=L; l>0u; --l) { --Xd; --Y; *Y = *Xd * *Xd; }

    //     //Smooth diffdiff of X
    //     ss = -9/2; prev_n = 0;
    //     for (size_t l=L; l>0u; --l, ++ss, ++Xd)
    //     {
    //         *Xd = 0.0;
    //         if (ss<0 || ss>(int)L-9)
    //         {
    //             for (int s=ss; s<ss+9; ++s)
    //             {
    //                 n = s; //This ensures extrapolation by signal reversal to any length
    //                 while (n<0 || n>=(int)L) { n = (n<0) ? -n-1 : (n<(int)L) ? n : 2*(int)L-1-n; }
    //                 X += n - prev_n;
    //                 *Xd += *X * DD[s-ss];
    //                 prev_n = n;
    //             }
    //         }
    //         else
    //         {
    //             X += ss - prev_n;
    //             for (size_t l=0u; l<9; ++l, ++X) { *Xd += *X * DD[l]; }
    //             X -= 8; prev_n = ss + 1;
    //         }
    //     }
    //     X += 4u;

    //     //Get X*Xdd (2nd term of smooth TKEO) and smooth TKEO
    //     Y += L;
    //     for (size_t l=L; l>0u; --l) { *--Y -= *--X * *--Xd; }
    // }
    // else
    // {
    //     const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
    //     const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
    //     const size_t V = N/L, G = V/B;

    //     if (K==1u && (G==1u || B==1u))
    //     {
    //         for (size_t v=V; v>0u; --v)
    //         {
    //             //Smooth diff of X
    //             ss = -11/2; prev_n = 0;
    //             for (size_t l=L; l>0u; --l, ++ss, ++Xd)
    //             {
    //                 *Xd = 0.0;
    //                 if (ss<0 || ss>(int)L-11)
    //                 {
    //                     for (int s=ss; s<ss+11; ++s)
    //                     {
    //                         n = s; //This ensures extrapolation by signal reversal to any length
    //                         while (n<0 || n>=(int)L) { n = (n<0) ? -n-1 : (n<(int)L) ? n : 2*(int)L-1-n; }
    //                         X += n - prev_n;
    //                         *Xd += *X * D[s-ss];
    //                         prev_n = n;
    //                     }
    //                 }
    //                 else
    //                 {
    //                     X += ss - prev_n;
    //                     for (size_t l=0u; l<11u; ++l, ++X) { *Xd += *X * D[l]; }
    //                     X -= 10; prev_n = ss + 1;
    //                 }
    //             }
    //             X -= L - 5u;

    //             //Get Xd^2 (1st term of smooth TKEO)
    //             Y += L;
    //             for (size_t l=L; l>0u; --l) { --Xd; --Y; *Y = *Xd * *Xd; }

    //             //Smooth diffdiff of X
    //             ss = -9/2; prev_n = 0;
    //             for (size_t l=L; l>0u; --l, ++ss, ++Xd)
    //             {
    //                 *Xd = 0.0;
    //                 if (ss<0 || ss>(int)L-9)
    //                 {
    //                     for (int s=ss; s<ss+9; ++s)
    //                     {
    //                         n = s; //This ensures extrapolation by signal reversal to any length
    //                         while (n<0 || n>=(int)L) { n = (n<0) ? -n-1 : (n<(int)L) ? n : 2*(int)L-1-n; }
    //                         X += n - prev_n;
    //                         *Xd += *X * DD[s-ss];
    //                         prev_n = n;
    //                     }
    //                 }
    //                 else
    //                 {
    //                     X += ss - prev_n;
    //                     for (size_t l=0u; l<9; ++l, ++X) { *Xd += *X * DD[l]; }
    //                     X -= 8; prev_n = ss + 1;
    //                 }
    //             }
    //             X -= L - 4u; Xd -= L;

    //             //Get X*Xdd (2nd term of smooth TKEO) and subtrac to get smooth TKEO
    //             for (size_t l=L; l>0u; --l, ++X, ++Xd, ++Y) { *Y -= *X * *Xd; }
    //             Xd -= L;
    //         }
    //     }
    //     else
    //     {
    //         for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
    //         {
    //             for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
    //             {
    //                 //Smooth diff of X
    //                 ss = -11/2; prev_n = 0;
    //                 for (size_t l=L; l>0u; --l, ++ss, ++Xd)
    //                 {
    //                     *Xd = 0.0;
    //                     if (ss<0 || ss>(int)L-11)
    //                     {
    //                         for (int s=ss; s<ss+11; ++s)
    //                         {
    //                             n = s; //This ensures extrapolation by signal reversal to any length
    //                             while (n<0 || n>=(int)L) { n = (n<0) ? -n-1 : (n<(int)L) ? n : 2*(int)L-1-n; }
    //                             X += (int)K * (n-prev_n);
    //                             *Xd += *X * D[s-ss];
    //                             prev_n = n;
    //                         }
    //                     }
    //                     else
    //                     {
    //                         X += (int)K * (ss-prev_n);
    //                         for (size_t l=0u; l<11u; ++l, X+=K) { *Xd += *X * D[l]; }
    //                         X -= 10u*K; prev_n = ss + 1;
    //                     }
    //                 }
    //                 X -= K * (L-5u);

    //                 //Get Xd^2 (1st term of smooth TKEO)
    //                 Y += K*L;
    //                 for (size_t l=L; l>0u; --l) { --Xd; Y-=K; *Y = *Xd * *Xd; }

    //                 //Smooth diffdiff of X
    //                 ss = -9/2; prev_n = 0;
    //                 for (size_t l=L; l>0u; --l, ++ss, ++Xd)
    //                 {
    //                     *Xd = 0.0;
    //                     if (ss<0 || ss>(int)L-9)
    //                     {
    //                         for (int s=ss; s<ss+9; ++s)
    //                         {
    //                             n = s; //This ensures extrapolation by signal reversal to any length
    //                             while (n<0 || n>=(int)L) { n = (n<0) ? -n-1 : (n<(int)L) ? n : 2*(int)L-1-n; }
    //                             X += (int)K * (n-prev_n);
    //                             *Xd += *X * DD[s-ss];
    //                             prev_n = n;
    //                         }
    //                     }
    //                     else
    //                     {
    //                         X += (int)K * (ss-prev_n);
    //                         for (size_t l=0u; l<9; ++l, X+=K) { *Xd += *X * DD[l]; }
    //                         X -= 8u*K; prev_n = ss + 1;
    //                     }
    //                 }
    //                 X -= K*(L-4u); Xd -= L;

    //                 //Get X*Xdd (2nd term of smooth TKEO) and subtrac to get smooth TKEO
    //                 for (size_t l=L; l>0u; --l, X+=K, ++Xd, Y+=K) { *Y -= *X * *Xd; }
    //                 Xd -= L;
    //             }
    //         }
    //     }
    // }
    free(Xd);

    return 0;
}


int tkeo_smooth_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t stp)
{
    if (dim>3u) { fprintf(stderr,"error in tkeo_smooth_c: dim must be in [0 3]\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in tkeo_smooth_c: stp must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = 1u + (Lx-1u)/stp;

    //Init Yr for intermediate output for real part
    float *Yr  = (float *)malloc(Ly*sizeof(float));
    if (!Yr) { fprintf(stderr,"error in tkeo_smooth_c: problem with malloc. "); perror("malloc"); return 1; }

    //Smooth TKEO of real part
    if (!tkeo_smooth_s(Yr,X,R,C,S,H,iscolmajor,dim,2u*stp))
    { fprintf(stderr,"error in tkeo_smooth_c: problem with tkeo_smooth_s for real part\n"); return 1; }

    //Smooth TKEO of imag part
    if (!tkeo_smooth_s(Y,&X[1],R,C,S,H,iscolmajor,dim,2u*stp))
    { fprintf(stderr,"error in tkeo_smooth_c: problem with tkeo_smooth_s for imag part\n"); return 1; }

    //Sum
    for (size_t n=N; n>0u; --n, ++Yr, ++Y) { *Y += *Yr; }
    Yr -= N;

    //Free
    free(Yr);

    return 0;
}


int tkeo_smooth_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t stp)
{
    if (dim>3u) { fprintf(stderr,"error in tkeo_smooth_z: dim must be in [0 3]\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in tkeo_smooth_z: stp must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = 1u + (Lx-1u)/stp;

    //Init Yr for intermediate output for real part
    double *Yr  = (double *)malloc(Ly*sizeof(double));
    if (!Yr) { fprintf(stderr,"error in tkeo_smooth_z: problem with malloc. "); perror("malloc"); return 1; }

    //Smooth TKEO of real part
    if (!tkeo_smooth_d(Yr,X,R,C,S,H,iscolmajor,dim,2u*stp))
    { fprintf(stderr,"error in tkeo_smooth_z: problem with tkeo_smooth_d for real part\n"); return 1; }

    //Smooth TKEO of imag part
    if (!tkeo_smooth_d(Y,&X[1],R,C,S,H,iscolmajor,dim,2u*stp))
    { fprintf(stderr,"error in tkeo_smooth_z: problem with tkeo_smooth_d for imag part\n"); return 1; }

    //Sum
    for (size_t n=N; n>0u; --n, ++Yr, ++Y) { *Y += *Yr; }
    Yr -= N;

    //Free
    free(Yr);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
