//The "unbiased" version uses N-l in the denominator instead of N.
//It is actually just "less biased", but is slower,
//has larger mean-squared error, and doesn't match FFT estimate.

//Recall that FFT is order N*log(N), instead of N*N for sig2ac.
//Since 2 FFTs (and set up time), this means use this if 2*log(N) < N.
//However, this is the case for N>2, so only set up time would prevent using this.
//On quick local-system test, this is worth using for N>?.

#include <stdio.h>
#include <fftw3.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int sig2ac_fft_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int L, const char unbiased);
int sig2ac_fft_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int L, const char unbiased);
int sig2ac_fft_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int L, const char unbiased);
int sig2ac_fft_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int L, const char unbiased);


int sig2ac_fft_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int L, const char unbiased)
{
    if (dim>3u) { fprintf(stderr,"error in sig2ac_fft_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L>Lx) { fprintf(stderr,"error in sig2ac_fft_s: L (num lags in output) must be <= Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        //Get nfft
        const size_t nfft = Lx + L;
        if (nfft>16384u) { nfft += nfft%2u; }
        else { size_t f = 1u; while (f<nfft) { f *= 2u; } nfft = fu; }
        const size_t F = nfft/2u + 1u;

        //Initialize fftw
        double *X1, *Y1;
        fftw_plan fplan, iplan;
        X1 = fftwf_alloc_real(nfft);
        Y1 = fftwf_alloc_real(nfft);
        fplan = fftwf_plan_r2r_1d(nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
        iplan = fftwf_plan_r2r_1d(nfft,Y1,X1,FFTW_R2HC,FFTW_ESTIMATE);
        if (!fplan || !iplan) { fprintf(stderr,"error in sig2ac_fft_s: problem creating fftw plan"); return 1; }
        for (size_t nf=0u; nf<nfft; ++nf) { X1[nf] = 0.0f; }

        if (Lx==N)
        {
            cblas_dcopy((int)R,X,1,X1,1);
            fftw_execute(fplan);
            for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
            for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
            fftw_execute(iplan);
            cblas_dcopy(L,X1,1,Y,1);

            if (corr)
            {
                const float y0 = *Y;
                *Y++ = 1.0f;
                for (size_t l=1u; l<L; ++l, ++Y) { *Y /= y0; }
            }
            else if (unbiased)
            {
                const float den = (float)(nfft*(Lx-l));
                for (size_t l=0u; l<L; ++l, ++Y) { *Y /= den; }
            }
            else //xcorr leaves this blank
            {
                const float den = (float)(nfft*Lx);
                for (size_t l=0u; l<L; ++l, ++Y) { *Y /= den; }
            }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;
            float sm, y0;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, X+=L-1u)
                {
                    for (size_t l=0u; l<L-1u; ++l, X-=Lx-l+1u, ++Y)
                    {
                        sm = 0.0f;
                        for (size_t n=0u; n<Lx-l; ++n, ++X) { sm += *X * *(X+l); }
                        *Y = sm;
                    }
                    sm = 0.0f;
                    for (size_t n=0u; n<Lx-L+1u; ++n, ++X) { sm += *X * *(X+L-1u); }
                    *Y = sm; Y -= L-1u;

                    if (corr)
                    {
                        y0 = *Y; *Y++ = 1.0f;
                        for (size_t l=1u; l<L; ++l, ++Y) { *Y /= y0; }
                    }
                    else if (unbiased)
                    {
                        const float den = (float)(nfft*(Lx-l));
                        for (size_t l=0u; l<L; ++l, ++Y) { *Y /= den; }
                    }
                    else
                    {
                        const float den = (float)(nfft*Lx);
                        for (size_t l=0u; l<L; ++l, ++Y) { *Y /= den; }
                    }
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, ++X, Y-=K*L-1u)
                    {
                        for (size_t l=0u; l<L; ++l, X-=K*(Lx-l+1u), Y+=K)
                        {
                            sm = 0.0f;
                            for (size_t n=0u; n<Lx-l; ++n, X+=K) { sm += *X * *(X+l*K); }
                            *Y = sm;
                        }
                        Y -= K*L;

                        if (corr)
                        {
                            y0 = *Y; *Y = 1.0f; Y += K;
                            for (size_t l=1u; l<L; ++l, Y+=K) { *Y /= y0; }
                        }
                        else if (unbiased)
                        {
                            const float den = (float)(nfft*(Lx-l));
                            for (size_t l=0u; l<L; ++l, Y+=K) { *Y /= den; }
                        }
                        else
                        {
                            const float den = (float)(nfft*Lx);
                            for (size_t l=0u; l<L; ++l, Y+=K) { *Y /= den; }
                        }
                    }
                }
            }
        }
        fftw_destroy_plan(fplan); fftw_destroy_plan(iplan);
        fftw_free(X1); fftw_free(Y1);
    }

    return 0;
}


// int sig2ac_fft_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int L, const char unbiased)
// {
//     if (dim>3) { fprintf(stderr,"error in sig2ac_fft_d: dim must be in [0 3]\n"); return 1; }

//     const double z = 0.0;
//     int r, c, l, f, nfft, F;
//     double *X1, *Y1;
//     fftw_plan fplan, iplan;

//     //Checks
//     if (R<1) { fprintf(stderr,"error in sig2ac_fft_d: nrows X must be positive\n"); return 1; }
//     if (C<1) { fprintf(stderr,"error in sig2ac_fft_d: ncols X must be positive\n"); return 1; }
//     if (dim==0 && L>R) { fprintf(stderr,"error in sig2ac_fft_d: nlags must be < nrows X for dim==0\n"); return 1; }
//     if (dim==1 && L>C) { fprintf(stderr,"error in sig2ac_fft_d: nlags must be < ncols X for dim==1\n"); return 1; }

//     //Get nfft
//     nfft = (dim==0) ? R+L : C+L;
//     if (nfft>16384) { nfft += nfft%2; }
//     else { f = 1; while (f<nfft) { f *= 2; } nfft = f; }
//     F = nfft/2 + 1;

//     //Initialize fftw
//     X1 = fftw_alloc_real((size_t)nfft);
//     Y1 = fftw_alloc_real((size_t)nfft);
//     fplan = fftw_plan_r2r_1d(nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
//     iplan = fftw_plan_r2r_1d(nfft,Y1,X1,FFTW_R2HC,FFTW_ESTIMATE);
//     if (!fplan || !iplan) { fprintf(stderr,"error in sig2ac_fft_d: problem creating fftw plan"); return 1; }

//     if (dim==0)
//     {
//         if (iscolmajor)
//         {
//             for (size_t c=0; c<C; c++)
//             {
//                 cblas_dcopy(nfft-R,&z,0,&X1[R],1); //zero-pad
//                 cblas_dcopy((int)R,&X[c*R],1,&X1[0],1);
//                 //m = cblas_ddot((int)R,&X1[0],1,&o,0) / R;
//                 //cblas_daxpy((int)R,-m,&o,0,&X1[0],1);
//                 fftw_execute(fplan);
//                 for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
//                 for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
//                 fftw_execute(iplan);
//                 cblas_dcopy(L,&X1[0],1,&Y[c*L],1);
//             }
//             cblas_dscal(L*C,1.0/nfft,Y,1);
//             if (unbiased) { for (l=1; l<L; l++) { cblas_dscal((int)C,(double)R/(double)(R-l),&Y[l],L); } }
//         }
//         else
//         {
//             for (size_t c=0; c<C; c++)
//             {
//                 cblas_dcopy(nfft-R,&z,0,&X1[R],1); //zero-pad
//                 cblas_dcopy((int)R,&X[c],(int)C,&X1[0],1);
//                 //m = cblas_ddot((int)R,&X1[0],1,&o,0) / R;
//                 //cblas_daxpy((int)R,-m,&o,0,&X1[0],1);
//                 fftw_execute(fplan);
//                 for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
//                 for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
//                 fftw_execute(iplan);
//                 cblas_dcopy(L,&X1[0],1,&Y[c],1);
//             }
//             cblas_dscal(L*C,1.0/nfft,Y,1);
//             if (unbiased) { for (l=1; l<L; l++) { cblas_dscal((int)C,(double)R/(double)(R-l),&Y[l*C],1); } }
//         }
//     }
//     else if (dim==1)
//     {
//         if (iscolmajor)
//         {
//             for (size_t r=0; r<R; r++)
//             {
//                 cblas_dcopy(nfft-C,&z,0,&X1[C],1); //zero-pad
//                 cblas_dcopy((int)C,&X[r],(int)R,&X1[0],1);
//                 //m = cblas_ddot((int)C,&X1[0],1,&o,0) / C;
//                 //cblas_daxpy((int)C,-m,&o,0,&X1[0],1);
//                 fftw_execute(fplan);
//                 for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
//                 for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
//                 fftw_execute(iplan);
//                 cblas_dcopy(L,&X1[0],1,&Y[r],(int)R);
//             }
//             cblas_dscal(L*R,1.0/nfft,Y,1);
//             if (unbiased) { for (l=1; l<L; l++) { cblas_dscal((int)R,(double)C/(double)(C-l),&Y[l*R],1); } }
//         }
//         else
//         {
//             for (size_t r=0; r<R; r++)
//             {
//                 cblas_dcopy(nfft-C,&z,0,&X1[C],1); //zero-pad
//                 cblas_dcopy((int)C,&X[r*C],1,&X1[0],1);
//                 //m = cblas_ddot((int)C,&X1[0],1,&o,0) / C;
//                 //cblas_daxpy((int)C,-m,&o,0,&X1[0],1);
//                 fftw_execute(fplan);
//                 for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
//                 for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
//                 fftw_execute(iplan);
//                 cblas_dcopy(L,&X1[0],1,&Y[r*L],1);
//             }
//             cblas_dscal(L*R,1.0/nfft,Y,1);
//             if (unbiased) { for (l=1; l<L; l++) { cblas_dscal((int)R,(double)C/(double)(C-l),&Y[l],L); } }
//         }
//     }
//     else
//     {
//         fprintf(stderr,"error in sig2ac_fft_d: dim must be 0 or 1.\n"); return 1;
//     }

//     fftw_destroy_plan(fplan); fftw_destroy_plan(iplan); fftw_free(X1); fftw_free(Y1);
//     return 0;
// }


// int sig2ac_fft_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int L, const char unbiased)
// {
//     if (dim>3) { fprintf(stderr,"error in sig2ac_fft_c: dim must be in [0 3]\n"); return 1; }

//     const float z[2] = {0.0f,0.0f}; //, o[2] = {1.0f,0.0f};
//     //_Complex float m;
//     int r, c, l, f, nfft;
//     fftwf_complex *X1, *Y1;
//     fftwf_plan fplan, iplan;

//     //Checks
//     if (R<1) { fprintf(stderr,"error in sig2ac_fft_c: nrows X must be positive\n"); return 1; }
//     if (C<1) { fprintf(stderr,"error in sig2ac_fft_c: ncols X must be positive\n"); return 1; }
//     if (dim==0 && L>R) { fprintf(stderr,"error in sig2ac_fft_c: nlags must be < nrows X for dim==0\n"); return 1; }
//     if (dim==1 && L>C) { fprintf(stderr,"error in sig2ac_fft_c: nlags must be < ncols X for dim==1\n"); return 1; }

//     //Get nfft
//     nfft = (dim==0) ? R+L : C+L;
//     if (nfft>16384) { nfft += nfft%2; }
//     else { f = 1; while (f<nfft) { f *= 2; } nfft = f; }

//     //Initialize fftw
//     X1 = fftwf_alloc_complex((size_t)nfft);
//     Y1 = fftwf_alloc_complex((size_t)nfft);
//     fplan = fftwf_plan_dft_1d(nfft,X1,Y1,FFTW_FORWARD,FFTW_ESTIMATE);
//     iplan = fftwf_plan_dft_1d(nfft,Y1,X1,FFTW_BACKWARD,FFTW_ESTIMATE);
//     if (!fplan || !iplan) { fprintf(stderr,"error in sig2ac_fft_c: problem creating fftw plan"); return 1; }

//     if (dim==0)
//     {
//         if (iscolmajor)
//         {
//             for (size_t c=0; c<C; c++)
//             {
//                 cblas_ccopy(nfft-R,&z[0],0,(float *)&X1[R],1); //zero-pad
//                 cblas_ccopy((int)R,&X[2*c*R],1,(float *)&X1[0],1);
//                 //m = -cblas_cdotu((int)R,(float *)&X1[0],1,&o[0],0) / R;
//                 //cblas_caxpy((int)R,(float *)&m,&o[0],0,(float *)&X1[0],1);
//                 fftwf_execute(fplan);
//                 #ifdef __cplusplus
//                     for (f=0; f<nfft; f++) { Y1[f][0] = Y1[f][0]*Y1[f][0] + Y1[f][1]*Y1[f][1]; Y1[f][1] = 0.0f; }
//                 #else
//                     for (f=0; f<nfft; f++) { Y1[f] *= conjf(Y1[f]); }
//                 #endif
//                 fftwf_execute(iplan);
//                 cblas_ccopy(L,(float *)&X1[0],1,&Y[2*c*L],1);
//             }
//             cblas_csscal(L*C,1.0f/nfft,Y,1);
//             if (unbiased) { for (l=1; l<L; l++) { cblas_csscal((int)C,(float)R/(float)(R-l),&Y[2*l],L); } }
//         }
//         else
//         {
//             for (size_t c=0; c<C; c++)
//             {
//                 cblas_ccopy(nfft-R,&z[0],0,(float *)&X1[R],1); //zero-pad
//                 cblas_ccopy((int)R,&X[2*c],2*C,(float *)&X1[0],1);
//                 //m = -cblas_cdotu((int)R,(float *)&X1[0],1,&o[0],0) / R;
//                 //cblas_caxpy((int)R,(float *)&m,&o[0],0,(float *)&X1[0],1);
//                 fftwf_execute(fplan);
//                 #ifdef __cplusplus
//                     for (f=0; f<nfft; f++) { Y1[f][0] = Y1[f][0]*Y1[f][0] + Y1[f][1]*Y1[f][1]; Y1[f][1] = 0.0f; }
//                 #else
//                     for (f=0; f<nfft; f++) { Y1[f] *= conjf(Y1[f]); }
//                 #endif
//                 fftwf_execute(iplan);
//                 cblas_ccopy(L,(float *)&X1[0],1,&Y[2*c],(int)C);
//             }
//             cblas_csscal(L*C,1.0f/nfft,Y,1);
//             if (unbiased) { for (l=1; l<L; l++) { cblas_csscal((int)C,(float)R/(float)(R-l),&Y[2*l*C],1); } }
//         }
//     }
//     else if (dim==1)
//     {
//         if (iscolmajor)
//         {
//             for (size_t r=0; r<R; r++)
//             {
//                 cblas_ccopy(nfft-C,&z[0],0,(float *)&X1[C],1); //zero-pad
//                 cblas_ccopy((int)C,&X[2*r],2*R,(float *)&X1[0],1);
//                 //m = -cblas_cdotu((int)C,(float *)&X1[0],1,&o[0],0) / C;
//                 //cblas_caxpy((int)C,(float *)&m,&o[0],0,(float *)&X1[0],1);
//                 fftwf_execute(fplan);
//                 #ifdef __cplusplus
//                     for (f=0; f<nfft; f++) { Y1[f][0] = Y1[f][0]*Y1[f][0] + Y1[f][1]*Y1[f][1]; Y1[f][1] = 0.0f; }
//                 #else
//                     for (f=0; f<nfft; f++) { Y1[f] *= conjf(Y1[f]); }
//                 #endif
//                 fftwf_execute(iplan);
//                 cblas_ccopy(L,(float *)&X1[0],1,&Y[2*r],(int)R);
//             }
//             cblas_csscal(L*R,1.0f/nfft,Y,1);
//             if (unbiased) { for (l=1; l<L; l++) { cblas_csscal((int)R,(float)C/(float)(C-l),&Y[2*l*R],1); } }
//         }
//         else
//         {
//             for (size_t r=0; r<R; r++)
//             {
//                 cblas_ccopy(nfft-C,&z[0],0,(float *)&X1[C],1); //zero-pad
//                 cblas_ccopy((int)C,&X[2*r*C],1,(float *)&X1[0],1);
//                 //m = -cblas_cdotu((int)C,(float *)&X1[0],1,&o[0],0) / C;
//                 //cblas_caxpy((int)C,(float *)&m,&o[0],0,(float *)&X1[0],1);
//                 fftwf_execute(fplan);
//                 #ifdef __cplusplus
//                     for (f=0; f<nfft; f++) { Y1[f][0] = Y1[f][0]*Y1[f][0] + Y1[f][1]*Y1[f][1]; Y1[f][1] = 0.0f; }
//                 #else
//                     for (f=0; f<nfft; f++) { Y1[f] *= conjf(Y1[f]); }
//                 #endif
//                 fftwf_execute(iplan);
//                 cblas_ccopy(L,(float *)&X1[0],1,&Y[2*r*L],1);
//             }
//             cblas_csscal(L*R,1.0f/nfft,Y,1);
//             if (unbiased) { for (l=1; l<L; l++) { cblas_csscal((int)R,(float)C/(float)(C-l),&Y[2*l],L); } }
//         }
//     }
//     else
//     {
//         fprintf(stderr,"error in sig2ac_fft_c: dim must be 0 or 1.\n"); return 1;
//     }

//     fftwf_destroy_plan(fplan); fftwf_destroy_plan(iplan); fftwf_free(X1); fftwf_free(Y1);
//     return 0;
// }


// int sig2ac_fft_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int L, const char unbiased)
// {
//     if (dim>3) { fprintf(stderr,"error in sig2ac_fft_z: dim must be in [0 3]\n"); return 1; }

//     const double z[2] = {0.0,0.0};
//     int r, c, l, f, nfft;
//     fftw_complex *X1, *Y1;
//     fftw_plan fplan, iplan;

//     //Checks
//     if (R<1) { fprintf(stderr,"error in sig2ac_fft_z: nrows X must be positive\n"); return 1; }
//     if (C<1) { fprintf(stderr,"error in sig2ac_fft_z: ncols X must be positive\n"); return 1; }
//     if (dim==0 && L>R) { fprintf(stderr,"error in sig2ac_fft_z: nlags must be < nrows X for dim==0\n"); return 1; }
//     if (dim==1 && L>C) { fprintf(stderr,"error in sig2ac_fft_z: nlags must be < ncols X for dim==1\n"); return 1; }

//     //Get nfft
//     nfft = (dim==0) ? R+L : C+L;
//     if (nfft>16384) { nfft += nfft%2; }
//     else { f = 1; while (f<nfft) { f *= 2; } nfft = f; }

//     //Initialize fftw
//     X1 = fftw_alloc_complex((size_t)nfft);
//     Y1 = fftw_alloc_complex((size_t)nfft);
//     fplan = fftw_plan_dft_1d(nfft,X1,Y1,FFTW_FORWARD,FFTW_ESTIMATE);
//     iplan = fftw_plan_dft_1d(nfft,Y1,X1,FFTW_BACKWARD,FFTW_ESTIMATE);
//     if (!fplan || !iplan) { fprintf(stderr,"error in sig2ac_fft_z: problem creating fftw plan"); return 1; }

//     if (dim==0)
//     {
//         if (iscolmajor)
//         {
//             for (size_t c=0; c<C; c++)
//             {
//                 cblas_zcopy(nfft-R,&z[0],0,(double *)&X1[R],1); //zero-pad
//                 cblas_zcopy((int)R,&X[2*c*R],1,(double *)&X1[0],1);
//                 //m = -cblas_zdotu((int)R,(double *)&X1[0],1,&o[0],0) / R;
//                 //cblas_zaxpy((int)R,(double *)&m,&o[0],0,(double *)&X1[0],1);
//                 fftw_execute(fplan);
//                 #ifdef __cplusplus
//                     for (f=0; f<nfft; f++) { Y1[f][0] = Y1[f][0]*Y1[f][0] + Y1[f][1]*Y1[f][1]; Y1[f][1] = 0.0; }
//                 #else
//                     for (f=0; f<nfft; f++) { Y1[f] *= conj(Y1[f]); }
//                 #endif
//                 fftw_execute(iplan);
//                 cblas_zcopy(L,(double *)&X1[0],1,&Y[2*c*L],1);
//             }
//             cblas_zdscal(L*C,1.0/nfft,Y,1);
//             if (unbiased) { for (l=1; l<L; l++) { cblas_zdscal((int)C,(double)R/(double)(R-l),&Y[2*l],L); } }
//         }
//         else
//         {
//             for (size_t c=0; c<C; c++)
//             {
//                 cblas_zcopy(nfft-R,&z[0],0,(double *)&X1[R],1); //zero-pad
//                 cblas_zcopy((int)R,&X[2*c],2*C,(double *)&X1[0],1);
//                 //m = -cblas_zdotu((int)R,(double *)&X1[0],1,&o[0],0) / R;
//                 //cblas_zaxpy((int)R,(double *)&m,&o[0],0,(double *)&X1[0],1);
//                 fftw_execute(fplan);
//                 #ifdef __cplusplus
//                     for (f=0; f<nfft; f++) { Y1[f][0] = Y1[f][0]*Y1[f][0] + Y1[f][1]*Y1[f][1]; Y1[f][1] = 0.0; }
//                 #else
//                     for (f=0; f<nfft; f++) { Y1[f] *= conj(Y1[f]); }
//                 #endif
//                 fftw_execute(iplan);
//                 cblas_zcopy(L,(double *)&X1[0],1,&Y[2*c],(int)C);
//             }
//             cblas_zdscal(L*C,1.0/nfft,Y,1);
//             if (unbiased) { for (l=1; l<L; l++) { cblas_zdscal((int)C,(double)R/(double)(R-l),&Y[2*l*C],1); } }
//         }
//     }
//     else if (dim==1)
//     {
//         if (iscolmajor)
//         {
//             for (size_t r=0; r<R; r++)
//             {
//                 cblas_zcopy(nfft-C,&z[0],0,(double *)&X1[C],1); //zero-pad
//                 cblas_zcopy((int)C,&X[2*r],2*R,(double *)&X1[0],1);
//                 //m = -cblas_zdotu((int)C,(double *)&X1[0],1,&o[0],0) / C;
//                 //cblas_zaxpy((int)C,(double *)&m,&o[0],0,(double *)&X1[0],1);
//                 fftw_execute(fplan);
//                 #ifdef __cplusplus
//                     for (f=0; f<nfft; f++) { Y1[f][0] = Y1[f][0]*Y1[f][0] + Y1[f][1]*Y1[f][1]; Y1[f][1] = 0.0; }
//                 #else
//                     for (f=0; f<nfft; f++) { Y1[f] *= conj(Y1[f]); }
//                 #endif
//                 fftw_execute(iplan);
//                 cblas_zcopy(L,(double *)&X1[0],1,&Y[2*r],(int)R);
//             }
//             cblas_zdscal(L*R,1.0/nfft,Y,1);
//             if (unbiased) { for (l=1; l<L; l++) { cblas_zdscal((int)R,(double)C/(double)(C-l),&Y[2*l*R],1); } }
//         }
//         else
//         {
//             for (size_t r=0; r<R; r++)
//             {
//                 cblas_zcopy(nfft-C,&z[0],0,(double *)&X1[C],1); //zero-pad
//                 cblas_zcopy((int)C,&X[2*r*C],1,(double *)&X1[0],1);
//                 //m = -cblas_zdotu((int)C,(double *)&X1[0],1,&o[0],0) / C;
//                 //cblas_zaxpy((int)C,(double *)&m,&o[0],0,(double *)&X1[0],1);
//                 fftw_execute(fplan);
//                 #ifdef __cplusplus
//                     for (f=0; f<nfft; f++) { Y1[f][0] = Y1[f][0]*Y1[f][0] + Y1[f][1]*Y1[f][1]; Y1[f][1] = 0.0; }
//                 #else
//                     for (f=0; f<nfft; f++) { Y1[f] *= conj(Y1[f]); }
//                 #endif
//                 fftw_execute(iplan);
//                 cblas_zcopy(L,(double *)&X1[0],1,&Y[2*r*L],1);
//             }
//             cblas_zdscal(L*R,1.0/nfft,Y,1);
//             if (unbiased) { for (l=1; l<L; l++) { cblas_zdscal((int)R,(double)C/(double)(C-l),&Y[2*l],L); } }
//         }
//     }
//     else
//     {
//         fprintf(stderr,"error in sig2ac_fft_z: dim must be 0 or 1.\n"); return 1;
//     }

//     fftw_destroy_plan(fplan); fftw_destroy_plan(iplan); fftw_free(X1); fftw_free(Y1);
//     return 0;
// }


#ifdef __cplusplus
}
}
#endif
