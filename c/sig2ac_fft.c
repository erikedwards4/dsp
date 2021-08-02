//The "unbiased" version uses N-l in the denominator instead of N.
//It is actually just "less biased", but is slower,
//has larger mean-squared error, and doesn't match FFT estimate.

//Recall that FFT is order N*log(N), instead of N*N for sig2ac.
//Since 2 FFTs (and set up time), this means use this if 2*log(N) < N.
//However, this is the case for N>2, so only set up time would prevent using this.
//On quick local-system test, this is worth using for N>?.

#include <stdio.h>
#include <fftw3.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int sig2ac_fft_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t L, const int unbiased, const int corr);
int sig2ac_fft_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t L, const int unbiased, const int corr);
int sig2ac_fft_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t L, const int unbiased, const int corr);
int sig2ac_fft_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t L, const int unbiased, const int corr);


int sig2ac_fft_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t L, const int unbiased, const int corr)
{
    if (dim>3u) { fprintf(stderr,"error in sig2ac_fft_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L>Lx) { fprintf(stderr,"error in sig2ac_fft_s: L (num lags in output) must be <= Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        //Get nfft
        size_t nfft = Lx + L;
        if (nfft>16384u) { nfft += nfft%2u; }
        else { size_t f = 1u; while (f<nfft) { f *= 2u; } nfft = f; }
        const size_t F = nfft/2u + 1u;

        //Initialize fftw
        float *X1, *Y1;
        fftwf_plan fplan, iplan;
        X1 = fftwf_alloc_real(nfft);
        Y1 = fftwf_alloc_real(nfft);
        fplan = fftwf_plan_r2r_1d((int)nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
        iplan = fftwf_plan_r2r_1d((int)nfft,Y1,X1,FFTW_R2HC,FFTW_ESTIMATE);
        if (!fplan || !iplan) { fprintf(stderr,"error in sig2ac_fft_s: problem creating fftw plan"); return 1; }
        for (size_t nf=0u; nf<nfft; ++nf) { X1[nf] = 0.0f; }

        if (Lx==N)
        {
            for (size_t l=0u; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
            X1 -= Lx;
            fftwf_execute(fplan);
            
            //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

            //Fastest
            for (size_t f=0u; f<nfft; ++f, ++Y1) { *Y1 *= *Y1; }
            Y1 -= nfft - 1u;
            for (size_t f=1u; f<F-1u; ++f, ++Y1) { *Y1 += *(Y1+nfft-2u*f); *(Y1+nfft-2u*f) = *Y1; }
            Y1 -= F - 1u;

            //2nd fastest
            // for (size_t f=0u; f<nfft; ++f, ++Y1) { *Y1 *= *Y1; }
            // --Y1;
            // for (size_t f=1u; f<F-1u; ++f, --Y1) { *Y1 += *(Y1-nfft+2u*f); *(Y1-nfft+2u*f) = *Y1; }
            // Y1 -= nfft - F + 1u;

            //3rd fastest
            // for (size_t f=0u; f<nfft; ++f) { Y1[f] *= Y1[f]; }
            // for (size_t f=1u; f<F-1u; ++f) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }

            //Slowest (along with other methods that copy into a 2nd array)
            // for (size_t f=0u; f<F; ++f, ++Y1, ++Y2) { *Y2 = *Y1 * *Y1; }
            // Y2 -= 2u;
            // for (size_t f=1u; f<F-1u; ++f, ++Y1, --Y2) { *Y2 += *Y1 * *Y1; }
            // Y1 -= nfft;

            //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
            
            fftwf_execute(iplan);
            for (size_t l=0u; l<L; ++l, ++X1, ++Y) { *Y = *X1; }
            X1 -= L; Y -= L;

            if (corr)
            {
                const float y0 = *Y;
                *Y++ = 1.0f;
                for (size_t l=1u; l<L; ++l, ++Y) { *Y /= y0; }
            }
            else if (unbiased)
            {
                for (size_t l=0u; l<L; ++l, ++Y) { *Y /= (float)(nfft*(Lx-l)); }
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

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v)
                {
                    for (size_t l=0u; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= Lx;
                    fftwf_execute(fplan);
                    for (size_t f=0u; f<nfft; ++f, ++Y1) { *Y1 *= *Y1; }
                    Y1 -= nfft - 1u;
                    for (size_t f=1u; f<F-1u; ++f, ++Y1) { *Y1 += *(Y1+nfft-2u*f); *(Y1+nfft-2u*f) = *Y1; }
                    Y1 -= F - 1u;
                    fftwf_execute(iplan);
                    for (size_t l=0u; l<L; ++l, ++X1, ++Y) { *Y = *X1; }
                    X1 -= L; Y -= L;

                    if (corr)
                    {
                        const float y0 = *Y;
                        *Y++ = 1.0f;
                        for (size_t l=1u; l<L; ++l, ++Y) { *Y /= y0; }
                    }
                    else if (unbiased)
                    {
                        for (size_t l=0u; l<L; ++l, ++Y) { *Y /= (float)(nfft*(Lx-l)); }
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
                    for (size_t b=0u; b<B; ++b, X-=K*Lx-1u, Y-=K*L-1u)
                    {
                        for (size_t l=0u; l<Lx; ++l, X+=K, ++X1) { *X1 = *X; }
                        X1 -= Lx;
                        fftwf_execute(fplan);
                        for (size_t f=0u; f<nfft; ++f, ++Y1) { *Y1 *= *Y1; }
                        Y1 -= nfft - 1u;
                        for (size_t f=1u; f<F-1u; ++f, ++Y1) { *Y1 += *(Y1+nfft-2u*f); *(Y1+nfft-2u*f) = *Y1; }
                        Y1 -= F - 1u;
                        fftwf_execute(iplan);
                        for (size_t l=0u; l<L; ++l, ++X1, Y+=K) { *Y = *X1; }
                        X1 -= L; Y -= K*L;

                        if (corr)
                        {
                            const float y0 = *Y;
                            *Y = 1.0f; Y += K;
                            for (size_t l=1u; l<L; ++l, Y+=K) { *Y /= y0; }
                        }
                        else if (unbiased)
                        {
                            for (size_t l=0u; l<L; ++l, Y+=K) { *Y /= (float)(nfft*(Lx-l)); }
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
        fftwf_destroy_plan(fplan); fftwf_destroy_plan(iplan);
        fftwf_free(X1); fftwf_free(Y1);
    }

    return 0;
}


int sig2ac_fft_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t L, const int unbiased, const int corr)
{
    if (dim>3u) { fprintf(stderr,"error in sig2ac_fft_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L>Lx) { fprintf(stderr,"error in sig2ac_fft_d: L (num lags in output) must be <= Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        //Get nfft
        size_t nfft = Lx + L;
        if (nfft>16384u) { nfft += nfft%2u; }
        else { size_t f = 1u; while (f<nfft) { f *= 2u; } nfft = f; }
        const size_t F = nfft/2u + 1u;

        //Initialize fftw
        double *X1, *Y1;
        fftw_plan fplan, iplan;
        X1 = fftw_alloc_real(nfft);
        Y1 = fftw_alloc_real(nfft);
        fplan = fftw_plan_r2r_1d((int)nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
        iplan = fftw_plan_r2r_1d((int)nfft,Y1,X1,FFTW_R2HC,FFTW_ESTIMATE);
        if (!fplan || !iplan) { fprintf(stderr,"error in sig2ac_fft_d: problem creating fftw plan"); return 1; }
        for (size_t nf=0u; nf<nfft; ++nf) { X1[nf] = 0.0; }

        if (Lx==N)
        {
            for (size_t l=0u; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
            X1 -= Lx;
            fftw_execute(fplan);
            for (size_t f=0u; f<nfft; ++f, ++Y1) { *Y1 *= *Y1; }
            Y1 -= nfft - 1u;
            for (size_t f=1u; f<F-1u; ++f, ++Y1) { *Y1 += *(Y1+nfft-2u*f); *(Y1+nfft-2u*f) = *Y1; }
            Y1 -= F - 1u;
            fftw_execute(iplan);
            for (size_t l=0u; l<L; ++l, ++X1, ++Y) { *Y = *X1; }
            X1 -= L; Y -= L;

            if (corr)
            {
                const double y0 = *Y;
                *Y++ = 1.0;
                for (size_t l=1u; l<L; ++l, ++Y) { *Y /= y0; }
            }
            else if (unbiased)
            {
                for (size_t l=0u; l<L; ++l, ++Y) { *Y /= (double)(nfft*(Lx-l)); }
            }
            else //xcorr leaves this blank
            {
                const double den = (double)(nfft*Lx);
                for (size_t l=0u; l<L; ++l, ++Y) { *Y /= den; }
            }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v)
                {
                    for (size_t l=0u; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= Lx;
                    fftw_execute(fplan);
                    for (size_t f=0u; f<nfft; ++f, ++Y1) { *Y1 *= *Y1; }
                    Y1 -= nfft - 1u;
                    for (size_t f=1u; f<F-1u; ++f, ++Y1) { *Y1 += *(Y1+nfft-2u*f); *(Y1+nfft-2u*f) = *Y1; }
                    Y1 -= F - 1u;
                    fftw_execute(iplan);
                    for (size_t l=0u; l<L; ++l, ++X1, ++Y) { *Y = *X1; }
                    X1 -= L; Y -= L;

                    if (corr)
                    {
                        const double y0 = *Y;
                        *Y++ = 1.0;
                        for (size_t l=1u; l<L; ++l, ++Y) { *Y /= y0; }
                    }
                    else if (unbiased)
                    {
                        for (size_t l=0u; l<L; ++l, ++Y) { *Y /= (double)(nfft*(Lx-l)); }
                    }
                    else
                    {
                        const double den = (double)(nfft*Lx);
                        for (size_t l=0u; l<L; ++l, ++Y) { *Y /= den; }
                    }
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=K*Lx-1u, Y-=K*L-1u)
                    {
                        for (size_t l=0u; l<Lx; ++l, X+=K, ++X1) { *X1 = *X; }
                        X1 -= Lx;
                        fftw_execute(fplan);
                        for (size_t f=0u; f<nfft; ++f, ++Y1) { *Y1 *= *Y1; }
                        Y1 -= nfft - 1u;
                        for (size_t f=1u; f<F-1u; ++f, ++Y1) { *Y1 += *(Y1+nfft-2u*f); *(Y1+nfft-2u*f) = *Y1; }
                        Y1 -= F - 1u;
                        fftw_execute(iplan);
                        for (size_t l=0u; l<L; ++l, ++X1, Y+=K) { *Y = *X1; }
                        X1 -= L; Y -= K*L;

                        if (corr)
                        {
                            const double y0 = *Y;
                            *Y = 1.0; Y += K;
                            for (size_t l=1u; l<L; ++l, Y+=K) { *Y /= y0; }
                        }
                        else if (unbiased)
                        {
                            for (size_t l=0u; l<L; ++l, Y+=K) { *Y /= (double)(nfft*(Lx-l)); }
                        }
                        else
                        {
                            const double den = (double)(nfft*Lx);
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


int sig2ac_fft_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t L, const int unbiased, const int corr)
{
    if (dim>3u) { fprintf(stderr,"error in sig2ac_fft_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L>Lx) { fprintf(stderr,"error in sig2ac_fft_c: L (num lags in output) must be <= Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        //Get nfft
        size_t nfft = Lx + L;
        if (nfft>16384u) { nfft += nfft%2u; }
        else { size_t f = 1u; while (f<nfft) { f *= 2u; } nfft = f; }

        //Initialize fftw
        float *X1, *Y1, yr, yi;
        X1 = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        Y1 = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        fftwf_plan fplan, iplan;
        fplan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)X1,(fftwf_complex *)Y1,FFTW_FORWARD,FFTW_ESTIMATE);
        iplan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)Y1,(fftwf_complex *)X1,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!fplan || !iplan) { fprintf(stderr,"error in sig2ac_fft_c: problem creating fftw plan"); return 1; }
        for (size_t nf=0u; nf<2u*nfft; ++nf) { X1[nf] = 0.0f; }

        if (Lx==N)
        {
            for (size_t l=0u; l<2u*Lx; ++l, ++X, ++X1) { *X1 = *X; }
            X1 -= 2u*Lx;
            fftwf_execute(fplan);
            for (size_t f=0u; f<nfft; ++f)
            {
                yr = *Y1; yi = *(Y1+1);
                *Y1++ = yr*yr + yi*yi;
                *Y1++ = 0.0f;
            }
            Y1 -= 2u*nfft;
            fftwf_execute(iplan);
            for (size_t l=0u; l<2u*L; ++l, ++X1, ++Y) { *Y = *X1; }
            X1 -= 2u*L; Y -= 2u*L;

            if (corr)
            {
                const float y0r = *Y, y0i = *(Y+1), y0a = y0r*y0r + y0i*y0i;
                *Y++ = 1.0f; *Y++ = 0.0f;
                for (size_t l=1u; l<L; ++l)
                {
                    yr = *Y; yi = *(Y+1);
                    *Y++ = (yr*y0r+yi*y0i) / y0a;
                    *Y++ = (yi*y0r-yr*y0i) / y0a;
                }
            }
            else if (unbiased)
            {
                float den;
                for (size_t l=0u; l<L; ++l) { den = (float)(nfft*(Lx-l)); *Y++ /= den; *Y++ /= den; }
            }
            else
            {
                const float den = (float)(nfft*Lx);
                for (size_t l=0u; l<2u*L; ++l, ++Y) { *Y /= den; }
            }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v)
                {
                    for (size_t l=0u; l<2u*Lx; ++l, ++X, ++X1) { *X1 = *X; }
                    for (size_t l=2u*Lx; l<2u*nfft; ++l, ++X1) { *X1 = 0.0f; }
                    X1 -= 2u*nfft;
                    fftwf_execute(fplan);
                    for (size_t f=0u; f<nfft; ++f)
                    {
                        yr = *Y1; yi = *(Y1+1);
                        *Y1++ = yr*yr + yi*yi;
                        *Y1++ = 0.0f;
                    }
                    Y1 -= 2u*nfft;
                    fftwf_execute(iplan);
                    for (size_t l=0u; l<2u*L; ++l, ++X1, ++Y) { *Y = *X1; }
                    X1 -= 2u*L; Y -= 2u*L;

                    if (corr)
                    {
                        const float y0r = *Y, y0i = *(Y+1), y0a = y0r*y0r + y0i*y0i;
                        *Y++ = 1.0f; *Y++ = 0.0f;
                        for (size_t l=1u; l<L; ++l)
                        {
                            yr = *Y; yi = *(Y+1);
                            *Y++ = (yr*y0r+yi*y0i) / y0a;
                            *Y++ = (yi*y0r-yr*y0i) / y0a;
                        }
                    }
                    else if (unbiased)
                    {
                        float den;
                        for (size_t l=0u; l<L; ++l) { den = (float)(nfft*(Lx-l)); *Y++ /= den; *Y++ /= den; }
                    }
                    else
                    {
                        const float den = (float)(nfft*Lx);
                        for (size_t l=0u; l<2u*L; ++l, ++Y) { *Y /= den; }
                    }
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=2u*K*Lx-2u, Y-=2u*K*L-2u)
                    {
                        for (size_t l=0u; l<Lx; ++l, X+=2u*K) { *X1++ = *X; *X1++ = *(X+1); }
                        for (size_t l=2u*Lx; l<2u*nfft; ++l, ++X1) { *X1 = 0.0f; }
                        X1 -= 2u*nfft;
                        fftwf_execute(fplan);
                        for (size_t f=0u; f<nfft; ++f)
                        {
                            yr = *Y1; yi = *(Y1+1);
                            *Y1++ = yr*yr + yi*yi;
                            *Y1++ = 0.0f;
                        }
                        Y1 -= 2u*nfft;
                        fftwf_execute(iplan);
                        for (size_t l=0u; l<L; ++l, Y+=2u*K) { *Y = *X1++; *(Y+1) = *X1++; }
                        X1 -= 2u*L; Y -= 2u*K*L;

                        if (corr)
                        {
                            const float y0r = *Y, y0i = *(Y+1), y0a = y0r*y0r + y0i*y0i;
                            *Y = 1.0f; *(Y+1) = 0.0f; Y += 2u*K;
                            for (size_t l=1u; l<L; ++l, Y+=2u*K)
                            {
                                yr = *Y; yi = *(Y+1);
                                *Y = (yr*y0r+yi*y0i) / y0a;
                                *(Y+1) = (yi*y0r-yr*y0i) / y0a;
                            }
                        }
                        else if (unbiased)
                        {
                            float den;
                            for (size_t l=0u; l<L; ++l, Y+=2u*K) { den = (float)(nfft*(Lx-l)); *Y /= den; *(Y+1) /= den; }
                        }
                        else
                        {
                            const float den = (float)(nfft*Lx);
                            for (size_t l=0u; l<L; ++l, Y+=2u*K) { *Y /= den; *(Y+1) /= den; }
                        }
                    }
                }
            }
        }
        fftwf_destroy_plan(fplan); fftwf_destroy_plan(iplan);
        fftwf_free(X1); fftwf_free(Y1);
    }

    return 0;
}


int sig2ac_fft_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t L, const int unbiased, const int corr)
{
    if (dim>3u) { fprintf(stderr,"error in sig2ac_fft_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L>Lx) { fprintf(stderr,"error in sig2ac_fft_z: L (num lags in output) must be <= Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        //Get nfft
        size_t nfft = Lx + L;
        if (nfft>16384u) { nfft += nfft%2u; }
        else { size_t f = 1u; while (f<nfft) { f *= 2u; } nfft = f; }

        //Initialize fftw
        double *X1, *Y1, yr, yi;
        X1 = (double *)fftw_malloc(2u*nfft*sizeof(double));
        Y1 = (double *)fftw_malloc(2u*nfft*sizeof(double));
        fftw_plan fplan, iplan;
        fplan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)X1,(fftw_complex *)Y1,FFTW_FORWARD,FFTW_ESTIMATE);
        iplan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)Y1,(fftw_complex *)X1,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!fplan || !iplan) { fprintf(stderr,"error in sig2ac_fft_z: problem creating fftw plan"); return 1; }
        for (size_t nf=0u; nf<2u*nfft; ++nf) { X1[nf] = 0.0; }

        if (Lx==N)
        {
            for (size_t l=0u; l<2u*Lx; ++l, ++X, ++X1) { *X1 = *X; }
            X1 -= 2u*Lx;
            fftw_execute(fplan);
            for (size_t f=0u; f<nfft; ++f)
            {
                yr = *Y1; yi = *(Y1+1);
                *Y1++ = yr*yr + yi*yi;
                *Y1++ = 0.0;
            }
            Y1 -= 2u*nfft;
            fftw_execute(iplan);
            for (size_t l=0u; l<2u*L; ++l, ++X1, ++Y) { *Y = *X1; }
            X1 -= 2u*L; Y -= 2u*L;

            if (corr)
            {
                const double y0r = *Y, y0i = *(Y+1), y0a = y0r*y0r + y0i*y0i;
                *Y++ = 1.0; *Y++ = 0.0;
                for (size_t l=1u; l<L; ++l)
                {
                    yr = *Y; yi = *(Y+1);
                    *Y++ = (yr*y0r+yi*y0i) / y0a;
                    *Y++ = (yi*y0r-yr*y0i) / y0a;
                }
            }
            else if (unbiased)
            {
                double den;
                for (size_t l=0u; l<L; ++l) { den = (double)(nfft*(Lx-l)); *Y++ /= den; *Y++ /= den; }
            }
            else
            {
                const double den = (double)(nfft*Lx);
                for (size_t l=0u; l<2u*L; ++l, ++Y) { *Y /= den; }
            }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v)
                {
                    for (size_t l=0u; l<2u*Lx; ++l, ++X, ++X1) { *X1 = *X; }
                    for (size_t l=2u*Lx; l<2u*nfft; ++l, ++X1) { *X1 = 0.0; }
                    X1 -= 2u*nfft;
                    fftw_execute(fplan);
                    for (size_t f=0u; f<nfft; ++f)
                    {
                        yr = *Y1; yi = *(Y1+1);
                        *Y1++ = yr*yr + yi*yi;
                        *Y1++ = 0.0;
                    }
                    Y1 -= 2u*nfft;
                    fftw_execute(iplan);
                    for (size_t l=0u; l<2u*L; ++l, ++X1, ++Y) { *Y = *X1; }
                    X1 -= 2u*L; Y -= 2u*L;

                    if (corr)
                    {
                        const double y0r = *Y, y0i = *(Y+1), y0a = y0r*y0r + y0i*y0i;
                        *Y++ = 1.0; *Y++ = 0.0;
                        for (size_t l=1u; l<L; ++l)
                        {
                            yr = *Y; yi = *(Y+1);
                            *Y++ = (yr*y0r+yi*y0i) / y0a;
                            *Y++ = (yi*y0r-yr*y0i) / y0a;
                        }
                    }
                    else if (unbiased)
                    {
                        double den;
                        for (size_t l=0u; l<L; ++l) { den = (double)(nfft*(Lx-l)); *Y++ /= den; *Y++ /= den; }
                    }
                    else
                    {
                        const double den = (double)(nfft*Lx);
                        for (size_t l=0u; l<2u*L; ++l, ++Y) { *Y /= den; }
                    }
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=2u*K*Lx-2u, Y-=2u*K*L-2u)
                    {
                        for (size_t l=0u; l<Lx; ++l, X+=2u*K) { *X1++ = *X; *X1++ = *(X+1); }
                        for (size_t l=2u*Lx; l<2u*nfft; ++l, ++X1) { *X1 = 0.0; }
                        X1 -= 2u*nfft;
                        fftw_execute(fplan);
                        for (size_t f=0u; f<nfft; ++f)
                        {
                            yr = *Y1; yi = *(Y1+1);
                            *Y1++ = yr*yr + yi*yi;
                            *Y1++ = 0.0;
                        }
                        Y1 -= 2u*nfft;
                        fftw_execute(iplan);
                        for (size_t l=0u; l<L; ++l, Y+=2u*K) { *Y = *X1++; *(Y+1) = *X1++; }
                        X1 -= 2u*L; Y -= 2u*K*L;

                        if (corr)
                        {
                            const double y0r = *Y, y0i = *(Y+1), y0a = y0r*y0r + y0i*y0i;
                            *Y = 1.0; *(Y+1) = 0.0; Y += 2u*K;
                            for (size_t l=1u; l<L; ++l, Y+=2u*K)
                            {
                                yr = *Y; yi = *(Y+1);
                                *Y = (yr*y0r+yi*y0i) / y0a;
                                *(Y+1) = (yi*y0r-yr*y0i) / y0a;
                            }
                        }
                        else if (unbiased)
                        {
                            double den;
                            for (size_t l=0u; l<L; ++l, Y+=2u*K) { den = (double)(nfft*(Lx-l)); *Y /= den; *(Y+1) /= den; }
                        }
                        else
                        {
                            const double den = (double)(nfft*Lx);
                            for (size_t l=0u; l<L; ++l, Y+=2u*K) { *Y /= den; *(Y+1) /= den; }
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


#ifdef __cplusplus
}
}
#endif
