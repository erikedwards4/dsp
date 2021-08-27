//Gets autocovariance (AC) for lags 0 to Ly-1 for each vector in X.
//By default, the means of X are NOT subtracted,
//and Y is NOT normalized by lag 0, 
//which matches the behavior of Octave's xcorr.

//The mnz option allows the mean to be zeroed first for each vec in X.
//However, this also required the removal of the const qualifier for input *X.

//The "biased" version uses N in the denominator,
//which is unlike Octave's xcorr, which leaves it unnormalized.

//The "unbiased" version uses N-l in the denominator instead of N.
//It is actually just "less biased", but is slower,
//has larger mean-squared error, and doesn't match FFT estimate.

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int sig2ac_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Ly, const int mnz, const int unbiased, const int corr);
int sig2ac_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Ly, const int mnz, const int unbiased, const int corr);
int sig2ac_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Ly, const int mnz, const int unbiased, const int corr);
int sig2ac_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Ly, const int mnz, const int unbiased, const int corr);


int sig2ac_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Ly, const int mnz, const int unbiased, const int corr)
{
    if (dim>3u) { fprintf(stderr,"error in sig2ac_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Ly>Lx) { fprintf(stderr,"error in sig2ac_s: Ly (num lags in output) must be <= Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        if (Lx==N)
        {
            if (mnz)
            {
                float mn = 0.0f;
                for (size_t l=Lx; l>0u; --l, ++X) { mn += *X; }
                mn /= (float)Lx;
                for (size_t l=Lx; l>0u; --l) { *--X -= mn; }
            }

            //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

            if (Ly>850u && Lx>80000u) //I also leave this solution since it is closer to real-time use
            {
                float x;
                for (size_t l=Ly; l>0u; --l, ++Y) { *Y = 0.0f; }
                Y -= Ly;
                for (size_t n=0u; n<N; ++n, ++X)
                {
                    x = *X;
                    for (size_t l=0u; l<Ly && l<N-n; ++l) { Y[l] += x * X[l]; }
                }
            }
            else if (Ly%2u)
            {
                float sm;
                for (size_t l=0u; l<Ly; ++l, X-=Lx-l+1u, ++Y)
                {
                    sm = 0.0f;
                    for (size_t n=Lx-l; n>0u; --n, ++X) { sm += *X * *(X+l); }
                    *Y = sm;
                }
                Y -= Ly;
            }
            else //tested, and this back-and-forth is slightly faster for even Ly only
            {
                float sm;
                for (size_t l=0u; l<Ly; ++l, ++Y)
                {
                    sm = 0.0f;
                    for (size_t n=Lx-l; n>0u; --n, ++X) { sm += *X * *(X+l); }
                    *Y++ = sm;
                    sm = 0.0f; ++l; --X;
                    for (size_t n=Lx-l; n>0u; --n) { --X; sm += *X * *(X+l); }
                    *Y = sm;
                }
                Y -= Ly;
            }

            //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(double)(toc.tv_sec-tic.tv_sec)*1e3+(double)(toc.tv_nsec-tic.tv_nsec)/1e6);

            if (corr)
            {
                const float y0 = *Y;
                *Y++ = 1.0f;
                for (size_t l=Ly; l>1u; --l, ++Y) { *Y /= y0; }
            }
            else if (unbiased)
            {
                for (size_t l=0u; l<Ly; ++l, ++Y) { *Y /= (float)(Lx-l); }
            }
            else //xcorr leaves this blank
            {
                for (size_t l=Ly; l>0u; --l, ++Y) { *Y /= (float)Lx; }
            }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;
            float sm;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=Ly-1u)
                {
                    if (mnz)
                    {
                        float mn = 0.0f;
                        for (size_t l=Lx; l>0u; --l, ++X) { mn += *X; }
                        mn /= (float)Lx;
                        for (size_t l=Lx; l>0u; --l) { *--X -= mn; }
                    }

                    for (size_t l=0u; l<Ly-1u; ++l, X-=Lx-l+1u, ++Y)
                    {
                        sm = 0.0f;
                        for (size_t n=Lx-l; n>0u; --n, ++X) { sm += *X * *(X+l); }
                        *Y = sm;
                    }
                    sm = 0.0f;
                    for (size_t n=0u; n<Lx-Ly+1u; ++n, ++X) { sm += *X * *(X+Ly-1u); }
                    *Y = sm; Y -= Ly-1u;

                    if (corr)
                    {
                        const float y0 = *Y;
                        *Y++ = 1.0f;
                        for (size_t l=Ly; l>1u; --l, ++Y) { *Y /= y0; }
                    }
                    else if (unbiased)
                    {
                        for (size_t l=0u; l<Ly; ++l, ++Y) { *Y /= (float)(Lx-l); }
                    }
                    else
                    {
                        for (size_t l=Ly; l>0u; --l, ++Y) { *Y /= (float)Lx; }
                    }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(Ly-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, Y-=K*Ly-1u)
                    {
                        if (mnz)
                        {
                            float mn = 0.0f;
                            for (size_t l=Lx; l>0u; --l, X+=K) { mn += *X; }
                            mn /= (float)Lx;
                            for (size_t l=Lx; l>0u; --l) { X-=K; *X -= mn; }
                        }

                        for (size_t l=0u; l<Ly; ++l, X-=K*(Lx-l+1u), Y+=K)
                        {
                            sm = 0.0f;
                            for (size_t n=0u; n<Lx-l; ++n, X+=K) { sm += *X * *(X+l*K); }
                            *Y = sm;
                        }
                        Y -= K*Ly;

                        if (corr)
                        {
                            const float y0 = *Y;
                            *Y = 1.0f; Y += K;
                            for (size_t l=Ly; l>1u; --l, Y+=K) { *Y /= y0; }
                        }
                        else if (unbiased)
                        {
                            for (size_t l=0u; l<Ly; ++l, Y+=K) { *Y /= (float)(Lx-l); }
                        }
                        else
                        {
                            for (size_t l=Ly; l>0u; --l, Y+=K) { *Y /= (float)Lx; }
                        }
                    }
                }
            }
        }
    }

    return 0;
}


int sig2ac_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Ly, const int mnz, const int unbiased, const int corr)
{
    if (dim>3u) { fprintf(stderr,"error in sig2ac_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Ly>Lx) { fprintf(stderr,"error in sig2ac_d: Ly (num lags in output) must be <= Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        if (Lx==N)
        {
            if (mnz)
            {
                double mn = 0.0;
                for (size_t l=Lx; l>0u; --l, ++X) { mn += *X; }
                mn /= (double)Lx;
                for (size_t l=Lx; l>0u; --l) { *--X -= mn; }
            }

            if (Ly>850u && Lx>80000u) //I also leave this solution since it is closer to real-time use
            {
                double x;
                for (size_t l=Ly; l>0u; --l, ++Y) { *Y = 0.0; }
                Y -= Ly;
                for (size_t n=0u; n<N; ++n, ++X)
                {
                    x = *X;
                    for (size_t l=0u; l<Ly && l<N-n; ++l) { Y[l] += x * X[l]; }
                }
            }
            else if (Ly%2u)
            {
                double sm;
                for (size_t l=0u; l<Ly; ++l, X-=Lx-l+1u, ++Y)
                {
                    sm = 0.0;
                    for (size_t n=Lx-l; n>0u; --n, ++X) { sm += *X * *(X+l); }
                    *Y = sm;
                }
                Y -= Ly;
            }
            else //tested, and this back-and-forth is slightly faster for even Ly only
            {
                double sm;
                for (size_t l=0u; l<Ly; ++l, ++Y)
                {
                    sm = 0.0;
                    for (size_t n=Lx-l; n>0u; --n, ++X) { sm += *X * *(X+l); }
                    *Y++ = sm;
                    sm = 0.0; ++l; --X;
                    for (size_t n=Lx-l; n>0u; --n) { --X; sm += *X * *(X+l); }
                    *Y = sm;
                }
                Y -= Ly;
            }

            if (corr)
            {
                const double y0 = *Y;
                *Y++ = 1.0;
                for (size_t l=Ly; l>1u; --l, ++Y) { *Y /= y0; }
            }
            else if (unbiased)
            {
                for (size_t l=0u; l<Ly; ++l, ++Y) { *Y /= (double)(Lx-l); }
            }
            else //xcorr leaves this blank
            {
                for (size_t l=Ly; l>0u; --l, ++Y) { *Y /= (double)Lx; }
            }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;
            double sm;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=Ly-1u)
                {
                    if (mnz)
                    {
                        double mn = 0.0;
                        for (size_t l=Lx; l>0u; --l, ++X) { mn += *X; }
                        mn /= (double)Lx;
                        for (size_t l=Lx; l>0u; --l) { *--X -= mn; }
                    }

                    for (size_t l=0u; l<Ly-1u; ++l, X-=Lx-l+1u, ++Y)
                    {
                        sm = 0.0;
                        for (size_t n=Lx-l; n>0u; --n, ++X) { sm += *X * *(X+l); }
                        *Y = sm;
                    }
                    sm = 0.0;
                    for (size_t n=0u; n<Lx-Ly+1u; ++n, ++X) { sm += *X * *(X+Ly-1u); }
                    *Y = sm; Y -= Ly-1u;

                    if (corr)
                    {
                        const double y0 = *Y;
                        *Y++ = 1.0;
                        for (size_t l=Ly; l>1u; --l, ++Y) { *Y /= y0; }
                    }
                    else if (unbiased)
                    {
                        for (size_t l=0u; l<Ly; ++l, ++Y) { *Y /= (double)(Lx-l); }
                    }
                    else
                    {
                        for (size_t l=Ly; l>0u; --l, ++Y) { *Y /= (double)Lx; }
                    }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(Ly-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, Y-=K*Ly-1u)
                    {
                        if (mnz)
                        {
                            double mn = 0.0;
                            for (size_t l=Lx; l>0u; --l, X+=K) { mn += *X; }
                            mn /= (double)Lx;
                            for (size_t l=Lx; l>0u; --l) { X-=K; *X -= mn; }
                        }

                        for (size_t l=0u; l<Ly; ++l, X-=K*(Lx-l+1u), Y+=K)
                        {
                            sm = 0.0;
                            for (size_t n=0u; n<Lx-l; ++n, X+=K) { sm += *X * *(X+l*K); }
                            *Y = sm;
                        }
                        Y -= K*Ly;

                        if (corr)
                        {
                            const double y0 = *Y;
                            *Y = 1.0; Y += K;
                            for (size_t l=Ly; l>1u; --l, Y+=K) { *Y /= y0; }
                        }
                        else if (unbiased)
                        {
                            for (size_t l=0u; l<Ly; ++l, Y+=K) { *Y /= (double)(Lx-l); }
                        }
                        else
                        {
                            for (size_t l=Ly; l>0u; --l, Y+=K) { *Y /= (double)Lx; }
                        }
                    }
                }
            }
        }
    }

    return 0;
}


int sig2ac_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Ly, const int mnz, const int unbiased, const int corr)
{
    if (dim>3u) { fprintf(stderr,"error in sig2ac_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Ly>Lx) { fprintf(stderr,"error in sig2ac_c: Ly (num lags in output) must be <= Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        if (Lx==N)
        {
            if (mnz)
            {
                float mnr = 0.0f, mni = 0.0f;
                for (size_t l=Lx; l>0u; --l) { mnr += *X++; mni += *X++; }
                mnr /= (float)Lx; mni /= (float)Lx;
                for (size_t l=Lx; l>0u; --l) { *--X -= mni; *--X -= mnr; }
            }

            if (Ly>850u && Lx>80000u) //I also leave this solution since it is closer to real-time use
            {
                float xr, xi;
                for (size_t l=2u*Ly; l>0u; --l, ++Y) { *Y = 0.0f; }
                Y -= 2u*Ly;
                for (size_t n=0u; n<N; ++n)
                {
                    xr = *X++; xi = *X++;
                    for (size_t l=0u; l<Ly && l<N-n; ++l)
                    {
                        Y[2u*l] += xr*X[2u*l] + xi*X[2u*l+1u];
                        Y[2u*l+1u] += xi*X[2u*l] - xr*X[2u*l+1u];
                    }
                }
            }
            else
            {
                float smr, smi;
                for (size_t l=0u; l<Ly; ++l, X-=2u*(Lx-l+1u))
                {
                    smr = smi = 0.0f;
                    for (size_t n=Lx-l; n>0u; --n, X+=2)
                    {
                        smr += *X**(X+2u*l) + *(X+1)**(X+2u*l+1u);
                        smi += *(X+1)**(X+2u*l) - *X**(X+2u*l+1u);
                    }
                    *Y++ = smr; *Y++ = smi;
                }
                Y -= 2u*Ly;
            }
            
            if (corr)
            {
                const float y0r = *Y, y0i = *(Y+1), y0a = y0r*y0r + y0i*y0i;
                float yr, yi;
                *Y++ = 1.0f; *Y++ = 0.0f;
                for (size_t l=Ly; l>1u; ++l)
                {
                    yr = *Y; yi = *(Y+1);
                    *Y++ = (yr*y0r+yi*y0i) / y0a;
                    *Y++ = (yi*y0r-yr*y0i) / y0a;
                }
            }
            else if (unbiased)
            {
                for (size_t l=0u; l<Ly; ++l) { *Y++ /= (float)(Lx-l); *Y++ /= (float)(Lx-l); }
            }
            else
            {
                for (size_t l=2u*Ly; l>0u; --l, ++Y) { *Y /= (float)Lx; }
            }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;
            float smr, smi;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=2u*Ly-2u)
                {
                    if (mnz)
                    {
                        float mnr = 0.0f, mni = 0.0f;
                        for (size_t l=Lx; l>0u; --l) { mnr += *X++; mni += *X++; }
                        mnr /= (float)Lx; mni /= (float)Lx;
                        for (size_t l=Lx; l>0u; --l) { *--X -= mni; *--X -= mnr; }
                    }

                    for (size_t l=0u; l<Ly-1u; ++l, X-=2u*(Lx-l+1u))
                    {
                        smr = smi = 0.0f;
                        for (size_t n=Lx-l; n>0u; --n, X+=2)
                        {
                            smr += *X**(X+2u*l) + *(X+1u)**(X+2u*l+1u);
                            smi += *(X+1u)**(X+2u*l) - *X**(X+2u*l+1u);
                        }
                        *Y++ = smr; *Y++ = smi;
                    }
                    smr = smi = 0.0f;
                    for (size_t n=Lx-Ly+1u; n>0u; --n, X+=2)
                    {
                        smr += *X**(X+2u*Ly-2u) + *(X+1u)**(X+2u*Ly-1u);
                        smi += *(X+1u)**(X+2u*Ly-2u) - *X**(X+2u*Ly-1u);
                    }
                    *Y = smr; *(Y+1) = smi; Y -= 2u*Ly-2u;

                    if (corr)
                    {
                        const float y0r = *Y, y0i = *(Y+1), y0a = y0r*y0r + y0i*y0i;
                        float yr, yi;
                        *Y++ = 1.0f; *Y++ = 0.0f;
                        for (size_t l=Ly; l>1u; ++l)
                        {
                            yr = *Y; yi = *(Y+1);
                            *Y++ = (yr*y0r+yi*y0i) / y0a;
                            *Y++ = (yi*y0r-yr*y0i) / y0a;
                        }
                    }
                    else if (unbiased)
                    {
                        for (size_t l=0u; l<Ly; ++l) { *Y++ /= (float)(Lx-l); *Y++ /= (float)(Lx-l); }
                    }
                    else
                    {
                        for (size_t l=2u*Ly; l>0u; --l, ++Y) { *Y /= (float)Lx; }
                    }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2, Y-=2u*K*Ly-2u)
                    {
                        if (mnz)
                        {
                            float mnr = 0.0f, mni = 0.0f;
                            for (size_t l=Lx; l>0u; --l, X+=2u*K) { mnr += *X; mni += *(X+1); }
                            mnr /= (float)Lx; mni /= (float)Lx;
                            for (size_t l=Lx; l>0u; --l) { X-=2u*K; *X -= mnr; *(X+1) -= mni; }
                        }

                        for (size_t l=0u; l<Ly; ++l, X-=2u*K*(Lx-l+1u), Y+=2u*K)
                        {
                            smr = smi = 0.0f;
                            for (size_t n=Lx-l; n>0u; --n, X+=2u*K)
                            {
                                smr += *X**(X+2u*l*K) + *(X+1)**(X+2u*l*K+1u);
                                smi += *(X+1)**(X+2u*l*K) - *X**(X+2u*l*K+1u);
                            }
                            *Y = smr; *(Y+1) = smi;
                        }
                        Y -= 2u*K*Ly;

                        if (corr)
                        {
                            const float y0r = *Y, y0i = *(Y+1), y0a = y0r*y0r + y0i*y0i;
                            float yr, yi;
                            *Y = 1.0f; *(Y+1) = 0.0f; Y += 2u*K;
                            for (size_t l=Ly; l>1u; --l, Y+=2u*K)
                            {
                                yr = *Y; yi = *(Y+1);
                                *Y = (yr*y0r+yi*y0i) / y0a;
                                *(Y+1) = (yi*y0r-yr*y0i) / y0a;
                            }
                        }
                        else if (unbiased)
                        {
                            for (size_t l=0u; l<Ly; ++l, Y+=2u*K) { *Y /= (float)(Lx-l); *(Y+1) /= (float)(Lx-l); }
                        }
                        else
                        {
                            for (size_t l=Ly; l>0u; --l, Y+=2u*K) { *Y /= (float)Lx; *(Y+1) /= (float)Lx; }
                        }
                    }
                }
            }
        }
    }

    return 0;
}


int sig2ac_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Ly, const int mnz, const int unbiased, const int corr)
{
    if (dim>3u) { fprintf(stderr,"error in sig2ac_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Ly>Lx) { fprintf(stderr,"error in sig2ac_z: Ly (num lags in output) must be <= Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        if (Lx==N)
        {
            if (mnz)
            {
                double mnr = 0.0, mni = 0.0;
                for (size_t l=Lx; l>0u; --l) { mnr += *X++; mni += *X++; }
                mnr /= (double)Lx; mni /= (double)Lx;
                for (size_t l=Lx; l>0u; --l) { *--X -= mni; *--X -= mnr; }
            }

            if (Ly>850u && Lx>80000u) //I also leave this solution since it is closer to real-time use
            {
                double xr, xi;
                for (size_t l=2u*Ly; l>0u; --l, ++Y) { *Y = 0.0; }
                Y -= 2u*Ly;
                for (size_t n=0u; n<N; ++n)
                {
                    xr = *X++; xi = *X++;
                    for (size_t l=0u; l<Ly && l<N-n; ++l)
                    {
                        Y[2u*l] += xr*X[2u*l] + xi*X[2u*l+1u];
                        Y[2u*l+1u] += xi*X[2u*l] - xr*X[2u*l+1u];
                    }
                }
            }
            else
            {
                double smr, smi;
                for (size_t l=0u; l<Ly; ++l, X-=2u*(Lx-l+1u))
                {
                    smr = smi = 0.0;
                    for (size_t n=Lx-l; n>0u; --n, X+=2)
                    {
                        smr += *X**(X+2u*l) + *(X+1)**(X+2u*l+1u);
                        smi += *(X+1)**(X+2u*l) - *X**(X+2u*l+1u);
                    }
                    *Y++ = smr; *Y++ = smi;
                }
                Y -= 2u*Ly;
            }
            
            if (corr)
            {
                const double y0r = *Y, y0i = *(Y+1), y0a = y0r*y0r + y0i*y0i;
                double yr, yi;
                *Y++ = 1.0; *Y++ = 0.0;
                for (size_t l=Ly; l>1u; ++l)
                {
                    yr = *Y; yi = *(Y+1);
                    *Y++ = (yr*y0r+yi*y0i) / y0a;
                    *Y++ = (yi*y0r-yr*y0i) / y0a;
                }
            }
            else if (unbiased)
            {
                for (size_t l=0u; l<Ly; ++l) { *Y++ /= (double)(Lx-l); *Y++ /= (double)(Lx-l); }
            }
            else
            {
                for (size_t l=2u*Ly; l>0u; --l, ++Y) { *Y /= (double)Lx; }
            }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;
            double smr, smi;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=2u*Ly-2u)
                {
                    if (mnz)
                    {
                        double mnr = 0.0, mni = 0.0;
                        for (size_t l=Lx; l>0u; --l) { mnr += *X++; mni += *X++; }
                        mnr /= (double)Lx; mni /= (double)Lx;
                        for (size_t l=Lx; l>0u; --l) { *--X -= mni; *--X -= mnr; }
                    }

                    for (size_t l=0u; l<Ly-1u; ++l, X-=2u*(Lx-l+1u))
                    {
                        smr = smi = 0.0;
                        for (size_t n=Lx-l; n>0u; --n, X+=2)
                        {
                            smr += *X**(X+2u*l) + *(X+1u)**(X+2u*l+1u);
                            smi += *(X+1u)**(X+2u*l) - *X**(X+2u*l+1u);
                        }
                        *Y++ = smr; *Y++ = smi;
                    }
                    smr = smi = 0.0;
                    for (size_t n=Lx-Ly+1u; n>0u; --n, X+=2)
                    {
                        smr += *X**(X+2u*Ly-2u) + *(X+1u)**(X+2u*Ly-1u);
                        smi += *(X+1u)**(X+2u*Ly-2u) - *X**(X+2u*Ly-1u);
                    }
                    *Y = smr; *(Y+1) = smi; Y -= 2u*Ly-2u;

                    if (corr)
                    {
                        const double y0r = *Y, y0i = *(Y+1), y0a = y0r*y0r + y0i*y0i;
                        double yr, yi;
                        *Y++ = 1.0; *Y++ = 0.0;
                        for (size_t l=Ly; l>1u; ++l)
                        {
                            yr = *Y; yi = *(Y+1);
                            *Y++ = (yr*y0r+yi*y0i) / y0a;
                            *Y++ = (yi*y0r-yr*y0i) / y0a;
                        }
                    }
                    else if (unbiased)
                    {
                        for (size_t l=0u; l<Ly; ++l) { *Y++ /= (double)(Lx-l); *Y++ /= (double)(Lx-l); }
                    }
                    else
                    {
                        for (size_t l=2u*Ly; l>0u; --l, ++Y) { *Y /= (double)Lx; }
                    }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2, Y-=2u*K*Ly-2u)
                    {
                        if (mnz)
                        {
                            double mnr = 0.0, mni = 0.0;
                            for (size_t l=Lx; l>0u; --l, X+=2u*K) { mnr += *X; mni += *(X+1); }
                            mnr /= (double)Lx; mni /= (double)Lx;
                            for (size_t l=Lx; l>0u; --l) { X-=2u*K; *X -= mnr; *(X+1) -= mni; }
                        }

                        for (size_t l=0u; l<Ly; ++l, X-=2u*K*(Lx-l+1u), Y+=2u*K)
                        {
                            smr = smi = 0.0;
                            for (size_t n=Lx-l; n>0u; --n, X+=2u*K)
                            {
                                smr += *X**(X+2u*l*K) + *(X+1)**(X+2u*l*K+1u);
                                smi += *(X+1)**(X+2u*l*K) - *X**(X+2u*l*K+1u);
                            }
                            *Y = smr; *(Y+1) = smi;
                        }
                        Y -= 2u*K*Ly;

                        if (corr)
                        {
                            const double y0r = *Y, y0i = *(Y+1), y0a = y0r*y0r + y0i*y0i;
                            double yr, yi;
                            *Y = 1.0; *(Y+1) = 0.0; Y += 2u*K;
                            for (size_t l=Ly; l>1u; --l, Y+=2u*K)
                            {
                                yr = *Y; yi = *(Y+1);
                                *Y = (yr*y0r+yi*y0i) / y0a;
                                *(Y+1) = (yi*y0r-yr*y0i) / y0a;
                            }
                        }
                        else if (unbiased)
                        {
                            for (size_t l=0u; l<Ly; ++l, Y+=2u*K) { *Y /= (double)(Lx-l); *(Y+1) /= (double)(Lx-l); }
                        }
                        else
                        {
                            for (size_t l=Ly; l>0u; --l, Y+=2u*K) { *Y /= (double)Lx; *(Y+1) /= (double)Lx; }
                        }
                    }
                }
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
