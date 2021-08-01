//Gets autocovariance (AC) for lags 0 to L-1 for each vector in X.
//The means are NOT subtracted, and lag 0 of Y is NOT normalized,
//which matches the behavior of Octave's xcorr.

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

int sig2ac_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t L, const int unbiased, const int corr);
int sig2ac_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t L, const int unbiased, const int corr);
int sig2ac_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t L, const int unbiased, const int corr);
int sig2ac_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t L, const int unbiased, const int corr);


int sig2ac_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t L, const int unbiased, const int corr)
{
    if (dim>3u) { fprintf(stderr,"error in sig2ac_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L>Lx) { fprintf(stderr,"error in sig2ac_s: L (num lags in output) must be <= Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        if (Lx==N)
        {
            if (L>850u && Lx>80000u) //I also leave this solution since it is closer to real-time use
            {
                float x;
                for (size_t l=0u; l<L; ++l, ++Y) { *Y = 0.0f; }
                Y -= L;
                for (size_t n=0u; n<N; ++n, ++X)
                {
                    x = *X;
                    for (size_t l=0u; l<L && l<N-n; ++l) { Y[l] += x * X[l]; }
                }
            }
            else
            {
                float sm;
                for (size_t l=0u; l<L; ++l, X-=Lx-l+1u, ++Y)
                {
                    sm = 0.0f;
                    for (size_t n=0u; n<Lx-l; ++n, ++X) { sm += *X * *(X+l); }
                    *Y = sm;
                }
                Y -= L;
            }

            if (corr)
            {
                const float y0 = *Y;
                *Y++ = 1.0f;
                for (size_t l=1u; l<L; ++l, ++Y) { *Y /= y0; }
            }
            else if (unbiased)
            {
                for (size_t l=0u; l<L; ++l, ++Y) { *Y /= (float)(Lx-l); }
            }
            else //xcorr leaves this blank
            {
                for (size_t l=0u; l<L; ++l, ++Y) { *Y /= (float)Lx; }
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
                        for (size_t l=0u; l<L; ++l, ++Y) { *Y /= (float)(Lx-l); }
                    }
                    else
                    {
                        for (size_t l=0u; l<L; ++l, ++Y) { *Y /= (float)Lx; }
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
                            for (size_t l=0u; l<L; ++l, Y+=K) { *Y /= (float)(Lx-l); }
                        }
                        else
                        {
                            for (size_t l=0u; l<L; ++l, Y+=K) { *Y /= (float)Lx; }
                        }
                    }
                }
            }
        }
    }

    return 0;
}


int sig2ac_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t L, const int unbiased, const int corr)
{
    if (dim>3u) { fprintf(stderr,"error in sig2ac_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L>Lx) { fprintf(stderr,"error in sig2ac_d: L (num lags in output) must be <= Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        if (Lx==N)
        {
            if (L>850u && Lx>80000u) //I also leave this solution since it is closer to real-time use
            {
                double x;
                for (size_t l=0u; l<L; ++l, ++Y) { *Y = 0.0; }
                Y -= L;
                for (size_t n=0u; n<N; ++n, ++X)
                {
                    x = *X;
                    for (size_t l=0u; l<L && l<N-n; ++l) { Y[l] += x * X[l]; }
                }
            }
            else
            {
                double sm;
                for (size_t l=0u; l<L; ++l, X-=Lx-l+1u, ++Y)
                {
                    sm = 0.0;
                    for (size_t n=0u; n<Lx-l; ++n, ++X) { sm += *X * *(X+l); }
                    *Y = sm;
                }
                Y -= L;
            }

            if (corr)
            {
                const double y0 = *Y;
                *Y++ = 1.0;
                for (size_t l=1u; l<L; ++l, ++Y) { *Y /= y0; }
            }
            else if (unbiased)
            {
                for (size_t l=0u; l<L; ++l, ++Y) { *Y /= (double)(Lx-l); }
            }
            else //xcorr leaves this blank
            {
                for (size_t l=0u; l<L; ++l, ++Y) { *Y /= (double)Lx; }
            }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;
            double sm, y0;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, X+=L-1u)
                {
                    for (size_t l=0u; l<L-1u; ++l, X-=Lx-l+1u, ++Y)
                    {
                        sm = 0.0;
                        for (size_t n=0u; n<Lx-l; ++n, ++X) { sm += *X * *(X+l); }
                        *Y = sm;
                    }
                    sm = 0.0;
                    for (size_t n=0u; n<Lx-L+1u; ++n, ++X) { sm += *X * *(X+L-1u); }
                    *Y = sm; Y -= L-1u;

                    if (corr)
                    {
                        y0 = *Y; *Y++ = 1.0;
                        for (size_t l=1u; l<L; ++l, ++Y) { *Y /= y0; }
                    }
                    else if (unbiased)
                    {
                        for (size_t l=0u; l<L; ++l, ++Y) { *Y /= (double)(Lx-l); }
                    }
                    else
                    {
                        for (size_t l=0u; l<L; ++l, ++Y) { *Y /= (double)Lx; }
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
                            sm = 0.0;
                            for (size_t n=0u; n<Lx-l; ++n, X+=K) { sm += *X * *(X+l*K); }
                            *Y = sm;
                        }
                        Y -= K*L;

                        if (corr)
                        {
                            y0 = *Y; *Y = 1.0; Y += K;
                            for (size_t l=1u; l<L; ++l, Y+=K) { *Y /= y0; }
                        }
                        else if (unbiased)
                        {
                            for (size_t l=0u; l<L; ++l, Y+=K) { *Y /= (double)(Lx-l); }
                        }
                        else
                        {
                            for (size_t l=0u; l<L; ++l, Y+=K) { *Y /= (double)Lx; }
                        }
                    }
                }
            }
        }
    }

    return 0;
}


int sig2ac_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t L, const int unbiased, const int corr)
{
    if (dim>3u) { fprintf(stderr,"error in sig2ac_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L>Lx) { fprintf(stderr,"error in sig2ac_c: L (num lags in output) must be <= Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        if (Lx==N)
        {
            if (L>850u && Lx>80000u) //I also leave this solution since it is closer to real-time use
            {
                float xr, xi;
                for (size_t l=0u; l<2u*L; ++l, ++Y) { *Y = 0.0f; }
                Y -= 2u*L;
                for (size_t n=0u; n<N; ++n)
                {
                    xr = *X++; xi = *X++;
                    for (size_t l=0u; l<L && l<N-n; ++l)
                    {
                        Y[2u*l] += xr*X[2u*l] + xi*X[2u*l+1u];
                        Y[2u*l+1u] += xi*X[2u*l] - xr*X[2u*l+1u];
                    }
                }
            }
            else
            {
                float smr, smi;
                for (size_t l=0u; l<L; ++l, X-=2u*(Lx-l+1u))
                {
                    smr = smi = 0.0f;
                    for (size_t n=0u; n<Lx-l; ++n, X+=2)
                    {
                        smr += *X**(X+2u*l) + *(X+1)**(X+2u*l+1u);
                        smi += *(X+1)**(X+2u*l) - *X**(X+2u*l+1u);
                    }
                    *Y++ = smr; *Y++ = smi;
                }
                Y -= 2u*L;
            }
            
            if (corr)
            {
                const float y0r = *Y, y0i = *(Y+1), y0a = y0r*y0r + y0i*y0i;
                float yr, yi;
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
                for (size_t l=0u; l<L; ++l) { *Y++ /= (float)(Lx-l); *Y++ /= (float)(Lx-l); }
            }
            else
            {
                for (size_t l=0u; l<2u*L; ++l, ++Y) { *Y /= (float)Lx; }
            }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;
            float smr, smi, y0r, y0i, y0a, yr, yi;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, X+=2u*L-2u)
                {
                    for (size_t l=0u; l<L-1u; ++l, X-=2u*(Lx-l+1u))
                    {
                        smr = smi = 0.0f;
                        for (size_t n=0u; n<Lx-l; ++n, X+=2)
                        {
                            smr += *X**(X+2u*l) + *(X+1u)**(X+2u*l+1u);
                            smi += *(X+1u)**(X+2u*l) - *X**(X+2u*l+1u);
                        }
                        *Y++ = smr; *Y++ = smi;
                    }
                    smr = smi = 0.0f;
                    for (size_t n=0u; n<Lx-L+1u; ++n, X+=2)
                    {
                        smr += *X**(X+2u*L-2u) + *(X+1u)**(X+2u*L-1u);
                        smi += *(X+1u)**(X+2u*L-2u) - *X**(X+2u*L-1u);
                    }
                    *Y = smr; *(Y+1) = smi; Y -= 2u*L-2u;

                    if (corr)
                    {
                        y0r = *Y; y0i = *(Y+1);
                        y0a = y0r*y0r + y0i*y0i;
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
                        for (size_t l=0u; l<L; ++l) { *Y++ /= (float)(Lx-l); *Y++ /= (float)(Lx-l); }
                    }
                    else
                    {
                        for (size_t l=0u; l<2u*L; ++l, ++Y) { *Y /= (float)Lx; }
                    }
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X+=2, Y-=2u*K*L-2u)
                    {
                        for (size_t l=0u; l<L; ++l, X-=2u*K*(Lx-l+1u), Y+=2u*K)
                        {
                            smr = smi = 0.0f;
                            for (size_t n=0u; n<Lx-l; ++n, X+=2u*K)
                            {
                                smr += *X**(X+2u*l*K) + *(X+1)**(X+2u*l*K+1u);
                                smi += *(X+1)**(X+2u*l*K) - *X**(X+2u*l*K+1u);
                            }
                            *Y = smr; *(Y+1) = smi;
                        }
                        Y -= 2u*K*L;

                        if (corr)
                        {
                            y0r = *Y; y0i = *(Y+1);
                            y0a = y0r*y0r + y0i*y0i;
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
                            for (size_t l=0u; l<L; ++l, Y+=2u*K) { *Y /= (float)(Lx-l); *(Y+1) /= (float)(Lx-l); }
                        }
                        else
                        {
                            for (size_t l=0u; l<L; ++l, Y+=2u*K) { *Y /= (float)Lx; *(Y+1) /= (float)Lx; }
                        }
                    }
                }
            }
        }
    }

    return 0;
}


int sig2ac_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t L, const int unbiased, const int corr)
{
    if (dim>3u) { fprintf(stderr,"error in sig2ac_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L>Lx) { fprintf(stderr,"error in sig2ac_z: L (num lags in output) must be <= Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        if (Lx==N)
        {
            if (L>850u && Lx>80000u) //I also leave this solution since it is closer to real-time use
            {
                double xr, xi;
                for (size_t l=0u; l<2u*L; ++l, ++Y) { *Y = 0.0; }
                Y -= 2u*L;
                for (size_t n=0u; n<N; ++n)
                {
                    xr = *X++; xi = *X++;
                    for (size_t l=0u; l<L && l<N-n; ++l)
                    {
                        Y[2u*l] += xr*X[2u*l] + xi*X[2u*l+1u];
                        Y[2u*l+1u] += xi*X[2u*l] - xr*X[2u*l+1u];
                    }
                }
            }
            else
            {
                double smr, smi;
                for (size_t l=0u; l<L; ++l, X-=2u*(Lx-l+1u))
                {
                    smr = smi = 0.0;
                    for (size_t n=0u; n<Lx-l; ++n, X+=2)
                    {
                        smr += *X**(X+2u*l) + *(X+1)**(X+2u*l+1u);
                        smi += *(X+1)**(X+2u*l) - *X**(X+2u*l+1u);
                    }
                    *Y++ = smr; *Y++ = smi;
                }
                Y -= 2u*L;
            }
            
            if (corr)
            {
                const double y0r = *Y, y0i = *(Y+1), y0a = y0r*y0r + y0i*y0i;
                double yr, yi;
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
                for (size_t l=0u; l<L; ++l) { *Y++ /= (double)(Lx-l); *Y++ /= (double)(Lx-l); }
            }
            else
            {
                for (size_t l=0u; l<2u*L; ++l, ++Y) { *Y /= (double)Lx; }
            }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;
            double smr, smi, y0r, y0i, y0a, yr, yi;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, X+=2u*L-2u)
                {
                    for (size_t l=0u; l<L-1u; ++l, X-=2u*(Lx-l+1u))
                    {
                        smr = smi = 0.0;
                        for (size_t n=0u; n<Lx-l; ++n, X+=2)
                        {
                            smr += *X**(X+2u*l) + *(X+1u)**(X+2u*l+1u);
                            smi += *(X+1u)**(X+2u*l) - *X**(X+2u*l+1u);
                        }
                        *Y++ = smr; *Y++ = smi;
                    }
                    smr = smi = 0.0;
                    for (size_t n=0u; n<Lx-L+1u; ++n, X+=2)
                    {
                        smr += *X**(X+2u*L-2u) + *(X+1u)**(X+2u*L-1u);
                        smi += *(X+1u)**(X+2u*L-2u) - *X**(X+2u*L-1u);
                    }
                    *Y = smr; *(Y+1) = smi; Y -= 2u*L-2u;

                    if (corr)
                    {
                        y0r = *Y; y0i = *(Y+1);
                        y0a = y0r*y0r + y0i*y0i;
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
                        for (size_t l=0u; l<L; ++l) { *Y++ /= (double)(Lx-l); *Y++ /= (double)(Lx-l); }
                    }
                    else
                    {
                        for (size_t l=0u; l<2u*L; ++l, ++Y) { *Y /= (double)Lx; }
                    }
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X+=2, Y-=2u*K*L-2u)
                    {
                        for (size_t l=0u; l<L; ++l, X-=2u*K*(Lx-l+1u), Y+=2u*K)
                        {
                            smr = smi = 0.0;
                            for (size_t n=0u; n<Lx-l; ++n, X+=2u*K)
                            {
                                smr += *X**(X+2u*l*K) + *(X+1)**(X+2u*l*K+1u);
                                smi += *(X+1)**(X+2u*l*K) - *X**(X+2u*l*K+1u);
                            }
                            *Y = smr; *(Y+1) = smi;
                        }
                        Y -= 2u*K*L;

                        if (corr)
                        {
                            y0r = *Y; y0i = *(Y+1);
                            y0a = y0r*y0r + y0i*y0i;
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
                            for (size_t l=0u; l<L; ++l, Y+=2u*K) { *Y /= (double)(Lx-l); *(Y+1) /= (double)(Lx-l); }
                        }
                        else
                        {
                            for (size_t l=0u; l<L; ++l, Y+=2u*K) { *Y /= (double)Lx; *(Y+1) /= (double)Lx; }
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
