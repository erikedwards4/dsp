//Gets the zero-crossing rate (ZCR) for each vec in X along dim.
//See zcs.c for comments on getting the zero-crossings (ZCs).

//The rate is the count of ZCs within a window of length L, divided by L.
//For causal, the current and previous L-1 samps are used.
//For noncausal, the current samp is the center of the window.

//The output Y has the same size (and therefore sample rate) as the input X.
//Y is in units of ZCs/sample, so multiply by sample rate to get ZCs/sec.
//Y is thus in [0.0 1.0], since at most 1 ZC can occur per sample.

//Edge windows can be handled by at least three methods:
//1. Unbiased: divide by num samps observed in the edge windows (but unstable at very edge).
//2. Biased: divided by Lw for all windows (gives nice transition from 0 to full-window rate).
//3. Extrapolate: just continue the nearest full-window result (fastest and easiest).

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int zcr_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const int causal);
int zcr_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const int causal);
int zcr_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const int causal);
int zcr_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const int causal);


int zcr_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const int causal)
{
    if (dim>3u) { fprintf(stderr,"error in zcr_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lw>L) { fprintf(stderr,"error in zcr_s: Lw (winlength) must be <= L (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        const float Lwf = (float)Lw;
        float y;
        int *Z, s, sp, sm;
        if (!(Z=(int *)malloc(L*sizeof(int)))) { fprintf(stderr,"error in zcr_s: problem with malloc. "); perror("malloc"); return 1; }

        if (L==N)
        {
            if (going==0)
            {
                sp = (*X++<0.0f); *Z++ = sm = 0;
                for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X<0.0f); sm += (s!=sp); *Z = sm; sp = s; }
                for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                y = (float)(*Z) / Lwf;
                for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (float)(*Z) / Lwf; }
                Z -= L;
                if (!causal)
                {
                    const size_t L2 = Lw - Lw/2u - 1u;
                    Y -= L;
                    for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                    y = *(Y-1);
                    for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                }
            }
            else if (going==1)
            {
                sp = (*X++>=0.0f); *Z++ = sm = 0;
                for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X>=0.0f); sm += s*(s!=sp); *Z = sm; sp = s; }
                for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                y = (float)(*Z) / Lwf;
                for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (float)(*Z) / Lwf; }
                Z -= L;
                if (!causal)
                {
                    const size_t L2 = Lw - Lw/2u - 1u;
                    Y -= L;
                    for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                    y = *(Y-1);
                    for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                }
            }
            else if (going==-1)
            {
                sp = (*X++<0.0f); *Z++ = sm = 0;
                for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X<0.0f); sm += s*(s!=sp); *Z = sm; sp = s; }
                for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                y = (float)(*Z) / Lwf;
                for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (float)(*Z) / Lwf; }
                Z -= L;
                if (!causal)
                {
                    const size_t L2 = Lw - Lw/2u - 1u;
                    Y -= L;
                    for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                    y = *(Y-1);
                    for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                }
            }
            else { fprintf(stderr,"error in zcr_s: going must be in {-1,0,1}\n"); return 1; }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                if (going==0)
                {
                    for (size_t v=0u; v<V; ++v)
                    {
                        sp = (*X++<0.0f); *Z++ = sm = 0;
                        for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X<0.0f); sm += (s!=sp); *Z = sm; sp = s; }
                        for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                        y = (float)(*Z) / Lwf;
                        for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                        for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (float)(*Z) / Lwf; }
                        Z -= L;
                        if (!causal)
                        {
                            const size_t L2 = Lw - Lw/2u - 1u;
                            Y -= L;
                            for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                            y = *(Y-1);
                            for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                        }
                    }
                }
                else if (going==1)
                {
                    for (size_t v=0u; v<V; ++v)
                    {
                        sp = (*X++>=0.0f); *Z++ = sm = 0;
                        for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X>=0.0f); sm += s*(s!=sp); *Z = sm; sp = s; }
                        for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                        y = (float)(*Z) / Lwf;
                        for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                        for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (float)(*Z) / Lwf; }
                        Z -= L;
                        if (!causal)
                        {
                            const size_t L2 = Lw - Lw/2u - 1u;
                            Y -= L;
                            for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                            y = *(Y-1);
                            for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                        }
                    }
                }
                else if (going==-1)
                {
                    for (size_t v=0u; v<V; ++v)
                    {
                        sp = (*X++<0.0f); *Z++ = sm = 0;
                        for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X<0.0f); sm += s*(s!=sp); *Z = sm; sp = s; }
                        for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                        y = (float)(*Z) / Lwf;
                        for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                        for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (float)(*Z) / Lwf; }
                        Z -= L;
                        if (!causal)
                        {
                            const size_t L2 = Lw - Lw/2u - 1u;
                            Y -= L;
                            for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                            y = *(Y-1);
                            for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                        }
                    }
                }
                else { fprintf(stderr,"error in zcr_s: going must be in {-1,0,1}\n"); return 1; }
            }
            else
            {
                if (going==0)
                {
                    for (size_t g=0u; g<G; ++g, X+=B*(L-1u), Y+=B*(L-1u))
                    {
                        for (size_t b=0u; b<B; ++b, X-=K*L-1u, Y-=K*L-1u)
                        {
                            sp = (*X<0.0f); X+=K; *Z++ = sm = 0;
                            for (size_t l=1u; l<L; ++l, X+=K, ++Z) { s = (*X<0.0f); sm += (s!=sp); *Z = sm; sp = s; }
                            for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                            y = (float)(*Z) / Lwf;
                            for (size_t l=0u; l<Lw-1u; ++l, Y+=K) { *Y = y; }
                            for (size_t l=Lw-1u; l<L; ++l, Y+=K, ++Z) { *Y = (float)(*Z) / Lwf; }
                            Z -= L;
                            if (!causal)
                            {
                                const size_t L2 = Lw - Lw/2u - 1u;
                                Y -= K*L;
                                for (size_t l=0u; l<L-L2; ++l, Y+=K) { *Y = *(Y+K*L2); }
                                y = *(Y-K);
                                for (size_t l=0u; l<L2; ++l, Y+=K) { *Y = y; }
                            }
                        }
                    }
                }
                else if (going==1)
                {
                    for (size_t g=0u; g<G; ++g, X+=B*(L-1u), Y+=B*(L-1u))
                    {
                        for (size_t b=0u; b<B; ++b, X-=K*L-1u, Y-=K*L-1u)
                        {
                            sp = (*X>=0.0f); X+=K; *Z++ = sm = 0;
                            for (size_t l=1u; l<L; ++l, X+=K, ++Z) { s = (*X>=0.0f); sm += s*(s!=sp); *Z = sm; sp = s; }
                            for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                            y = (float)(*Z) / Lwf;
                            for (size_t l=0u; l<Lw-1u; ++l, Y+=K) { *Y = y; }
                            for (size_t l=Lw-1u; l<L; ++l, Y+=K, ++Z) { *Y = (float)(*Z) / Lwf; }
                            Z -= L;
                            if (!causal)
                            {
                                const size_t L2 = Lw - Lw/2u - 1u;
                                Y -= K*L;
                                for (size_t l=0u; l<L-L2; ++l, Y+=K) { *Y = *(Y+K*L2); }
                                y = *(Y-K);
                                for (size_t l=0u; l<L2; ++l, Y+=K) { *Y = y; }
                            }
                        }
                    }
                }
                else if (going==-1)
                {
                    for (size_t g=0u; g<G; ++g, X+=B*(L-1u), Y+=B*(L-1u))
                    {
                        for (size_t b=0u; b<B; ++b, X-=K*L-1u, Y-=K*L-1u)
                        {
                            sp = (*X<0.0f); X+=K; *Z++ = sm = 0;
                            for (size_t l=1u; l<L; ++l, X+=K, ++Z) { s = (*X<0.0f); sm += s*(s!=sp); *Z = sm; sp = s; }
                            for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                            y = (float)(*Z) / Lwf;
                            for (size_t l=0u; l<Lw-1u; ++l, Y+=K) { *Y = y; }
                            for (size_t l=Lw-1u; l<L; ++l, Y+=K, ++Z) { *Y = (float)(*Z) / Lwf; }
                            Z -= L;
                            if (!causal)
                            {
                                const size_t L2 = Lw - Lw/2u - 1u;
                                Y -= K*L;
                                for (size_t l=0u; l<L-L2; ++l, Y+=K) { *Y = *(Y+K*L2); }
                                y = *(Y-K);
                                for (size_t l=0u; l<L2; ++l, Y+=K) { *Y = y; }
                            }
                        }
                    }
                }
                else { fprintf(stderr,"error in zcr_s: going must be in {-1,0,1}\n"); return 1; }
            }
        }
        free(Z);
    }

    return 0;
}


int zcr_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const int causal)
{
    if (dim>3u) { fprintf(stderr,"error in zcr_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lw>L) { fprintf(stderr,"error in zcr_d: Lw (winlength) must be <= L (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        const double Lwf = (double)Lw;
        double y;
        int *Z, s, sp, sm;
        if (!(Z=(int *)malloc(L*sizeof(int)))) { fprintf(stderr,"error in zcr_d: problem with malloc. "); perror("malloc"); return 1; }

        if (L==N)
        {
            if (going==0)
            {
                sp = (*X++<0.0); *Z++ = sm = 0;
                for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X<0.0); sm += (s!=sp); *Z = sm; sp = s; }
                for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                y = (double)(*Z) / Lwf;
                for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (double)(*Z) / Lwf; }
                Z -= L;
                if (!causal)
                {
                    const size_t L2 = Lw - Lw/2u - 1u;
                    Y -= L;
                    for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                    y = *(Y-1);
                    for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                }
            }
            else if (going==1)
            {
                sp = (*X++>=0.0); *Z++ = sm = 0;
                for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X>=0.0); sm += s*(s!=sp); *Z = sm; sp = s; }
                for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                y = (double)(*Z) / Lwf;
                for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (double)(*Z) / Lwf; }
                Z -= L;
                if (!causal)
                {
                    const size_t L2 = Lw - Lw/2u - 1u;
                    Y -= L;
                    for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                    y = *(Y-1);
                    for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                }
            }
            else if (going==-1)
            {
                sp = (*X++<0.0); *Z++ = sm = 0;
                for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X<0.0); sm += s*(s!=sp); *Z = sm; sp = s; }
                for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                y = (double)(*Z) / Lwf;
                for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (double)(*Z) / Lwf; }
                Z -= L;
                if (!causal)
                {
                    const size_t L2 = Lw - Lw/2u - 1u;
                    Y -= L;
                    for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                    y = *(Y-1);
                    for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                }
            }
            else { fprintf(stderr,"error in zcr_d: going must be in {-1,0,1}\n"); return 1; }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                if (going==0)
                {
                    for (size_t v=0u; v<V; ++v)
                    {
                        sp = (*X++<0.0); *Z++ = sm = 0;
                        for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X<0.0); sm += (s!=sp); *Z = sm; sp = s; }
                        for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                        y = (double)(*Z) / Lwf;
                        for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                        for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (double)(*Z) / Lwf; }
                        Z -= L;
                        if (!causal)
                        {
                            const size_t L2 = Lw - Lw/2u - 1u;
                            Y -= L;
                            for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                            y = *(Y-1);
                            for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                        }
                    }
                }
                else if (going==1)
                {
                    for (size_t v=0u; v<V; ++v)
                    {
                        sp = (*X++>=0.0); *Z++ = sm = 0;
                        for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X>=0.0); sm += s*(s!=sp); *Z = sm; sp = s; }
                        for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                        y = (double)(*Z) / Lwf;
                        for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                        for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (double)(*Z) / Lwf; }
                        Z -= L;
                        if (!causal)
                        {
                            const size_t L2 = Lw - Lw/2u - 1u;
                            Y -= L;
                            for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                            y = *(Y-1);
                            for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                        }
                    }
                }
                else if (going==-1)
                {
                    for (size_t v=0u; v<V; ++v)
                    {
                        sp = (*X++<0.0); *Z++ = sm = 0;
                        for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X<0.0); sm += s*(s!=sp); *Z = sm; sp = s; }
                        for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                        y = (double)(*Z) / Lwf;
                        for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                        for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (double)(*Z) / Lwf; }
                        Z -= L;
                        if (!causal)
                        {
                            const size_t L2 = Lw - Lw/2u - 1u;
                            Y -= L;
                            for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                            y = *(Y-1);
                            for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                        }
                    }
                }
                else { fprintf(stderr,"error in zcr_d: going must be in {-1,0,1}\n"); return 1; }
            }
            else
            {
                if (going==0)
                {
                    for (size_t g=0u; g<G; ++g, X+=B*(L-1u), Y+=B*(L-1u))
                    {
                        for (size_t b=0u; b<B; ++b, X-=K*L-1u, Y-=K*L-1u)
                        {
                            sp = (*X<0.0); X+=K; *Z++ = sm = 0;
                            for (size_t l=1u; l<L; ++l, X+=K, ++Z) { s = (*X<0.0); sm += (s!=sp); *Z = sm; sp = s; }
                            for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                            y = (double)(*Z) / Lwf;
                            for (size_t l=0u; l<Lw-1u; ++l, Y+=K) { *Y = y; }
                            for (size_t l=Lw-1u; l<L; ++l, Y+=K, ++Z) { *Y = (double)(*Z) / Lwf; }
                            Z -= L;
                            if (!causal)
                            {
                                const size_t L2 = Lw - Lw/2u - 1u;
                                Y -= K*L;
                                for (size_t l=0u; l<L-L2; ++l, Y+=K) { *Y = *(Y+K*L2); }
                                y = *(Y-K);
                                for (size_t l=0u; l<L2; ++l, Y+=K) { *Y = y; }
                            }
                        }
                    }
                }
                else if (going==1)
                {
                    for (size_t g=0u; g<G; ++g, X+=B*(L-1u), Y+=B*(L-1u))
                    {
                        for (size_t b=0u; b<B; ++b, X-=K*L-1u, Y-=K*L-1u)
                        {
                            sp = (*X>=0.0); X+=K; *Z++ = sm = 0;
                            for (size_t l=1u; l<L; ++l, X+=K, ++Z) { s = (*X>=0.0); sm += s*(s!=sp); *Z = sm; sp = s; }
                            for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                            y = (double)(*Z) / Lwf;
                            for (size_t l=0u; l<Lw-1u; ++l, Y+=K) { *Y = y; }
                            for (size_t l=Lw-1u; l<L; ++l, Y+=K, ++Z) { *Y = (double)(*Z) / Lwf; }
                            Z -= L;
                            if (!causal)
                            {
                                const size_t L2 = Lw - Lw/2u - 1u;
                                Y -= K*L;
                                for (size_t l=0u; l<L-L2; ++l, Y+=K) { *Y = *(Y+K*L2); }
                                y = *(Y-K);
                                for (size_t l=0u; l<L2; ++l, Y+=K) { *Y = y; }
                            }
                        }
                    }
                }
                else if (going==-1)
                {
                    for (size_t g=0u; g<G; ++g, X+=B*(L-1u), Y+=B*(L-1u))
                    {
                        for (size_t b=0u; b<B; ++b, X-=K*L-1u, Y-=K*L-1u)
                        {
                            sp = (*X<0.0); X+=K; *Z++ = sm = 0;
                            for (size_t l=1u; l<L; ++l, X+=K, ++Z) { s = (*X<0.0); sm += s*(s!=sp); *Z = sm; sp = s; }
                            for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                            y = (double)(*Z) / Lwf;
                            for (size_t l=0u; l<Lw-1u; ++l, Y+=K) { *Y = y; }
                            for (size_t l=Lw-1u; l<L; ++l, Y+=K, ++Z) { *Y = (double)(*Z) / Lwf; }
                            Z -= L;
                            if (!causal)
                            {
                                const size_t L2 = Lw - Lw/2u - 1u;
                                Y -= K*L;
                                for (size_t l=0u; l<L-L2; ++l, Y+=K) { *Y = *(Y+K*L2); }
                                y = *(Y-K);
                                for (size_t l=0u; l<L2; ++l, Y+=K) { *Y = y; }
                            }
                        }
                    }
                }
                else { fprintf(stderr,"error in zcr_d: going must be in {-1,0,1}\n"); return 1; }
            }
        }
        free(Z);
    }

    return 0;
}


int zcr_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const int causal)
{
    if (dim>3u) { fprintf(stderr,"error in zcr_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lw>L) { fprintf(stderr,"error in zcr_c: Lw (winlength) must be <= L (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        const float Lwf = (float)Lw;
        float y;
        int *Z, s, sp, sm;
        if (!(Z=(int *)malloc(L*sizeof(int)))) { fprintf(stderr,"error in zcr_c: problem with malloc. "); perror("malloc"); return 1; }
        
        //Comment this out to use real-valued ZCs
        ++X;

        if (L==N)
        {
            if (going==0)
            {
                sp = (*X<0.0f); X+=2; *Z++ = sm = 0;
                for (size_t l=1u; l<L; ++l, X+=2, ++Z) { s = (*X<0.0f); sm += (s!=sp); *Z = sm; sp = s; }
                for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                y = (float)(*Z) / Lwf;
                for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (float)(*Z) / Lwf; }
                Z -= L;
                if (!causal)
                {
                    const size_t L2 = Lw - Lw/2u - 1u;
                    Y -= L;
                    for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                    y = *(Y-1);
                    for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                }
            }
            else if (going==1)
            {
                sp = (*X>=0.0f); X+=2; *Z++ = sm = 0;
                for (size_t l=1u; l<L; ++l, X+=2, ++Z) { s = (*X>=0.0f); sm += s*(s!=sp); *Z = sm; sp = s; }
                for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                y = (float)(*Z) / Lwf;
                for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (float)(*Z) / Lwf; }
                Z -= L;
                if (!causal)
                {
                    const size_t L2 = Lw - Lw/2u - 1u;
                    Y -= L;
                    for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                    y = *(Y-1);
                    for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                }
            }
            else if (going==-1)
            {
                sp = (*X<0.0f); X+=2; *Z++ = sm = 0;
                for (size_t l=1u; l<L; ++l, X+=2, ++Z) { s = (*X<0.0f); sm += s*(s!=sp); *Z = sm; sp = s; }
                for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                y = (float)(*Z) / Lwf;
                for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (float)(*Z) / Lwf; }
                Z -= L;
                if (!causal)
                {
                    const size_t L2 = Lw - Lw/2u - 1u;
                    Y -= L;
                    for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                    y = *(Y-1);
                    for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                }
            }
            else { fprintf(stderr,"error in zcr_c: going must be in {-1,0,1}\n"); return 1; }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                if (going==0)
                {
                    for (size_t v=0u; v<V; ++v)
                    {
                        sp = (*X<0.0f); X+=2; *Z++ = sm = 0;
                        for (size_t l=1u; l<L; ++l, X+=2, ++Z) { s = (*X<0.0f); sm += (s!=sp); *Z = sm; sp = s; }
                        for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                        y = (float)(*Z) / Lwf;
                        for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                        for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (float)(*Z) / Lwf; }
                        Z -= L;
                        if (!causal)
                        {
                            const size_t L2 = Lw - Lw/2u - 1u;
                            Y -= L;
                            for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                            y = *(Y-1);
                            for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                        }
                    }
                }
                else if (going==1)
                {
                    for (size_t v=0u; v<V; ++v)
                    {
                        sp = (*X>=0.0f); X+=2; *Z++ = sm = 0;
                        for (size_t l=1u; l<L; ++l, X+=2, ++Z) { s = (*X>=0.0f); sm += s*(s!=sp); *Z = sm; sp = s; }
                        for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                        y = (float)(*Z) / Lwf;
                        for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                        for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (float)(*Z) / Lwf; }
                        Z -= L;
                        if (!causal)
                        {
                            const size_t L2 = Lw - Lw/2u - 1u;
                            Y -= L;
                            for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                            y = *(Y-1);
                            for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                        }
                    }
                }
                else if (going==-1)
                {
                    for (size_t v=0u; v<V; ++v)
                    {
                        sp = (*X<0.0f); X+=2; *Z++ = sm = 0;
                        for (size_t l=1u; l<L; ++l, X+=2, ++Z) { s = (*X<0.0f); sm += s*(s!=sp); *Z = sm; sp = s; }
                        for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                        y = (float)(*Z) / Lwf;
                        for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                        for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (float)(*Z) / Lwf; }
                        Z -= L;
                        if (!causal)
                        {
                            const size_t L2 = Lw - Lw/2u - 1u;
                            Y -= L;
                            for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                            y = *(Y-1);
                            for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                        }
                    }
                }
                else { fprintf(stderr,"error in zcr_c: going must be in {-1,0,1}\n"); return 1; }
            }
            else
            {
                if (going==0)
                {
                    for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u), Y+=B*(L-1u))
                    {
                        for (size_t b=0u; b<B; ++b, X-=2u*K*L-2u, Y-=K*L-1u)
                        {
                            sp = (*X<0.0f); X+=2u*K; *Z++ = sm = 0;
                            for (size_t l=1u; l<L; ++l, X+=2u*K, ++Z) { s = (*X<0.0f); sm += (s!=sp); *Z = sm; sp = s; }
                            for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                            y = (float)(*Z) / Lwf;
                            for (size_t l=0u; l<Lw-1u; ++l, Y+=K) { *Y = y; }
                            for (size_t l=Lw-1u; l<L; ++l, Y+=K, ++Z) { *Y = (float)(*Z) / Lwf; }
                            Z -= L;
                            if (!causal)
                            {
                                const size_t L2 = Lw - Lw/2u - 1u;
                                Y -= K*L;
                                for (size_t l=0u; l<L-L2; ++l, Y+=K) { *Y = *(Y+K*L2); }
                                y = *(Y-K);
                                for (size_t l=0u; l<L2; ++l, Y+=K) { *Y = y; }
                            }
                        }
                    }
                }
                else if (going==1)
                {
                    for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u), Y+=B*(L-1u))
                    {
                        for (size_t b=0u; b<B; ++b, X-=2u*K*L-2u, Y-=K*L-1u)
                        {
                            sp = (*X>=0.0f); X+=2u*K; *Z++ = sm = 0;
                            for (size_t l=1u; l<L; ++l, X+=2u*K, ++Z) { s = (*X>=0.0f); sm += s*(s!=sp); *Z = sm; sp = s; }
                            for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                            y = (float)(*Z) / Lwf;
                            for (size_t l=0u; l<Lw-1u; ++l, Y+=K) { *Y = y; }
                            for (size_t l=Lw-1u; l<L; ++l, Y+=K, ++Z) { *Y = (float)(*Z) / Lwf; }
                            Z -= L;
                            if (!causal)
                            {
                                const size_t L2 = Lw - Lw/2u - 1u;
                                Y -= K*L;
                                for (size_t l=0u; l<L-L2; ++l, Y+=K) { *Y = *(Y+K*L2); }
                                y = *(Y-K);
                                for (size_t l=0u; l<L2; ++l, Y+=K) { *Y = y; }
                            }
                        }
                    }
                }
                else if (going==-1)
                {
                    for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u), Y+=B*(L-1u))
                    {
                        for (size_t b=0u; b<B; ++b, X-=2u*K*L-2u, Y-=K*L-1u)
                        {
                            sp = (*X<0.0f); X+=2u*K; *Z++ = sm = 0;
                            for (size_t l=1u; l<L; ++l, X+=2u*K, ++Z) { s = (*X<0.0f); sm += s*(s!=sp); *Z = sm; sp = s; }
                            for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                            y = (float)(*Z) / Lwf;
                            for (size_t l=0u; l<Lw-1u; ++l, Y+=K) { *Y = y; }
                            for (size_t l=Lw-1u; l<L; ++l, Y+=K, ++Z) { *Y = (float)(*Z) / Lwf; }
                            Z -= L;
                            if (!causal)
                            {
                                const size_t L2 = Lw - Lw/2u - 1u;
                                Y -= K*L;
                                for (size_t l=0u; l<L-L2; ++l, Y+=K) { *Y = *(Y+K*L2); }
                                y = *(Y-K);
                                for (size_t l=0u; l<L2; ++l, Y+=K) { *Y = y; }
                            }
                        }
                    }
                }
                else { fprintf(stderr,"error in zcr_c: going must be in {-1,0,1}\n"); return 1; }
            }
        }
        free(Z);
    }

    return 0;
}


int zcr_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const int causal)
{
    if (dim>3u) { fprintf(stderr,"error in zcr_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lw>L) { fprintf(stderr,"error in zcr_z: Lw (winlength) must be <= L (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        const double Lwf = (double)Lw;
        double y;
        int *Z, s, sp, sm;
        if (!(Z=(int *)malloc(L*sizeof(int)))) { fprintf(stderr,"error in zcr_z: problem with malloc. "); perror("malloc"); return 1; }
        
        //Comment this out to use real-valued ZCs
        ++X;

        if (L==N)
        {
            if (going==0)
            {
                sp = (*X<0.0); X+=2; *Z++ = sm = 0;
                for (size_t l=1u; l<L; ++l, X+=2, ++Z) { s = (*X<0.0); sm += (s!=sp); *Z = sm; sp = s; }
                for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                y = (double)(*Z) / Lwf;
                for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (double)(*Z) / Lwf; }
                Z -= L;
                if (!causal)
                {
                    const size_t L2 = Lw - Lw/2u - 1u;
                    Y -= L;
                    for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                    y = *(Y-1);
                    for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                }
            }
            else if (going==1)
            {
                sp = (*X>=0.0); X+=2; *Z++ = sm = 0;
                for (size_t l=1u; l<L; ++l, X+=2, ++Z) { s = (*X>=0.0); sm += s*(s!=sp); *Z = sm; sp = s; }
                for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                y = (double)(*Z) / Lwf;
                for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (double)(*Z) / Lwf; }
                Z -= L;
                if (!causal)
                {
                    const size_t L2 = Lw - Lw/2u - 1u;
                    Y -= L;
                    for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                    y = *(Y-1);
                    for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                }
            }
            else if (going==-1)
            {
                sp = (*X<0.0); X+=2; *Z++ = sm = 0;
                for (size_t l=1u; l<L; ++l, X+=2, ++Z) { s = (*X<0.0); sm += s*(s!=sp); *Z = sm; sp = s; }
                for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                y = (double)(*Z) / Lwf;
                for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (double)(*Z) / Lwf; }
                Z -= L;
                if (!causal)
                {
                    const size_t L2 = Lw - Lw/2u - 1u;
                    Y -= L;
                    for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                    y = *(Y-1);
                    for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                }
            }
            else { fprintf(stderr,"error in zcr_z: going must be in {-1,0,1}\n"); return 1; }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                if (going==0)
                {
                    for (size_t v=0u; v<V; ++v)
                    {
                        sp = (*X<0.0); X+=2; *Z++ = sm = 0;
                        for (size_t l=1u; l<L; ++l, X+=2, ++Z) { s = (*X<0.0); sm += (s!=sp); *Z = sm; sp = s; }
                        for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                        y = (double)(*Z) / Lwf;
                        for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                        for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (double)(*Z) / Lwf; }
                        Z -= L;
                        if (!causal)
                        {
                            const size_t L2 = Lw - Lw/2u - 1u;
                            Y -= L;
                            for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                            y = *(Y-1);
                            for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                        }
                    }
                }
                else if (going==1)
                {
                    for (size_t v=0u; v<V; ++v)
                    {
                        sp = (*X>=0.0); X+=2; *Z++ = sm = 0;
                        for (size_t l=1u; l<L; ++l, X+=2, ++Z) { s = (*X>=0.0); sm += s*(s!=sp); *Z = sm; sp = s; }
                        for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                        y = (double)(*Z) / Lwf;
                        for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                        for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (double)(*Z) / Lwf; }
                        Z -= L;
                        if (!causal)
                        {
                            const size_t L2 = Lw - Lw/2u - 1u;
                            Y -= L;
                            for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                            y = *(Y-1);
                            for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                        }
                    }
                }
                else if (going==-1)
                {
                    for (size_t v=0u; v<V; ++v)
                    {
                        sp = (*X<0.0); X+=2; *Z++ = sm = 0;
                        for (size_t l=1u; l<L; ++l, X+=2, ++Z) { s = (*X<0.0); sm += s*(s!=sp); *Z = sm; sp = s; }
                        for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                        y = (double)(*Z) / Lwf;
                        for (size_t l=0u; l<Lw-1u; ++l, ++Y) { *Y = y; }
                        for (size_t l=Lw-1u; l<L; ++l, ++Y, ++Z) { *Y = (double)(*Z) / Lwf; }
                        Z -= L;
                        if (!causal)
                        {
                            const size_t L2 = Lw - Lw/2u - 1u;
                            Y -= L;
                            for (size_t l=0u; l<L-L2; ++l, ++Y) { *Y = *(Y+L2); }
                            y = *(Y-1);
                            for (size_t l=0u; l<L2; ++l, ++Y) { *Y = y; }
                        }
                    }
                }
                else { fprintf(stderr,"error in zcr_z: going must be in {-1,0,1}\n"); return 1; }
            }
            else
            {
                if (going==0)
                {
                    for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u), Y+=B*(L-1u))
                    {
                        for (size_t b=0u; b<B; ++b, X-=2u*K*L-2u, Y-=K*L-1u)
                        {
                            sp = (*X<0.0); X+=2u*K; *Z++ = sm = 0;
                            for (size_t l=1u; l<L; ++l, X+=2u*K, ++Z) { s = (*X<0.0); sm += (s!=sp); *Z = sm; sp = s; }
                            for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                            y = (double)(*Z) / Lwf;
                            for (size_t l=0u; l<Lw-1u; ++l, Y+=K) { *Y = y; }
                            for (size_t l=Lw-1u; l<L; ++l, Y+=K, ++Z) { *Y = (double)(*Z) / Lwf; }
                            Z -= L;
                            if (!causal)
                            {
                                const size_t L2 = Lw - Lw/2u - 1u;
                                Y -= K*L;
                                for (size_t l=0u; l<L-L2; ++l, Y+=K) { *Y = *(Y+K*L2); }
                                y = *(Y-K);
                                for (size_t l=0u; l<L2; ++l, Y+=K) { *Y = y; }
                            }
                        }
                    }
                }
                else if (going==1)
                {
                    for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u), Y+=B*(L-1u))
                    {
                        for (size_t b=0u; b<B; ++b, X-=2u*K*L-2u, Y-=K*L-1u)
                        {
                            sp = (*X>=0.0); X+=2u*K; *Z++ = sm = 0;
                            for (size_t l=1u; l<L; ++l, X+=2u*K, ++Z) { s = (*X>=0.0); sm += s*(s!=sp); *Z = sm; sp = s; }
                            for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                            y = (double)(*Z) / Lwf;
                            for (size_t l=0u; l<Lw-1u; ++l, Y+=K) { *Y = y; }
                            for (size_t l=Lw-1u; l<L; ++l, Y+=K, ++Z) { *Y = (double)(*Z) / Lwf; }
                            Z -= L;
                            if (!causal)
                            {
                                const size_t L2 = Lw - Lw/2u - 1u;
                                Y -= K*L;
                                for (size_t l=0u; l<L-L2; ++l, Y+=K) { *Y = *(Y+K*L2); }
                                y = *(Y-K);
                                for (size_t l=0u; l<L2; ++l, Y+=K) { *Y = y; }
                            }
                        }
                    }
                }
                else if (going==-1)
                {
                    for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u), Y+=B*(L-1u))
                    {
                        for (size_t b=0u; b<B; ++b, X-=2u*K*L-2u, Y-=K*L-1u)
                        {
                            sp = (*X<0.0); X+=2u*K; *Z++ = sm = 0;
                            for (size_t l=1u; l<L; ++l, X+=2u*K, ++Z) { s = (*X<0.0); sm += s*(s!=sp); *Z = sm; sp = s; }
                            for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                            y = (double)(*Z) / Lwf;
                            for (size_t l=0u; l<Lw-1u; ++l, Y+=K) { *Y = y; }
                            for (size_t l=Lw-1u; l<L; ++l, Y+=K, ++Z) { *Y = (double)(*Z) / Lwf; }
                            Z -= L;
                            if (!causal)
                            {
                                const size_t L2 = Lw - Lw/2u - 1u;
                                Y -= K*L;
                                for (size_t l=0u; l<L-L2; ++l, Y+=K) { *Y = *(Y+K*L2); }
                                y = *(Y-K);
                                for (size_t l=0u; l<L2; ++l, Y+=K) { *Y = y; }
                            }
                        }
                    }
                }
                else { fprintf(stderr,"error in zcr_z: going must be in {-1,0,1}\n"); return 1; }
            }
        }
        free(Z);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
