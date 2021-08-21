//This gets MCs as usual, and then averages over (rectangular) window of length Lw.
//See mcs.c and zcr.c for general comments and info.

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int mcr_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const int causal);
int mcr_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const int causal);


int mcr_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const int causal)
{
    if (dim>3u) { fprintf(stderr,"error in mcr_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lw>L) { fprintf(stderr,"error in mcr_s: Lw (winlength) must be <= L (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        const float Lwf = (float)Lw;
        float y, mn;
        int *Z, s, sp, sm;
        if (!(Z=(int *)malloc(L*sizeof(int)))) { fprintf(stderr,"error in mcr_s: problem with malloc. "); perror("malloc"); return 1; }

        if (L==N)
        {
            //Mean
            mn = 0.0f;
            for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
            mn /= (float)L; X -= L;

            //MCR
            if (going==0)
            {
                sp = (*X++<mn); *Z++ = sm = 0;
                for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X<mn); sm += (s!=sp); *Z = sm; sp = s; }
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
                sp = (*X++>=mn); *Z++ = sm = 0;
                for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X>=mn); sm += s*(s!=sp); *Z = sm; sp = s; }
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
                sp = (*X++<mn); *Z++ = sm = 0;
                for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X<mn); sm += s*(s!=sp); *Z = sm; sp = s; }
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
            else { fprintf(stderr,"error in mcr_s: going must be in {-1,0,1}\n"); return 1; }
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
                    for (size_t v=V; v>0u; --v)
                    {
                        mn = 0.0f;
                        for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
                        mn /= (float)L; X -= L;

                        sp = (*X++<mn); *Z++ = sm = 0;
                        for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X<mn); sm += (s!=sp); *Z = sm; sp = s; }
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
                    for (size_t v=V; v>0u; --v)
                    {
                        mn = 0.0f;
                        for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
                        mn /= (float)L; X -= L;

                        sp = (*X++>=mn); *Z++ = sm = 0;
                        for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X>=mn); sm += s*(s!=sp); *Z = sm; sp = s; }
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
                    for (size_t v=V; v>0u; --v)
                    {
                        mn = 0.0f;
                        for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
                        mn /= (float)L; X -= L;

                        sp = (*X++<mn); *Z++ = sm = 0;
                        for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X<mn); sm += s*(s!=sp); *Z = sm; sp = s; }
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
                else { fprintf(stderr,"error in mcr_s: going must be in {-1,0,1}\n"); return 1; }
            }
            else
            {
                if (going==0)
                {
                    for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
                    {
                        for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                        {
                            mn = 0.0f;
                            for (size_t l=0u; l<L; ++l, X+=K) { mn += *X; }
                            mn /= (float)L; X -= K*L;
                            sp = (*X<mn); X+=K; *Z++ = sm = 0;
                            for (size_t l=1u; l<L; ++l, X+=K, ++Z) { s = (*X<mn); sm += (s!=sp); *Z = sm; sp = s; }
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
                    for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
                    {
                        for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                        {
                            mn = 0.0f;
                            for (size_t l=0u; l<L; ++l, X+=K) { mn += *X; }
                            mn /= (float)L; X -= K*L;
                            sp = (*X>=mn); X+=K; *Z++ = sm = 0;
                            for (size_t l=1u; l<L; ++l, X+=K, ++Z) { s = (*X>=mn); sm += s*(s!=sp); *Z = sm; sp = s; }
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
                    for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
                    {
                        for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                        {
                            mn = 0.0f;
                            for (size_t l=0u; l<L; ++l, X+=K) { mn += *X; }
                            mn /= (float)L; X -= K*L;
                            sp = (*X<mn); X+=K; *Z++ = sm = 0;
                            for (size_t l=1u; l<L; ++l, X+=K, ++Z) { s = (*X<mn); sm += s*(s!=sp); *Z = sm; sp = s; }
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
                else { fprintf(stderr,"error in mcr_s: going must be in {-1,0,1}\n"); return 1; }
            }
        }
        free(Z);
    }

    return 0;
}


int mcr_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const int causal)
{
    if (dim>3u) { fprintf(stderr,"error in mcr_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lw>L) { fprintf(stderr,"error in mcr_d: Lw (winlength) must be <= L (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        const double Lwf = (double)Lw;
        double y, mn;
        int *Z, s, sp, sm;
        if (!(Z=(int *)malloc(L*sizeof(int)))) { fprintf(stderr,"error in mcr_d: problem with malloc. "); perror("malloc"); return 1; }

        if (L==N)
        {
            mn = 0.0;
            for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
            mn /= (double)L; X -= L;

            if (going==0)
            {
                sp = (*X++<mn); *Z++ = sm = 0;
                for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X<mn); sm += (s!=sp); *Z = sm; sp = s; }
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
                sp = (*X++>=mn); *Z++ = sm = 0;
                for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X>=mn); sm += s*(s!=sp); *Z = sm; sp = s; }
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
                sp = (*X++<mn); *Z++ = sm = 0;
                for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X<mn); sm += s*(s!=sp); *Z = sm; sp = s; }
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
            else { fprintf(stderr,"error in mcr_d: going must be in {-1,0,1}\n"); return 1; }
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
                    for (size_t v=V; v>0u; --v)
                    {
                        mn = 0.0;
                        for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
                        mn /= (double)L; X -= L;

                        sp = (*X++<mn); *Z++ = sm = 0;
                        for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X<mn); sm += (s!=sp); *Z = sm; sp = s; }
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
                    for (size_t v=V; v>0u; --v)
                    {
                        mn = 0.0;
                        for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
                        mn /= (double)L; X -= L;

                        sp = (*X++>=mn); *Z++ = sm = 0;
                        for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X>=mn); sm += s*(s!=sp); *Z = sm; sp = s; }
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
                    for (size_t v=V; v>0u; --v)
                    {
                        mn = 0.0;
                        for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
                        mn /= (double)L; X -= L;

                        sp = (*X++<mn); *Z++ = sm = 0;
                        for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X<mn); sm += s*(s!=sp); *Z = sm; sp = s; }
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
                else { fprintf(stderr,"error in mcr_d: going must be in {-1,0,1}\n"); return 1; }
            }
            else
            {
                if (going==0)
                {
                    for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
                    {
                        for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                        {
                            mn = 0.0;
                            for (size_t l=0u; l<L; ++l, X+=K) { mn += *X; }
                            mn /= (double)L; X -= K*L;
                            sp = (*X<mn); X+=K; *Z++ = sm = 0;
                            for (size_t l=1u; l<L; ++l, X+=K, ++Z) { s = (*X<mn); sm += (s!=sp); *Z = sm; sp = s; }
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
                    for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
                    {
                        for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                        {
                            mn = 0.0;
                            for (size_t l=0u; l<L; ++l, X+=K) { mn += *X; }
                            mn /= (double)L; X -= K*L;
                            sp = (*X>=mn); X+=K; *Z++ = sm = 0;
                            for (size_t l=1u; l<L; ++l, X+=K, ++Z) { s = (*X>=mn); sm += s*(s!=sp); *Z = sm; sp = s; }
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
                    for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
                    {
                        for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                        {
                            mn = 0.0;
                            for (size_t l=0u; l<L; ++l, X+=K) { mn += *X; }
                            mn /= (double)L; X -= K*L;
                            sp = (*X<mn); X+=K; *Z++ = sm = 0;
                            for (size_t l=1u; l<L; ++l, X+=K, ++Z) { s = (*X<mn); sm += s*(s!=sp); *Z = sm; sp = s; }
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
                else { fprintf(stderr,"error in mcr_d: going must be in {-1,0,1}\n"); return 1; }
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
