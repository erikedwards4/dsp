//Gets the level-crossing rate (LCR) for each vec in X along dim.
//If lvl==0.0, then the output is identical to zcr.
//See zcr.c for general comments and info.

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int lcr_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const float lvl, const int causal);
int lcr_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const double lvl, const int causal);


int lcr_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const float lvl, const int causal)
{
    if (dim>3u) { fprintf(stderr,"error in lcr_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lw>L) { fprintf(stderr,"error in lcr_s: Lw (winlength) must be <= L (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        const float Lwf = (float)Lw;
        float y;
        int *Z, s, sp, sm;
        if (!(Z=(int *)malloc(L*sizeof(int)))) { fprintf(stderr,"error in lcr_s: problem with malloc. "); perror("malloc"); return 1; }

        if (L==N)
        {
            if (going==0)
            {
                sp = (*X++<lvl); *Z++ = sm = 0;
                for (size_t l=L; l>1u; --l, ++X, ++Z) { s = (*X<lvl); sm += (s!=sp); *Z = sm; sp = s; }
                for (size_t l=L-Lw+1u; l>0u; --l) { --Z; *Z -= *(Z-Lw); }
                y = (float)(*Z) / Lwf;
                for (size_t l=Lw-1u; l>0u; --l, ++Y) { *Y = y; }
                for (size_t l=L-Lw+1u; l>0u; --l, ++Y, ++Z) { *Y = (float)(*Z) / Lwf; }
                Z -= L;
                if (!causal)
                {
                    const size_t L2 = Lw - Lw/2u - 1u;
                    Y -= L;
                    for (size_t l=L-L2; l>0u; --l, ++Y) { *Y = *(Y+L2); }
                    y = *(Y-1);
                    for (size_t l=L2; l>0u; --l, ++Y) { *Y = y; }
                }
            }
            else if (going==1)
            {
                sp = (*X++>=lvl); *Z++ = sm = 0;
                for (size_t l=L; l>1u; --l, ++X, ++Z) { s = (*X>=lvl); sm += s*(s!=sp); *Z = sm; sp = s; }
                for (size_t l=L-Lw+1u; l>0u; --l) { --Z; *Z -= *(Z-Lw); }
                y = (float)(*Z) / Lwf;
                for (size_t l=Lw-1u; l>0u; --l, ++Y) { *Y = y; }
                for (size_t l=L-Lw+1u; l>0u; --l, ++Y, ++Z) { *Y = (float)(*Z) / Lwf; }
                Z -= L;
                if (!causal)
                {
                    const size_t L2 = Lw - Lw/2u - 1u;
                    Y -= L;
                    for (size_t l=L-L2; l>0u; --l, ++Y) { *Y = *(Y+L2); }
                    y = *(Y-1);
                    for (size_t l=L2; l>0u; --l, ++Y) { *Y = y; }
                }
            }
            else if (going==-1)
            {
                sp = (*X++<lvl); *Z++ = sm = 0;
                for (size_t l=L; l>1u; --l, ++X, ++Z) { s = (*X<lvl); sm += s*(s!=sp); *Z = sm; sp = s; }
                for (size_t l=L-Lw+1u; l>0u; --l) { --Z; *Z -= *(Z-Lw); }
                y = (float)(*Z) / Lwf;
                for (size_t l=Lw-1u; l>0u; --l, ++Y) { *Y = y; }
                for (size_t l=L-Lw+1u; l>0u; --l, ++Y, ++Z) { *Y = (float)(*Z) / Lwf; }
                Z -= L;
                if (!causal)
                {
                    const size_t L2 = Lw - Lw/2u - 1u;
                    Y -= L;
                    for (size_t l=L-L2; l>0u; --l, ++Y) { *Y = *(Y+L2); }
                    y = *(Y-1);
                    for (size_t l=L2; l>0u; --l, ++Y) { *Y = y; }
                }
            }
            else { fprintf(stderr,"error in lcr_s: going must be in {-1,0,1}\n"); return 1; }
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
                        sp = (*X++<lvl); *Z++ = sm = 0;
                        for (size_t l=L; l>1u; --l, ++X, ++Z) { s = (*X<lvl); sm += (s!=sp); *Z = sm; sp = s; }
                        for (size_t l=L-Lw+1u; l>0u; --l) { --Z; *Z -= *(Z-Lw); }
                        y = (float)(*Z) / Lwf;
                        for (size_t l=Lw-1u; l>0u; --l, ++Y) { *Y = y; }
                        for (size_t l=L-Lw+1u; l>0u; --l, ++Y, ++Z) { *Y = (float)(*Z) / Lwf; }
                        Z -= L;
                        if (!causal)
                        {
                            const size_t L2 = Lw - Lw/2u - 1u;
                            Y -= L;
                            for (size_t l=L-L2; l>0u; --l, ++Y) { *Y = *(Y+L2); }
                            y = *(Y-1);
                            for (size_t l=L2; l>0u; --l, ++Y) { *Y = y; }
                        }
                    }
                }
                else if (going==1)
                {
                    for (size_t v=V; v>0u; --v)
                    {
                        sp = (*X++>=lvl); *Z++ = sm = 0;
                        for (size_t l=L; l>1u; --l, ++X, ++Z) { s = (*X>=lvl); sm += s*(s!=sp); *Z = sm; sp = s; }
                        for (size_t l=L-Lw+1u; l>0u; --l) { --Z; *Z -= *(Z-Lw); }
                        y = (float)(*Z) / Lwf;
                        for (size_t l=Lw-1u; l>0u; --l, ++Y) { *Y = y; }
                        for (size_t l=L-Lw+1u; l>0u; --l, ++Y, ++Z) { *Y = (float)(*Z) / Lwf; }
                        Z -= L;
                        if (!causal)
                        {
                            const size_t L2 = Lw - Lw/2u - 1u;
                            Y -= L;
                            for (size_t l=L-L2; l>0u; --l, ++Y) { *Y = *(Y+L2); }
                            y = *(Y-1);
                            for (size_t l=L2; l>0u; --l, ++Y) { *Y = y; }
                        }
                    }
                }
                else if (going==-1)
                {
                    for (size_t v=V; v>0u; --v)
                    {
                        sp = (*X++<lvl); *Z++ = sm = 0;
                        for (size_t l=L; l>1u; --l, ++X, ++Z) { s = (*X<lvl); sm += s*(s!=sp); *Z = sm; sp = s; }
                        for (size_t l=L-Lw+1u; l>0u; --l) { --Z; *Z -= *(Z-Lw); }
                        y = (float)(*Z) / Lwf;
                        for (size_t l=Lw-1u; l>0u; --l, ++Y) { *Y = y; }
                        for (size_t l=L-Lw+1u; l>0u; --l, ++Y, ++Z) { *Y = (float)(*Z) / Lwf; }
                        Z -= L;
                        if (!causal)
                        {
                            const size_t L2 = Lw - Lw/2u - 1u;
                            Y -= L;
                            for (size_t l=L-L2; l>0u; --l, ++Y) { *Y = *(Y+L2); }
                            y = *(Y-1);
                            for (size_t l=L2; l>0u; --l, ++Y) { *Y = y; }
                        }
                    }
                }
                else { fprintf(stderr,"error in lcr_s: going must be in {-1,0,1}\n"); return 1; }
            }
            else
            {
                if (going==0)
                {
                    for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
                    {
                        for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                        {
                            sp = (*X<lvl); X+=K; *Z++ = sm = 0;
                            for (size_t l=L; l>1u; --l, X+=K, ++Z) { s = (*X<lvl); sm += (s!=sp); *Z = sm; sp = s; }
                            for (size_t l=L-Lw+1u; l>0u; --l) { --Z; *Z -= *(Z-Lw); }
                            y = (float)(*Z) / Lwf;
                            for (size_t l=Lw-1u; l>0u; --l, Y+=K) { *Y = y; }
                            for (size_t l=L-Lw+1u; l>0u; --l, Y+=K, ++Z) { *Y = (float)(*Z) / Lwf; }
                            Z -= L;
                            if (!causal)
                            {
                                const size_t L2 = Lw - Lw/2u - 1u;
                                Y -= K*L;
                                for (size_t l=L-L2; l>0u; --l, Y+=K) { *Y = *(Y+K*L2); }
                                y = *(Y-K);
                                for (size_t l=L2; l>0u; --l, Y+=K) { *Y = y; }
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
                            sp = (*X>=lvl); X+=K; *Z++ = sm = 0;
                            for (size_t l=L; l>1u; --l, X+=K, ++Z) { s = (*X>=lvl); sm += s*(s!=sp); *Z = sm; sp = s; }
                            for (size_t l=L-Lw+1u; l>0u; --l) { --Z; *Z -= *(Z-Lw); }
                            y = (float)(*Z) / Lwf;
                            for (size_t l=Lw-1u; l>0u; --l, Y+=K) { *Y = y; }
                            for (size_t l=L-Lw+1u; l>0u; --l, Y+=K, ++Z) { *Y = (float)(*Z) / Lwf; }
                            Z -= L;
                            if (!causal)
                            {
                                const size_t L2 = Lw - Lw/2u - 1u;
                                Y -= K*L;
                                for (size_t l=L-L2; l>0u; --l, Y+=K) { *Y = *(Y+K*L2); }
                                y = *(Y-K);
                                for (size_t l=L2; l>0u; --l, Y+=K) { *Y = y; }
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
                            sp = (*X<lvl); X+=K; *Z++ = sm = 0;
                            for (size_t l=L; l>1u; --l, X+=K, ++Z) { s = (*X<lvl); sm += s*(s!=sp); *Z = sm; sp = s; }
                            for (size_t l=L-Lw+1u; l>0u; --l) { --Z; *Z -= *(Z-Lw); }
                            y = (float)(*Z) / Lwf;
                            for (size_t l=Lw-1u; l>0u; --l, Y+=K) { *Y = y; }
                            for (size_t l=L-Lw+1u; l>0u; --l, Y+=K, ++Z) { *Y = (float)(*Z) / Lwf; }
                            Z -= L;
                            if (!causal)
                            {
                                const size_t L2 = Lw - Lw/2u - 1u;
                                Y -= K*L;
                                for (size_t l=L-L2; l>0u; --l, Y+=K) { *Y = *(Y+K*L2); }
                                y = *(Y-K);
                                for (size_t l=L2; l>0u; --l, Y+=K) { *Y = y; }
                            }
                        }
                    }
                }
                else { fprintf(stderr,"error in lcr_s: going must be in {-1,0,1}\n"); return 1; }
            }
        }
        free(Z);
    }

    return 0;
}


int lcr_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const double lvl, const int causal)
{
    if (dim>3u) { fprintf(stderr,"error in lcr_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lw>L) { fprintf(stderr,"error in lcr_d: Lw (winlength) must be <= L (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        const double Lwf = (double)Lw;
        double y;
        int *Z, s, sp, sm;
        if (!(Z=(int *)malloc(L*sizeof(int)))) { fprintf(stderr,"error in lcr_d: problem with malloc. "); perror("malloc"); return 1; }

        if (L==N)
        {
            if (going==0)
            {
                sp = (*X++<lvl); *Z++ = sm = 0;
                for (size_t l=L; l>1u; --l, ++X, ++Z) { s = (*X<lvl); sm += (s!=sp); *Z = sm; sp = s; }
                for (size_t l=L-Lw+1u; l>0u; --l) { --Z; *Z -= *(Z-Lw); }
                y = (double)(*Z) / Lwf;
                for (size_t l=Lw-1u; l>0u; --l, ++Y) { *Y = y; }
                for (size_t l=L-Lw+1u; l>0u; --l, ++Y, ++Z) { *Y = (double)(*Z) / Lwf; }
                Z -= L;
                if (!causal)
                {
                    const size_t L2 = Lw - Lw/2u - 1u;
                    Y -= L;
                    for (size_t l=L-L2; l>0u; --l, ++Y) { *Y = *(Y+L2); }
                    y = *(Y-1);
                    for (size_t l=L2; l>0u; --l, ++Y) { *Y = y; }
                }
            }
            else if (going==1)
            {
                sp = (*X++>=lvl); *Z++ = sm = 0;
                for (size_t l=L; l>1u; --l, ++X, ++Z) { s = (*X>=lvl); sm += s*(s!=sp); *Z = sm; sp = s; }
                for (size_t l=L-Lw+1u; l>0u; --l) { --Z; *Z -= *(Z-Lw); }
                y = (double)(*Z) / Lwf;
                for (size_t l=Lw-1u; l>0u; --l, ++Y) { *Y = y; }
                for (size_t l=L-Lw+1u; l>0u; --l, ++Y, ++Z) { *Y = (double)(*Z) / Lwf; }
                Z -= L;
                if (!causal)
                {
                    const size_t L2 = Lw - Lw/2u - 1u;
                    Y -= L;
                    for (size_t l=L-L2; l>0u; --l, ++Y) { *Y = *(Y+L2); }
                    y = *(Y-1);
                    for (size_t l=L2; l>0u; --l, ++Y) { *Y = y; }
                }
            }
            else if (going==-1)
            {
                sp = (*X++<lvl); *Z++ = sm = 0;
                for (size_t l=L; l>1u; --l, ++X, ++Z) { s = (*X<lvl); sm += s*(s!=sp); *Z = sm; sp = s; }
                for (size_t l=L-Lw+1u; l>0u; --l) { --Z; *Z -= *(Z-Lw); }
                y = (double)(*Z) / Lwf;
                for (size_t l=Lw-1u; l>0u; --l, ++Y) { *Y = y; }
                for (size_t l=L-Lw+1u; l>0u; --l, ++Y, ++Z) { *Y = (double)(*Z) / Lwf; }
                Z -= L;
                if (!causal)
                {
                    const size_t L2 = Lw - Lw/2u - 1u;
                    Y -= L;
                    for (size_t l=L-L2; l>0u; --l, ++Y) { *Y = *(Y+L2); }
                    y = *(Y-1);
                    for (size_t l=L2; l>0u; --l, ++Y) { *Y = y; }
                }
            }
            else { fprintf(stderr,"error in lcr_d: going must be in {-1,0,1}\n"); return 1; }
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
                        sp = (*X++<lvl); *Z++ = sm = 0;
                        for (size_t l=L; l>1u; --l, ++X, ++Z) { s = (*X<lvl); sm += (s!=sp); *Z = sm; sp = s; }
                        for (size_t l=L-Lw+1u; l>0u; --l) { --Z; *Z -= *(Z-Lw); }
                        y = (double)(*Z) / Lwf;
                        for (size_t l=Lw-1u; l>0u; --l, ++Y) { *Y = y; }
                        for (size_t l=L-Lw+1u; l>0u; --l, ++Y, ++Z) { *Y = (double)(*Z) / Lwf; }
                        Z -= L;
                        if (!causal)
                        {
                            const size_t L2 = Lw - Lw/2u - 1u;
                            Y -= L;
                            for (size_t l=L-L2; l>0u; --l, ++Y) { *Y = *(Y+L2); }
                            y = *(Y-1);
                            for (size_t l=L2; l>0u; --l, ++Y) { *Y = y; }
                        }
                    }
                }
                else if (going==1)
                {
                    for (size_t v=V; v>0u; --v)
                    {
                        sp = (*X++>=lvl); *Z++ = sm = 0;
                        for (size_t l=L; l>1u; --l, ++X, ++Z) { s = (*X>=lvl); sm += s*(s!=sp); *Z = sm; sp = s; }
                        for (size_t l=L-Lw+1u; l>0u; --l) { --Z; *Z -= *(Z-Lw); }
                        y = (double)(*Z) / Lwf;
                        for (size_t l=Lw-1u; l>0u; --l, ++Y) { *Y = y; }
                        for (size_t l=L-Lw+1u; l>0u; --l, ++Y, ++Z) { *Y = (double)(*Z) / Lwf; }
                        Z -= L;
                        if (!causal)
                        {
                            const size_t L2 = Lw - Lw/2u - 1u;
                            Y -= L;
                            for (size_t l=L-L2; l>0u; --l, ++Y) { *Y = *(Y+L2); }
                            y = *(Y-1);
                            for (size_t l=L2; l>0u; --l, ++Y) { *Y = y; }
                        }
                    }
                }
                else if (going==-1)
                {
                    for (size_t v=V; v>0u; --v)
                    {
                        sp = (*X++<lvl); *Z++ = sm = 0;
                        for (size_t l=L; l>1u; --l, ++X, ++Z) { s = (*X<lvl); sm += s*(s!=sp); *Z = sm; sp = s; }
                        for (size_t l=L-Lw+1u; l>0u; --l) { --Z; *Z -= *(Z-Lw); }
                        y = (double)(*Z) / Lwf;
                        for (size_t l=Lw-1u; l>0u; --l, ++Y) { *Y = y; }
                        for (size_t l=L-Lw+1u; l>0u; --l, ++Y, ++Z) { *Y = (double)(*Z) / Lwf; }
                        Z -= L;
                        if (!causal)
                        {
                            const size_t L2 = Lw - Lw/2u - 1u;
                            Y -= L;
                            for (size_t l=L-L2; l>0u; --l, ++Y) { *Y = *(Y+L2); }
                            y = *(Y-1);
                            for (size_t l=L2; l>0u; --l, ++Y) { *Y = y; }
                        }
                    }
                }
                else { fprintf(stderr,"error in lcr_d: going must be in {-1,0,1}\n"); return 1; }
            }
            else
            {
                if (going==0)
                {
                    for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
                    {
                        for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                        {
                            sp = (*X<lvl); X+=K; *Z++ = sm = 0;
                            for (size_t l=L; l>1u; --l, X+=K, ++Z) { s = (*X<lvl); sm += (s!=sp); *Z = sm; sp = s; }
                            for (size_t l=L-Lw+1u; l>0u; --l) { --Z; *Z -= *(Z-Lw); }
                            y = (double)(*Z) / Lwf;
                            for (size_t l=Lw-1u; l>0u; --l, Y+=K) { *Y = y; }
                            for (size_t l=L-Lw+1u; l>0u; --l, Y+=K, ++Z) { *Y = (double)(*Z) / Lwf; }
                            Z -= L;
                            if (!causal)
                            {
                                const size_t L2 = Lw - Lw/2u - 1u;
                                Y -= K*L;
                                for (size_t l=L-L2; l>0u; --l, Y+=K) { *Y = *(Y+K*L2); }
                                y = *(Y-K);
                                for (size_t l=L2; l>0u; --l, Y+=K) { *Y = y; }
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
                            sp = (*X>=lvl); X+=K; *Z++ = sm = 0;
                            for (size_t l=L; l>1u; --l, X+=K, ++Z) { s = (*X>=lvl); sm += s*(s!=sp); *Z = sm; sp = s; }
                            for (size_t l=L-Lw+1u; l>0u; --l) { --Z; *Z -= *(Z-Lw); }
                            y = (double)(*Z) / Lwf;
                            for (size_t l=Lw-1u; l>0u; --l, Y+=K) { *Y = y; }
                            for (size_t l=L-Lw+1u; l>0u; --l, Y+=K, ++Z) { *Y = (double)(*Z) / Lwf; }
                            Z -= L;
                            if (!causal)
                            {
                                const size_t L2 = Lw - Lw/2u - 1u;
                                Y -= K*L;
                                for (size_t l=L-L2; l>0u; --l, Y+=K) { *Y = *(Y+K*L2); }
                                y = *(Y-K);
                                for (size_t l=L2; l>0u; --l, Y+=K) { *Y = y; }
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
                            sp = (*X<lvl); X+=K; *Z++ = sm = 0;
                            for (size_t l=L; l>1u; --l, X+=K, ++Z) { s = (*X<lvl); sm += s*(s!=sp); *Z = sm; sp = s; }
                            for (size_t l=L-Lw+1u; l>0u; --l) { --Z; *Z -= *(Z-Lw); }
                            y = (double)(*Z) / Lwf;
                            for (size_t l=Lw-1u; l>0u; --l, Y+=K) { *Y = y; }
                            for (size_t l=L-Lw+1u; l>0u; --l, Y+=K, ++Z) { *Y = (double)(*Z) / Lwf; }
                            Z -= L;
                            if (!causal)
                            {
                                const size_t L2 = Lw - Lw/2u - 1u;
                                Y -= K*L;
                                for (size_t l=L-L2; l>0u; --l, Y+=K) { *Y = *(Y+K*L2); }
                                y = *(Y-K);
                                for (size_t l=L2; l>0u; --l, Y+=K) { *Y = y; }
                            }
                        }
                    }
                }
                else { fprintf(stderr,"error in lcr_d: going must be in {-1,0,1}\n"); return 1; }
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
