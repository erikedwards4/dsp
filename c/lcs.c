//Level crossings, with Boolean (int) output.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int lcs_s (int *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int going, const float lvl);
int lcs_d (int *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int going, const double lvl);


int lcs_s (int *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int going, const float lvl)
{
    if (dim>3u) { fprintf(stderr,"error in lcs_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int s, sp; //sign info and previous sign info

    if (N==0u) {}
    else if (L==N)
    {
        if (going==0)
        {
            sp = (*X++<lvl); *Y++ = 0;
            for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X<lvl); *Y = (s!=sp); sp = s; }
        }
        else if (going==1)
        {
            sp = (*X++>=lvl); *Y++ = 0;
            for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X>=lvl); *Y = s*(s!=sp); sp = s; }
        }
        else if (going==-1)
        {
            sp = (*X++<lvl); *Y++ = 0;
            for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X<lvl); *Y = s*(s!=sp); sp = s; }
        }
        else { fprintf(stderr,"error in lcs_s: going must be in {-1,0,1}\n"); return 1; }
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
                    sp = (*X++<lvl); *Y++ = 0;
                    for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X<lvl); *Y = (s!=sp); sp = s; }
                }
            }
            else if (going==1)
            {
                for (size_t v=V; v>0u; --v)
                {
                    sp = (*X++>=lvl); *Y++ = 0;
                    for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X>=lvl); *Y = s*(s!=sp); sp = s; }
                }
            }
            else if (going==-1)
            {
                for (size_t v=V; v>0u; --v)
                {
                    sp = (*X++<lvl); *Y++ = 0;
                    for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X<lvl); *Y = s*(s!=sp); sp = s; }
                }
            }
            else { fprintf(stderr,"error in lcs_s: going must be in {-1,0,1}\n"); return 1; }
        }
        else
        {
            if (going==0)
            {
                for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                    {
                        sp = (*X<lvl); *Y = 0; X+=K; Y+=K;
                        for (size_t l=1u; l<L; ++l, X+=K, Y+=K) { s = (*X<lvl); *Y = (s!=sp); sp = s; }
                    }
                }
            }
            else if (going==1)
            {
                for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                    {
                        sp = (*X>=lvl); *Y = 0; X+=K; Y+=K;
                        for (size_t l=1u; l<L; ++l, X+=K, Y+=K) { s = (*X>=lvl); *Y = s*(s!=sp); sp = s; }
                    }
                }
            }
            else if (going==-1)
            {
                for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                    {
                        sp = (*X<lvl); *Y = 0; X+=K; Y+=K;
                        for (size_t l=1u; l<L; ++l, X+=K, Y+=K) { s = (*X<lvl); *Y = s*(s!=sp); sp = s; }
                    }
                }
            }
            else { fprintf(stderr,"error in lcs_s: going must be in {-1,0,1}\n"); return 1; }
        }
    }

    return 0;
}


int lcs_d (int *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int going, const double lvl)
{
    if (dim>3u) { fprintf(stderr,"error in lcs_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int s, sp; //sign info and previous sign info

    if (N==0u) {}
    else if (L==N)
    {
        if (going==0)
        {
            sp = (*X++<lvl); *Y++ = 0;
            for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X<lvl); *Y = (s!=sp); sp = s; }
        }
        else if (going==1)
        {
            sp = (*X++>=lvl); *Y++ = 0;
            for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X>=lvl); *Y = s*(s!=sp); sp = s; }
        }
        else if (going==-1)
        {
            sp = (*X++<lvl); *Y++ = 0;
            for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X<lvl); *Y = s*(s!=sp); sp = s; }
        }
        else { fprintf(stderr,"error in lcs_d: going must be in {-1,0,1}\n"); return 1; }
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
                    sp = (*X++<lvl); *Y++ = 0;
                    for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X<lvl); *Y = (s!=sp); sp = s; }
                }
            }
            else if (going==1)
            {
                for (size_t v=V; v>0u; --v)
                {
                    sp = (*X++>=lvl); *Y++ = 0;
                    for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X>=lvl); *Y = s*(s!=sp); sp = s; }
                }
            }
            else if (going==-1)
            {
                for (size_t v=V; v>0u; --v)
                {
                    sp = (*X++<lvl); *Y++ = 0;
                    for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X<lvl); *Y = s*(s!=sp); sp = s; }
                }
            }
            else { fprintf(stderr,"error in lcs_d: going must be in {-1,0,1}\n"); return 1; }
        }
        else
        {
            if (going==0)
            {
                for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                    {
                        sp = (*X<lvl); *Y = 0; X+=K; Y+=K;
                        for (size_t l=1u; l<L; ++l, X+=K, Y+=K) { s = (*X<lvl); *Y = (s!=sp); sp = s; }
                    }
                }
            }
            else if (going==1)
            {
                for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                    {
                        sp = (*X>=lvl); *Y = 0; X+=K; Y+=K;
                        for (size_t l=1u; l<L; ++l, X+=K, Y+=K) { s = (*X>=lvl); *Y = s*(s!=sp); sp = s; }
                    }
                }
            }
            else if (going==-1)
            {
                for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                    {
                        sp = (*X<lvl); *Y = 0; X+=K; Y+=K;
                        for (size_t l=1u; l<L; ++l, X+=K, Y+=K) { s = (*X<lvl); *Y = s*(s!=sp); sp = s; }
                    }
                }
            }
            else { fprintf(stderr,"error in lcs_d: going must be in {-1,0,1}\n"); return 1; }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
