//Gets reflection coefficients (RCs) from autoregressive (AR) coeffs for each vector in X.
//This is similar to poly2rc, but no need to normalize by X[0], since AR coeffs assume this.

//This currently has the same sign of Matlab output (see commented code below).
//However, this may not match the Octave tsa toolbox in sign convention.
//Right now, only the current sign convention allows round-trip with rc2poly.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ar2rc_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ar2rc_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ar2rc_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ar2rc_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);


int ar2rc_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ar2rc_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else
    {
        float sc, sc2, *y;
        if (!(y=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in ar2rc_s: problem with malloc. "); perror("malloc"); return 1; }
        
        if (L==N)
        {
            for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = -*X; }
            Y -= L;
            for (size_t l=L-1u; l>0u; --l)
            {
                for (size_t q=0u; q<=l; ++q, ++Y, ++y) { *y = *Y; }
                sc = *--y; Y -= l + 1u;
                for (size_t q=0u; q<l; ++q, ++Y) { --y; *Y += sc * *y; }
                sc2 = fmaf(sc,-sc,1.0f);
                //*Y = -*Y;  //to match other sign convention
                for (size_t q=0u; q<l; ++q) { --Y; *Y /= sc2; }
            }
            //*Y = -*Y;  //to match other sign convention
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, Y+=L)
                {
                    for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = -*X; }
                    Y -= L;
                    for (size_t l=L-1u; l>0u; --l)
                    {
                        for (size_t q=0u; q<=l; ++q, ++Y, ++y) { *y = *Y; }
                        sc = *--y; Y -= l + 1u;
                        for (size_t q=0u; q<l; ++q, ++Y) { --y; *Y += sc * *y; }
                        sc2 = fmaf(sc,-sc,1.0f);
                        //*Y = -*Y;  //to match other sign convention
                        for (size_t q=0u; q<l; ++q) { --Y; *Y /= sc2; }
                    }
                    //*Y = -*Y;  //to match other sign convention
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=B*(L-1u), Y+=B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=K*L-1u, ++Y)
                    {
                        for (size_t l=0u; l<L; ++l, X+=K, Y+=K) { *Y = -*X; }
                        Y -= L*K;
                        for (size_t l=L-1u; l>0u; --l)
                        {
                            for (size_t q=0u; q<=l; ++q, Y+=K, ++y) { *y = *Y; }
                            sc = *--y; Y -= K*(l+1u);
                            for (size_t q=0u; q<l; ++q, Y+=K) { --y; *Y += sc * *y; }
                            sc2 = fmaf(sc,-sc,1.0f);
                            //*Y = -*Y;  //to match other sign convention
                            for (size_t q=0u; q<l; ++q) { Y-=K; *Y /= sc2; }
                        }
                        //*Y = -*Y;  //to match other sign convention
                    }
                }
            }
        }
        free(y);
    }

    return 0;
}


int ar2rc_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ar2rc_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else
    {
        double sc, sc2, *y;
        if (!(y=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in ar2rc_d: problem with malloc. "); perror("malloc"); return 1; }
        
        if (L==N)
        {
            for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = -*X; }
            Y -= L;
            for (size_t l=L-1u; l>0u; --l)
            {
                for (size_t q=0u; q<=l; ++q, ++Y, ++y) { *y = *Y; }
                sc = *--y; Y -= l + 1u;
                for (size_t q=0u; q<l; ++q, ++Y) { --y; *Y += sc * *y; }
                sc2 = fma(sc,-sc,1.0);
                //*Y = -*Y;  //to match other sign convention
                for (size_t q=0u; q<l; ++q) { --Y; *Y /= sc2; }
            }
            //*Y = -*Y;  //to match other sign convention
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, Y+=L)
                {
                    for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = -*X; }
                    Y -= L;
                    for (size_t l=L-1u; l>0u; --l)
                    {
                        for (size_t q=0u; q<=l; ++q, ++Y, ++y) { *y = *Y; }
                        sc = *--y; Y -= l + 1u;
                        for (size_t q=0u; q<l; ++q, ++Y) { --y; *Y += sc * *y; }
                        sc2 = fma(sc,-sc,1.0);
                        //*Y = -*Y;  //to match other sign convention
                        for (size_t q=0u; q<l; ++q) { --Y; *Y /= sc2; }
                    }
                    //*Y = -*Y;  //to match other sign convention
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=B*(L-1u), Y+=B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=K*L-1u, ++Y)
                    {
                        for (size_t l=0u; l<L; ++l, X+=K, Y+=K) { *Y = -*X; }
                        Y -= L*K;
                        for (size_t l=L-1u; l>0u; --l)
                        {
                            for (size_t q=0u; q<=l; ++q, Y+=K, ++y) { *y = *Y; }
                            sc = *--y; Y -= K*(l+1u);
                            for (size_t q=0u; q<l; ++q, Y+=K) { --y; *Y += sc * *y; }
                            sc2 = fma(sc,-sc,1.0);
                            //*Y = -*Y;  //to match other sign convention
                            for (size_t q=0u; q<l; ++q) { Y-=K; *Y /= sc2; }
                        }
                        //*Y = -*Y;  //to match other sign convention
                    }
                }
            }
        }
        free(y);
    }

    return 0;
}


int ar2rc_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ar2rc_c: dim must be in [0 3]\n"); return 1; }
    
    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else
    {
        float scr, sci, sc2r, sc2i, sc2a, yr, yi, *y;
        if (!(y=(float *)malloc(2u*L*sizeof(float)))) { fprintf(stderr,"error in ar2rc_s: problem with malloc. "); perror("malloc"); return 1; }
        
        if (L==N)
        {
            for (size_t l=0u; l<2u*L; ++l, ++X, ++Y) { *Y = -*X; }
            Y -= 2u*L;

            for (size_t l=L-1u; l>0u; --l)
            {
                for (size_t q=0u; q<=l; ++q) { *y++ = *Y++; *y++ = *Y++; }
                sci = *--y; scr = *--y;
                Y -= 2u*l + 2u;
                for (size_t q=0u; q<l; ++q)
                {
                    y -= 2;
                    *Y++ += scr**y - sci**(y+1);
                    *Y++ += scr**(y+1) + sci**y;
                }
                sc2r = 1.0f - scr*scr + sci*sci;
                sc2i = -2.0f*scr*sci;
                sc2a = sc2r*sc2r + sc2i*sc2i;
                //*Y = -*Y; *(Y+1) = -*(Y+1);  //to match other sign convention
                for (size_t q=0u; q<l; ++q)
                {
                    yi = *--Y; yr = *--Y;
                    *Y = (yr*sc2r+yi*sc2i) / sc2a;
                    *(Y+1) = (yi*sc2r-yr*sc2i) / sc2a;
                }
            }
            //*Y = -*Y; *(Y+1) = -*(Y+1);  //to match other sign convention
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, Y+=2u*L)
                {
                    for (size_t l=0u; l<2u*L; ++l, ++X, ++Y) { *Y = -*X; }
                    Y -= 2u*L;

                    for (size_t l=L-1u; l>0u; --l)
                    {
                        for (size_t q=0u; q<=l; ++q) { *y++ = *Y++; *y++ = *Y++; }
                        sci = *--y; scr = *--y;
                        Y -= 2u*l + 2u;
                        for (size_t q=0u; q<l; ++q)
                        {
                            y -= 2;
                            *Y++ += scr**y - sci**(y+1);
                            *Y++ += scr**(y+1) + sci**y;
                        }
                        sc2r = 1.0f - scr*scr + sci*sci;
                        sc2i = -2.0f*scr*sci;
                        sc2a = sc2r*sc2r + sc2i*sc2i;
                        //*Y = -*Y; *(Y+1) = -*(Y+1);  //to match other sign convention
                        for (size_t q=0u; q<l; ++q)
                        {
                            yi = *--Y; yr = *--Y;
                            *Y = (yr*sc2r+yi*sc2i) / sc2a;
                            *(Y+1) = (yi*sc2r-yr*sc2i) / sc2a;
                        }
                    }
                    //*Y = -*Y; *(Y+1) = -*(Y+1);  //to match other sign convention
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=2u*K*L-2u, Y+=2u)
                    {
                        for (size_t l=0u; l<L; ++l, X+=2u*K, Y+=2u*K)
                        {
                            *Y = -*X; *(Y+1) = -*(X+1);
                        }
                        Y -= 2u*K*L;

                        for (size_t l=L-1u; l>0u; --l)
                        {
                            for (size_t q=0u; q<=l; ++q, Y+=2u*K) { *y++ = *Y; *y++ = *(Y+1); }
                            sci = *--y; scr = *--y;
                            Y -= 2u*K*(l+1u);
                            for (size_t q=0u; q<l; ++q, Y+=2u*K)
                            {
                                y -= 2;
                                *Y += scr**y - sci**(y+1);
                                *(Y+1) += scr**(y+1) + sci**y;
                            }
                            sc2r = 1.0f - scr*scr + sci*sci;
                            sc2i = -2.0f*scr*sci;
                            sc2a = sc2r*sc2r + sc2i*sc2i;
                            //*Y = -*Y; *(Y+1) = -*(Y+1);  //to match other sign convention
                            for (size_t q=0u; q<l; ++q)
                            {
                                Y -= 2u*K;
                                yr = *Y; yi = *(Y+1);
                                *Y = (yr*sc2r+yi*sc2i) / sc2a;
                                *(Y+1) = (yi*sc2r-yr*sc2i) / sc2a;
                            }
                        }
                        //*Y = -*Y; *(Y+1) = -*(Y+1);  //to match other sign convention
                    }
                }
            }
        }
        free(y);
    }

    return 0;
}


int ar2rc_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ar2rc_z: dim must be in [0 3]\n"); return 1; }
    
    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else
    {
        double scr, sci, sc2r, sc2i, sc2a, yr, yi, *y;
        if (!(y=(double *)malloc(2u*L*sizeof(double)))) { fprintf(stderr,"error in ar2rc_z: problem with malloc. "); perror("malloc"); return 1; }
        
        if (L==N)
        {
            for (size_t l=0u; l<2u*L; ++l, ++X, ++Y) { *Y = -*X; }
            Y -= 2u*L;

            for (size_t l=L-1u; l>0u; --l)
            {
                for (size_t q=0u; q<=l; ++q) { *y++ = *Y++; *y++ = *Y++; }
                sci = *--y; scr = *--y;
                Y -= 2u*l + 2u;
                for (size_t q=0u; q<l; ++q)
                {
                    y -= 2;
                    *Y++ += scr**y - sci**(y+1);
                    *Y++ += scr**(y+1) + sci**y;
                }
                sc2r = 1.0 - scr*scr + sci*sci;
                sc2i = -2.0*scr*sci;
                sc2a = sc2r*sc2r + sc2i*sc2i;
                //*Y = -*Y; *(Y+1) = -*(Y+1);  //to match other sign convention
                for (size_t q=0u; q<l; ++q)
                {
                    yi = *--Y; yr = *--Y;
                    *Y = (yr*sc2r+yi*sc2i) / sc2a;
                    *(Y+1) = (yi*sc2r-yr*sc2i) / sc2a;
                }
                //*Y = -*Y; *(Y+1) = -*(Y+1);  //to match other sign convention
            }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, Y+=2u*L)
                {
                    for (size_t l=0u; l<2u*L; ++l, ++X, ++Y) { *Y = -*X; }
                    Y -= 2u*L;

                    for (size_t l=L-1u; l>0u; --l)
                    {
                        for (size_t q=0u; q<=l; ++q) { *y++ = *Y++; *y++ = *Y++; }
                        sci = *--y; scr = *--y;
                        Y -= 2u*l + 2u;
                        for (size_t q=0u; q<l; ++q)
                        {
                            y -= 2;
                            *Y++ += scr**y - sci**(y+1);
                            *Y++ += scr**(y+1) + sci**y;
                        }
                        sc2r = 1.0 - scr*scr + sci*sci;
                        sc2i = -2.0*scr*sci;
                        sc2a = sc2r*sc2r + sc2i*sc2i;
                        //*Y = -*Y; *(Y+1) = -*(Y+1);  //to match other sign convention
                        for (size_t q=0u; q<l; ++q)
                        {
                            yi = *--Y; yr = *--Y;
                            *Y = (yr*sc2r+yi*sc2i) / sc2a;
                            *(Y+1) = (yi*sc2r-yr*sc2i) / sc2a;
                        }
                    }
                    //*Y = -*Y; *(Y+1) = -*(Y+1);  //to match other sign convention
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=2u*K*L-2u, Y+=2u)
                    {
                        for (size_t l=0u; l<L; ++l, X+=2u*K, Y+=2u*K)
                        {
                            *Y = -*X; *(Y+1) = -*(X+1);
                        }
                        Y -= 2u*K*L;

                        for (size_t l=L-1u; l>0u; --l)
                        {
                            for (size_t q=0u; q<=l; ++q, Y+=2u*K) { *y++ = *Y; *y++ = *(Y+1); }
                            sci = *--y; scr = *--y;
                            Y -= 2u*K*(l+1u);
                            for (size_t q=0u; q<l; ++q, Y+=2u*K)
                            {
                                y -= 2;
                                *Y += scr**y - sci**(y+1);
                                *(Y+1) += scr**(y+1) + sci**y;
                            }
                            sc2r = 1.0 - scr*scr + sci*sci;
                            sc2i = -2.0*scr*sci;
                            sc2a = sc2r*sc2r + sc2i*sc2i;
                            //*Y = -*Y; *(Y+1) = -*(Y+1);  //to match other sign convention
                            for (size_t q=0u; q<l; ++q)
                            {
                                Y -= 2u*K;
                                yr = *Y; yi = *(Y+1);
                                *Y = (yr*sc2r+yi*sc2i) / sc2a;
                                *(Y+1) = (yi*sc2r-yr*sc2i) / sc2a;
                            }
                        }
                        //*Y = -*Y; *(Y+1) = -*(Y+1);  //to match other sign convention
                    }
                }
            }
        }
        free(y);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
