//This computes "the" IDST (discrete sine transformation) along dim of matrix X.
//This uses direct matrix multiplication by the DST-I matrix.
//This is the same as DST-I except scaling.

//For complex input X, the output Y is complex, and the IDST is just the IDST of the
//real and imaginary parts separately (following Octave convention).

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "codee_dsp.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int idst_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idst_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndst<Lx) { fprintf(stderr,"error in idst_s: ndst must be >= Lx (length of vecs in X)\n"); return 1; }

    if (ndst==0u || N==0u) {}
    else if (ndst==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        //Scaling
        const float s = sc ? 2.0f/(float)(ndst+1u) : 2.0f/(float)(2u*ndst+2u);

        //Initialize DST-I matrix
        const size_t LN = Lx * ndst;
        const float P_N = (float)(M_PI/(double)(ndst+1u));
        float *DST, sm;
        DST = (float *)aligned_alloc(sizeof(float),LN*sizeof(float));
        if (!DST) { fprintf(stderr,"error in idst_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t n=0u; n<ndst; ++n)
        {
            for (size_t l=0u; l<Lx; ++l, ++DST)
            {
                *DST = s * sinf(P_N*(float)(l+1u)*(float)(n+1u));
            }
        }
        DST -= LN;
    
        if (Lx==N)
        {
            for (size_t n=ndst; n>0u; --n, X-=Lx, ++Y)
            {
                sm = 0.0f;
                for (size_t l=Lx; l>0u; --l, ++X, ++DST) { sm += *X * *DST; }
                *Y = sm;
            }
            DST -= LN;
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
                    for (size_t n=ndst; n>1u; --n, X-=Lx, ++Y)
                    {
                        sm = 0.0f;
                        for (size_t l=Lx; l>0u; --l, ++X, ++DST) { sm += *X * *DST; }
                        *Y = sm;
                    }
                    sm = 0.0f;
                    for (size_t l=Lx; l>0u; --l, ++X, ++DST) { sm += *X * *DST; }
                    *Y++ = sm;
                    DST -= LN;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(ndst-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, Y-=K*ndst-1u)
                    {
                        for (size_t n=ndst; n>0u; --n, X-=K*Lx, Y+=K)
                        {
                            sm = 0.0f;
                            for (size_t l=Lx; l>0u; --l, X+=K, ++DST) { sm += *X * *DST; }
                            *Y = sm;
                        }
                        DST -= LN;
                    }
                }
            }
        }
        free(DST);
    }

    return 0;
}


int idst_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idst_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndst<Lx) { fprintf(stderr,"error in idst_d: ndst must be >= Lx (length of vecs in X)\n"); return 1; }

    if (ndst==0u || N==0u) {}
    else if (ndst==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        //Scaling
        const double s = sc ? 2.0/(double)(ndst+1u) : 2.0/(double)(2u*ndst+2u);

        //Initialize DST-I matrix
        const size_t LN = Lx * ndst;
        const double P_N = M_PI/(double)(ndst+1u);
        double *DST, sm;
        DST = (double *)aligned_alloc(sizeof(double),LN*sizeof(double));
        if (!DST) { fprintf(stderr,"error in idst_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t n=0u; n<ndst; ++n)
        {
            for (size_t l=0u; l<Lx; ++l, ++DST)
            {
                *DST = s * sin(P_N*(double)(l+1u)*(double)(n+1u));
            }
        }
        DST -= LN;
    
        if (Lx==N)
        {
            for (size_t n=ndst; n>0u; --n, X-=Lx, ++Y)
            {
                sm = 0.0;
                for (size_t l=Lx; l>0u; --l, ++X, ++DST) { sm += *X * *DST; }
                *Y = sm;
            }
            DST -= LN;
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
                    for (size_t n=ndst; n>1u; --n, X-=Lx, ++Y)
                    {
                        sm = 0.0;
                        for (size_t l=Lx; l>0u; --l, ++X, ++DST) { sm += *X * *DST; }
                        *Y = sm;
                    }
                    sm = 0.0;
                    for (size_t l=Lx; l>0u; --l, ++X, ++DST) { sm += *X * *DST; }
                    *Y++ = sm;
                    DST -= LN;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(ndst-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, Y-=K*ndst-1u)
                    {
                        for (size_t n=ndst; n>0u; --n, X-=K*Lx, Y+=K)
                        {
                            sm = 0.0;
                            for (size_t l=Lx; l>0u; --l, X+=K, ++DST) { sm += *X * *DST; }
                            *Y = sm;
                        }
                        DST -= LN;
                    }
                }
            }
        }
        free(DST);
    }

    return 0;
}


int idst_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idst_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndst<Lx) { fprintf(stderr,"error in idst_c: ndst must be >= Lx (length of vecs in X)\n"); return 1; }

    if (ndst==0u || N==0u) {}
    else if (ndst==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        //Scaling
        const float s = sc ? 2.0f/(float)(ndst+1u) : 2.0f/(float)(2u*ndst+2u);

        //Initialize DST-I matrix
        const size_t LN = Lx * ndst;
        const float P_N = (float)(M_PI/(double)(ndst+1u));
        float *DST, smr, smi;
        DST = (float *)aligned_alloc(sizeof(float),LN*sizeof(float));
        if (!DST) { fprintf(stderr,"error in idst_c: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t n=0u; n<ndst; ++n)
        {
            for (size_t l=0u; l<Lx; ++l, ++DST)
            {
                *DST = s * sinf(P_N*(float)(l+1u)*(float)(n+1u));
            }
        }
        DST -= LN;
    
        if (Lx==N)
        {
            for (size_t n=ndst; n>0u; --n, X-=2u*Lx, ++Y)
            {
                smr = smi = 0.0f;
                for (size_t l=Lx; l>0u; --l, ++X, ++DST) { smr += *X * *DST; smi += *++X * *DST; }
                *Y = smr; *++Y = smi;
            }
            DST -= LN;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, ++Y)
                {
                    for (size_t n=ndst; n>1u; --n, X-=2u*Lx, ++Y)
                    {
                        smr = smi = 0.0f;
                        for (size_t l=Lx; l>0u; --l, ++X, ++DST) { smr += *X * *DST; smi += *++X * *DST; }
                        *Y = smr; *++Y = smi;
                    }
                    smr = smi = 0.0f;
                    for (size_t l=Lx; l>0u; --l, ++X, ++DST) { smr += *X * *DST; smi += *++X * *DST; }
                    *Y = smr; *++Y = smi;
                    DST -= LN;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(ndst-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u, Y-=2u*K*ndst-2u)
                    {
                        for (size_t n=ndst; n>0u; --n, X-=2u*K*Lx, Y+=2u*K)
                        {
                            smr = smi = 0.0f;
                            for (size_t l=Lx; l>0u; --l, X+=2u*K, ++DST) { smr += *X**DST; smi += *(X+1)**DST; }
                            *Y = smr; *(Y+1) = smi;
                        }
                        DST -= LN;
                    }
                }
            }
        }
        free(DST);
    }

    return 0;
}


int idst_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idst_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndst<Lx) { fprintf(stderr,"error in idst_z: ndst must be >= Lx (length of vecs in X)\n"); return 1; }

    if (ndst==0u || N==0u) {}
    else if (ndst==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        //Scaling
        const double s = sc ? 2.0/(double)(ndst+1u) : 2.0/(double)(2u*ndst+2u);

        //Initialize DST-I matrix
        const size_t LN = Lx * ndst;
        const double P_N = M_PI/(double)(ndst+1u);
        double *DST, smr, smi;
        DST = (double *)aligned_alloc(sizeof(double),LN*sizeof(double));
        if (!DST) { fprintf(stderr,"error in idst_z: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t n=0u; n<ndst; ++n)
        {
            for (size_t l=0u; l<Lx; ++l, ++DST)
            {
                *DST = s * sin(P_N*(double)(l+1u)*(double)(n+1u));
            }
        }
        DST -= LN;
    
        if (Lx==N)
        {
            for (size_t n=ndst; n>0u; --n, X-=2u*Lx, ++Y)
            {
                smr = smi = 0.0;
                for (size_t l=Lx; l>0u; --l, ++X, ++DST) { smr += *X * *DST; smi += *++X * *DST; }
                *Y = smr; *++Y = smi;
            }
            DST -= LN;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, ++Y)
                {
                    for (size_t n=ndst; n>1u; --n, X-=2u*Lx, ++Y)
                    {
                        smr = smi = 0.0;
                        for (size_t l=Lx; l>0u; --l, ++X, ++DST) { smr += *X * *DST; smi += *++X * *DST; }
                        *Y = smr; *++Y = smi;
                    }
                    smr = smi = 0.0;
                    for (size_t l=Lx; l>0u; --l, ++X, ++DST) { smr += *X * *DST; smi += *++X * *DST; }
                    *Y = smr; *++Y = smi;
                    DST -= LN;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(ndst-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u, Y-=2u*K*ndst-2u)
                    {
                        for (size_t n=ndst; n>0u; --n, X-=2u*K*Lx, Y+=2u*K)
                        {
                            smr = smi = 0.0;
                            for (size_t l=Lx; l>0u; --l, X+=2u*K, ++DST) { smr += *X**DST; smi += *(X+1)**DST; }
                            *Y = smr; *(Y+1) = smi;
                        }
                        DST -= LN;
                    }
                }
            }
        }
        free(DST);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
