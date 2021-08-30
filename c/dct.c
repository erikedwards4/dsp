//This computes "the" DCT (discrete cosine transformation) along dim of matrix X.
//This uses a matrix multiplication by the DCT-II matrix.

//For complex input X, the output Y is complex, and the DCT is just the DCT of the
//real and imaginary parts separately (following Octave convention).

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int dct_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int dct_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int dct_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int dct_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);


int dct_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in dct_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndct<Lx) { fprintf(stderr,"error in dct_s: ndct must be >= Lx (length of vecs in X)\n"); return 1; }

    if (ndct==0u || N==0u) {}
    else if (ndct==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        //Scaling
        const float s = sc ? 2.0f/sqrtf((float)(2u*ndct)) : 2.0f;
        const float dcsc = sc ? 1.0f/sqrtf((float)ndct) : 2.0f;

        //Initialize DCT-II matrix
        const size_t LN1 = Lx * (ndct-1u);
        const float P_N = (float)(M_PI/(double)ndct);
        float *DCT, sm;
        DCT = (float *)aligned_alloc(sizeof(float),LN1*sizeof(float));
        if (!DCT) { fprintf(stderr,"error in dct_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t n=1u; n<ndct; ++n)
        {
            for (size_t l=0u; l<Lx; ++l, ++DCT)
            {
                *DCT = s * cosf(P_N*(0.5f+(float)l)*(float)n);
            }
        }
        DCT -= LN1;
    
        if (Lx==N)
        {
            //DC term
            sm = 0.0f;
            for (size_t l=Lx; l>0u; --l, ++X) { sm += *X; }
            *Y = dcsc*sm; ++Y; X-=Lx;

            //Multiply by DCT matrix
            for (size_t n=ndct; n>1u; --n, X-=Lx, ++Y)
            {
                sm = 0.0f;
                for (size_t l=Lx; l>0u; --l, ++X, ++DCT) { sm += *X * *DCT; }
                *Y = sm;
            }
            DCT -= LN1;
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
                    //DC term
                    sm = 0.0f;
                    for (size_t l=Lx; l>0u; --l, ++X) { sm += *X; }
                    *Y = dcsc*sm; ++Y; X-=Lx;

                    //Multiply by DCT matrix
                    for (size_t n=ndct; n>2u; --n, X-=Lx, ++Y)
                    {
                        sm = 0.0f;
                        for (size_t l=Lx; l>0u; --l, ++X, ++DCT) { sm += *X * *DCT; }
                        *Y = sm;
                    }
                    sm = 0.0f;
                    for (size_t l=Lx; l>0u; --l, ++X, ++DCT) { sm += *X * *DCT; }
                    *Y++ = sm;
                    DCT -= LN1;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(ndct-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, Y-=K*ndct-1u)
                    {
                        //DC term
                        sm = 0.0f;
                        for (size_t l=Lx; l>0u; --l, X+=K) { sm += *X; }
                        *Y = dcsc*sm; Y+=K; X-=K*Lx;

                        //Multiply by DCT matrix
                        for (size_t n=ndct; n>1u; --n, X-=K*Lx, Y+=K)
                        {
                            sm = 0.0f;
                            for (size_t l=Lx; l>0u; --l, X+=K, ++DCT) { sm += *X * *DCT; }
                            *Y = sm;
                        }
                        DCT -= LN1;
                    }
                }
            }
        }
        free(DCT);
    }

    return 0;
}


int dct_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in dct_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndct<Lx) { fprintf(stderr,"error in dct_d: ndct must be >= Lx (length of vecs in X)\n"); return 1; }

    if (ndct==0u || N==0u) {}
    else if (ndct==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        //Scaling
        const double s = sc ? 2.0/sqrt((double)(2u*ndct)) : 2.0;
        const double dcsc = sc ? 1.0/sqrt((double)ndct) : 2.0;

        //Initialize DCT-II matrix
        const size_t LN1 = Lx * (ndct-1u);
        const double P_N = M_PI/(double)ndct;
        double *DCT, sm;
        DCT = (double *)aligned_alloc(sizeof(double),LN1*sizeof(double));
        if (!DCT) { fprintf(stderr,"error in dct_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t n=1u; n<ndct; ++n)
        {
            for (size_t l=0u; l<Lx; ++l, ++DCT)
            {
                *DCT = s * cos(P_N*(0.5+(double)l)*(double)n);
            }
        }
        DCT -= LN1;
    
        if (Lx==N)
        {
            //DC term
            sm = 0.0;
            for (size_t l=Lx; l>0u; --l, ++X) { sm += *X; }
            *Y = dcsc*sm; ++Y; X-=Lx;

            //Multiply by DCT matrix
            for (size_t n=ndct; n>1u; --n, X-=Lx, ++Y)
            {
                sm = 0.0;
                for (size_t l=Lx; l>0u; --l, ++X, ++DCT) { sm += *X * *DCT; }
                *Y = sm;
            }
            DCT -= LN1;
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
                    //DC term
                    sm = 0.0;
                    for (size_t l=Lx; l>0u; --l, ++X) { sm += *X; }
                    *Y = dcsc*sm; ++Y; X-=Lx;

                    //Multiply by DCT matrix
                    for (size_t n=ndct; n>2u; --n, X-=Lx, ++Y)
                    {
                        sm = 0.0;
                        for (size_t l=Lx; l>0u; --l, ++X, ++DCT) { sm += *X * *DCT; }
                        *Y = sm;
                    }
                    sm = 0.0;
                    for (size_t l=Lx; l>0u; --l, ++X, ++DCT) { sm += *X * *DCT; }
                    *Y++ = sm;
                    DCT -= LN1;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(ndct-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, Y-=K*ndct-1u)
                    {
                        //DC term
                        sm = 0.0;
                        for (size_t l=Lx; l>0u; --l, X+=K) { sm += *X; }
                        *Y = dcsc*sm; Y+=K; X-=K*Lx;

                        //Multiply by DCT matrix
                        for (size_t n=ndct; n>1u; --n, X-=K*Lx, Y+=K)
                        {
                            sm = 0.0;
                            for (size_t l=Lx; l>0u; --l, X+=K, ++DCT) { sm += *X * *DCT; }
                            *Y = sm;
                        }
                        DCT -= LN1;
                    }
                }
            }
        }
        free(DCT);
    }

    return 0;
}


int dct_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in dct_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndct<Lx) { fprintf(stderr,"error in dct_c: ndct must be >= Lx (length of vecs in X)\n"); return 1; }

    if (ndct==0u || N==0u) {}
    else if (ndct==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        //Scaling
        const float s = sc ? 2.0f/sqrtf((float)(2u*ndct)) : 2.0f;
        const float dcsc = sc ? 1.0f/sqrtf((float)ndct) : 2.0f;

        //Initialize DCT-II matrix
        const size_t LN1 = Lx * (ndct-1u);
        const float P_N = (float)(M_PI/(double)ndct);
        float *DCT, smr, smi;
        DCT = (float *)aligned_alloc(sizeof(float),LN1*sizeof(float));
        if (!DCT) { fprintf(stderr,"error in dct_c: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t n=1u; n<ndct; ++n)
        {
            for (size_t l=0u; l<Lx; ++l, ++DCT)
            {
                *DCT = s * cosf(P_N*(0.5f+(float)l)*(float)n);
            }
        }
        DCT -= LN1;
    
        if (Lx==N)
        {
            //DC term
            smr = smi = 0.0f;
            for (size_t l=Lx; l>0u; --l, ++X) { smr += *X; smi += *++X; }
            *Y++ = dcsc*smr; *Y++ = dcsc*smi; X -= 2u*Lx;

            //Multiply by DCT matrix
            for (size_t n=ndct; n>1u; --n, X-=2u*Lx, ++Y)
            {
                smr = smi = 0.0f;
                for (size_t l=Lx; l>0u; --l, ++X, ++DCT) { smr += *X * *DCT; smi += *++X * *DCT; }
                *Y = smr; *++Y = smi;
            }
            DCT -= LN1;
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
                    //DC term
                    smr = smi = 0.0f;
                    for (size_t l=Lx; l>0u; --l, ++X) { smr += *X; smi += *++X; }
                    *Y++ = dcsc*smr; *Y++ = dcsc*smi; X -= 2u*Lx;

                    //Multiply by DCT matrix
                    for (size_t n=ndct; n>2u; --n, X-=2u*Lx, ++Y)
                    {
                        smr = smi = 0.0f;
                        for (size_t l=Lx; l>0u; --l, ++X, ++DCT) { smr += *X * *DCT; smi += *++X * *DCT; }
                        *Y = smr; *++Y = smi;
                    }
                    smr = smi = 0.0f;
                    for (size_t l=Lx; l>0u; --l, ++X, ++DCT) { smr += *X * *DCT; smi += *++X * *DCT; }
                    *Y++ = smr; *Y++ = smi;
                    DCT -= LN1;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(ndct-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2, Y-=2u*K*ndct-2u)
                    {
                        //DC term
                        smr = smi = 0.0f;
                        for (size_t l=Lx; l>0u; --l, X+=2u*K) { smr += *X; smi += *(X+1); }
                        *Y = dcsc*smr; *(Y+1) = dcsc*smi; Y+=2u*K; X-=2u*K*Lx;

                        //Multiply by DCT matrix
                        for (size_t n=ndct; n>1u; --n, X-=2u*K*Lx, Y+=2u*K)
                        {
                            smr = smi = 0.0f;
                            for (size_t l=Lx; l>0u; --l, X+=2u*K, ++DCT) { smr += *X**DCT; smi += *(X+1)**DCT; }
                            *Y = smr; *(Y+1) = smi;
                        }
                        DCT -= LN1;
                    }
                }
            }
        }
        free(DCT);
    }

    return 0;
}


int dct_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in dct_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndct<Lx) { fprintf(stderr,"error in dct_z: ndct must be >= Lx (length of vecs in X)\n"); return 1; }

    if (ndct==0u || N==0u) {}
    else if (ndct==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        //Scaling
        const double s = sc ? 2.0/sqrt((double)(2u*ndct)) : 2.0;
        const double dcsc = sc ? 1.0/sqrt((double)ndct) : 2.0;

        //Initialize DCT-II matrix
        const size_t LN1 = Lx * (ndct-1u);
        const double P_N = M_PI/(double)ndct;
        double *DCT, smr, smi;
        DCT = (double *)aligned_alloc(sizeof(double),LN1*sizeof(double));
        if (!DCT) { fprintf(stderr,"error in dct_z: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t n=1u; n<ndct; ++n)
        {
            for (size_t l=0u; l<Lx; ++l, ++DCT)
            {
                *DCT = s * cos(P_N*(0.5+(double)l)*(double)n);
            }
        }
        DCT -= LN1;
    
        if (Lx==N)
        {
            //DC term
            smr = smi = 0.0;
            for (size_t l=Lx; l>0u; --l, ++X) { smr += *X; smi += *++X; }
            *Y++ = dcsc*smr; *Y++ = dcsc*smi; X -= 2u*Lx;

            //Multiply by DCT matrix
            for (size_t n=ndct; n>1u; --n, X-=2u*Lx, ++Y)
            {
                smr = smi = 0.0;
                for (size_t l=Lx; l>0u; --l, ++X, ++DCT) { smr += *X * *DCT; smi += *++X * *DCT; }
                *Y = smr; *++Y = smi;
            }
            DCT -= LN1;
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
                    //DC term
                    smr = smi = 0.0;
                    for (size_t l=Lx; l>0u; --l, ++X) { smr += *X; smi += *++X; }
                    *Y++ = dcsc*smr; *Y++ = dcsc*smi; X -= 2u*Lx;

                    //Multiply by DCT matrix
                    for (size_t n=ndct; n>2u; --n, X-=2u*Lx, ++Y)
                    {
                        smr = smi = 0.0;
                        for (size_t l=Lx; l>0u; --l, ++X, ++DCT) { smr += *X * *DCT; smi += *++X * *DCT; }
                        *Y = smr; *++Y = smi;
                    }
                    smr = smi = 0.0;
                    for (size_t l=Lx; l>0u; --l, ++X, ++DCT) { smr += *X * *DCT; smi += *++X * *DCT; }
                    *Y++ = smr; *Y++ = smi;
                    DCT -= LN1;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(ndct-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2, Y-=2u*K*ndct-2u)
                    {
                        //DC term
                        smr = smi = 0.0;
                        for (size_t l=Lx; l>0u; --l, X+=2u*K) { smr += *X; smi += *(X+1); }
                        *Y = dcsc*smr; *(Y+1) = dcsc*smi; Y+=2u*K; X-=2u*K*Lx;

                        //Multiply by DCT matrix
                        for (size_t n=ndct; n>1u; --n, X-=2u*K*Lx, Y+=2u*K)
                        {
                            smr = smi = 0.0;
                            for (size_t l=Lx; l>0u; --l, X+=2u*K, ++DCT) { smr += *X**DCT; smi += *(X+1)**DCT; }
                            *Y = smr; *(Y+1) = smi;
                        }
                        DCT -= LN1;
                    }
                }
            }
        }
        free(DCT);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
