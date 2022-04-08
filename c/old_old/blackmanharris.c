//Note that this is a continuation of a series of generalized symmetric cosine windows,
//with [Wikipedia]: a0 = 0.35875, a1 = 0.48829, a2 = 0.14128, a3 = 0.01168;

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifndef M_PIf
   #define M_PIf 3.141592653589793238462643383279502884f
#endif

#ifndef M_PI
   #define M_PI 3.141592653589793238462643383279502884
#endif

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int blackmanharris_s (float *X, const size_t L, const char normalize);
int blackmanharris_d (double *X, const size_t L, const char normalize);
int blackmanharris_c (float *X, const size_t L, const char normalize);
int blackmanharris_z (double *X, const size_t L, const char normalize);


int blackmanharris_s (float *X, const size_t L, const char normalize)
{
    const float p = 2.0f*M_PIf/(L-1.0f);
    size_t l = 0;

    while (l<L/2) { X[l] = 0.35875f - 0.48829f*cosf(p*l) + 0.14128f*cosf(2.0f*p*l) - 0.01168f*cosf(3.0f*p*l); l++; }
    if (L%2) { X[l] = 1.0f; l++; }
    while (l<L) { X[l] = X[L-l-1]; l++; }

    if (normalize)
    {
        const float o = 1.0f;
        float sm = cblas_sdot((int)L,X,1,&o,0);
        cblas_sscal((int)L,1.0f/sm,X,1);
        sm = cblas_sdot((int)L,X,1,&o,0);
        X[L/2] += 1.0f - sm;
    }

    return 0;
}


int blackmanharris_d (double *X, const size_t L, const char normalize)
{
    const double p = 2.0*M_PI/(L-1.0);
    size_t l = 0;

    while (l<L/2) { X[l] = 0.35875 - 0.48829*cos(p*l) + 0.14128*cos(2.0*p*l) - 0.01168*cos(3.0*p*l); l++; }
    if (L%2) { X[l] = 1.0; l++; }
    while (l<L) { X[l] = X[L-l-1]; l++; }

    if (normalize)
    {
        const double o = 1.0;
        double sm = cblas_ddot((int)L,X,1,&o,0);
        cblas_dscal((int)L,1.0/sm,X,1);
        sm = cblas_ddot((int)L,X,1,&o,0);
        X[L/2] += 1.0 - sm;
    }

    return 0;
}


int blackmanharris_c (float *X, const size_t L, const char normalize)
{
    const float z = 0.0f, p = 2.0f*M_PIf/(L-1.0f);
    size_t l = 0;

    cblas_scopy((int)L,&z,0,&X[1],2);

    while (l<L/2) { X[2*l] = 0.35875f - 0.48829f*cosf(p*l) + 0.14128f*cosf(2.0f*p*l) - 0.01168f*cosf(3.0f*p*l); l++; }
    if (L%2) { X[2*l] = 1.0f; l++; }
    while (l<L) { X[2*l] = X[2*(L-l-1)]; l++; }

    if (normalize)
    {
        const float o = 1.0f;
        float sm = cblas_sdot((int)L,X,2,&o,0);
        cblas_sscal((int)L,1.0f/sm,X,2);
        sm = cblas_sdot((int)L,X,2,&o,0);
        X[2*(L/2)] += 1.0f - sm;
    }

    return 0;
}


int blackmanharris_z (double *X, const size_t L, const char normalize)
{
    const double z = 0.0, p = 2.0*M_PI/(L-1.0);
    size_t l = 0;

    cblas_dcopy((int)L,&z,0,&X[1],2);

    while (l<L/2) { X[2*l] = 0.35875 - 0.48829*cos(p*l) + 0.14128*cos(2.0*p*l) - 0.01168*cos(3.0*p*l); l++; }
    if (L%2) { X[2*l] = 1.0; l++; }
    while (l<L) { X[2*l] = X[2*(L-l-1)]; l++; }

    if (normalize)
    {
        const double o = 1.0;
        double sm = cblas_ddot((int)L,X,2,&o,0);
        cblas_dscal((int)L,1.0/sm,X,2);
        sm = cblas_ddot((int)L,X,2,&o,0);
        X[2*(L/2)] += 1.0 - sm;
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
