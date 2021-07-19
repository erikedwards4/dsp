//There is an 'exact' blackman with:
//a0 = 7938/18608 ≈ 0.42659, a1 = 9240/18608 ≈ 0.49656, and a2 = 1430/18608 ≈ 0.076849
//However:
//a0 = 0.42, a1 = 0.5, a2 = 0.08
//are more-often used (they have improved 18 db/oct fall-off but do not null the sidelobes as well [Wikepedia]).
//However, this later window has negative values at the edges, so I include the exact version here as an option.

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

int blackman_s (float *X, const size_t L, const char normalize, const char exact);
int blackman_d (double *X, const size_t L, const char normalize, const char exact);
int blackman_c (float *X, const size_t L, const char normalize, const char exact);
int blackman_z (double *X, const size_t L, const char normalize, const char exact);


int blackman_s (float *X, const size_t L, const char normalize, const char exact)
{
    const float p = 2.0f*M_PIf/(L-1.0f);
    size_t l = 0;

    if (exact)
    {
        while (l<L/2) { X[l] = 7938.0f/18608.0f - (9240.0f/18608.0f)*cosf(p*l) + (1430.0f/18608.0f)*cosf(2.0f*p*l); l++; }
    }
    else
    {
        while (l<L/2) { X[l] = 0.42f - 0.5f*cosf(p*l) + 0.08f*cosf(2.0f*p*l); l++; }
    }
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


int blackman_d (double *X, const size_t L, const char normalize, const char exact)
{
    const double p = 2.0*M_PI/(L-1.0);
    size_t l = 0;

    if (exact)
    {
        while (l<L/2) { X[l] = 7938.0/18608.0 - (9240.0/18608.0)*cos(p*l) + (1430.0/18608.0)*cos(2.0*p*l); l++; }
    }
    else
    {
        while (l<L/2) { X[l] = 0.42 - 0.5*cos(p*l) + 0.08*cos(2.0*p*l); l++; }
    }
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


int blackman_c (float *X, const size_t L, const char normalize, const char exact)
{
    const float z = 0.0f, p = 2.0f*M_PIf/(L-1.0f);
    size_t l = 0;

    cblas_scopy((int)L,&z,0,&X[1],2);

    if (exact)
    {
        while (l<L/2) { X[2*l] = 7938.0f/18608.0f - (9240.0f/18608.0f)*cosf(p*l) + (1430.0f/18608.0f)*cosf(2.0f*p*l); l++; }
    }
    else
    {
        while (l<L/2) { X[2*l] = 0.42f - 0.5f*cosf(p*l) + 0.08f*cosf(2.0f*p*l); l++; }
    }
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


int blackman_z (double *X, const size_t L, const char normalize, const char exact)
{
    const double z = 0.0, p = 2.0*M_PI/(L-1.0);
    size_t l = 0;

    cblas_dcopy((int)L,&z,0,&X[1],2);

    if (exact)
    {
        while (l<L/2) { X[2*l] = 7938.0/18608.0 - (9240.0/18608.0)*cos(p*l) + (1430.0/18608.0)*cos(2.0*p*l); l++; }
    }
    else
    {
        while (l<L/2) { X[2*l] = 0.42 - 0.5*cos(p*l) + 0.08*cos(2.0*p*l); l++; }
    }
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
