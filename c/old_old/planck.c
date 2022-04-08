//Tukey (tapered cosine) window with length L and cosine fraction r

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int planck_s (float *X, const size_t L, const float epsilon, const char normalize);
int planck_d (double *X, const size_t L, const double epsilon, const char normalize);
int planck_c (float *X, const size_t L, const float epsilon, const char normalize);
int planck_z (double *X, const size_t L, const double epsilon, const char normalize);


int planck_s (float *X, const size_t L, const float epsilon, const char normalize)
{
    if (epsilon<0.0f || epsilon>0.5f) { fprintf(stderr,"error in planck_s: epsilon must be in [0.0 0.5] \n"); return 1; }
    if (L==0) { return 0; }
    if (L==1) { *X = 1.0f; return 0; }

    const float p = epsilon*L;
    size_t l = 0;

    X[l++] = 0.0f;
    while (l<p && l<L/2) { X[l] = 1.0f/(1.0f+expf(p/l-p/(p-l))); l++; }
    while (l<L/2) { X[l] = 1.0f; l++; }
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


int planck_d (double *X, const size_t L, const double epsilon, const char normalize)
{
    if (epsilon<0.0 || epsilon>0.5) { fprintf(stderr,"error in planck_d: epsilon must be in [0.0 0.5] \n"); return 1; }
    if (L==0) { return 0; }
    if (L==1) { *X = 1.0; return 0; }

    const double p = epsilon*L;
    size_t l = 0;

    X[l++] = 0.0;
    while (l<p && l<L/2) { X[l] = 1.0/(1.0+exp(p/l-p/(p-l))); l++; }
    while (l<L/2) { X[l] = 1.0; l++; }
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


int planck_c (float *X, const size_t L, const float epsilon, const char normalize)
{
    if (epsilon<0.0f || epsilon>0.5f) { fprintf(stderr,"error in planck_c: epsilon must be in [0.0 0.5] \n"); return 1; }
    if (L==0) { return 0; }
    if (L==1) { *X++ = 1.0f; *X = 0.0f; return 0; }

    const float z = 0.0f, p = epsilon*L;
    size_t l = 0;

    cblas_scopy((int)L,&z,0,&X[1],2);

    X[l++] = 0.0f;
    while (l<p && l<L/2) { X[2*l] = 1.0f/(1.0f+expf(p/l-p/(p-l))); l++; }
    while (l<L/2) { X[2*l] = 1.0f; l++; }
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


int planck_z (double *X, const size_t L, const double epsilon, const char normalize)
{
    if (epsilon<0.0 || epsilon>0.5) { fprintf(stderr,"error in planck_z: epsilon must be in [0.0 0.5] \n"); return 1; }
    if (L==0) { return 0; }
    if (L==1) { *X++ = 1.0; *X = 0.0; return 0; }
    
    const double z = 0.0, p = epsilon*L;
    size_t l = 0;

    cblas_dcopy((int)L,&z,0,&X[1],2);

    X[l++] = 0.0;
    while (l<p && l<L/2) { X[2*l] = 1.0/(1.0+exp(p/l-p/(p-l))); l++; }
    while (l<L/2) { X[2*l] = 1.0; l++; }
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
