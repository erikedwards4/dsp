//Gauss window with length L and stdev param a = (L-1)/(2*stdev)

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int gauss_s (float *X, const size_t L, const float r, const char normalize);
int gauss_d (double *X, const size_t L, const double r, const char normalize);
int gauss_c (float *X, const size_t L, const float r, const char normalize);
int gauss_z (double *X, const size_t L, const double r, const char normalize);


int gauss_s (float *X, const size_t L, const float a, const char normalize)
{
    if (a<=0.0f) { fprintf(stderr,"error in gauss_s: a must be positive \n"); return 1; }

    const float p = -2.0f*(a*a/((L-1)*(L-1)));
    const float m = 0.5f*(L-1);
    size_t l = 0;

    while (l<L/2) { X[l] = expf((l-m)*(l-m)*p); l++; }
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


int gauss_d (double *X, const size_t L, const double a, const char normalize)
{
    if (a<=0.0) { fprintf(stderr,"error in gauss_d: a must be positive \n"); return 1; }

    const double p = -2.0*(a*a/((L-1)*(L-1)));
    const double m = 0.5*(L-1);
    size_t l = 0;

    while (l<L/2) { X[l] = exp((l-m)*(l-m)*p); l++; }
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


int gauss_c (float *X, const size_t L, const float a, const char normalize)
{
    if (a<=0.0f) { fprintf(stderr,"error in gauss_c: a must be positive \n"); return 1; }

    const float p = -2.0f*(a*a/((L-1)*(L-1)));
    const float m = 0.5f*(L-1);
    const float z = 0.0f;
    size_t l = 0;

    cblas_scopy((int)L,&z,0,&X[1],2);

    while (l<L/2) { X[2*l] = expf((l-m)*(l-m)*p); l++; }
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


int gauss_z (double *X, const size_t L, const double a, const char normalize)
{
    if (a<=0.0) { fprintf(stderr,"error in gauss_z: a must be psoitive \n"); return 1; }
    
    const double p = -2.0*(a*a/((L-1)*(L-1)));
    const double m = 0.5*(L-1);
    const double z = 0.0;
    size_t l = 0;

    cblas_dcopy((int)L,&z,0,&X[1],2);

    while (l<L/2) { X[2*l] = exp((l-m)*(l-m)*p); l++; }
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
