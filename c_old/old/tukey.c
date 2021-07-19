//Tukey (tapered cosine) window with length L and cosine fraction r

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

int tukey_s (float *X, const size_t L, const float r, const char normalize);
int tukey_d (double *X, const size_t L, const double r, const char normalize);
int tukey_c (float *X, const size_t L, const float r, const char normalize);
int tukey_z (double *X, const size_t L, const double r, const char normalize);


int tukey_s (float *X, const size_t L, const float r, const char normalize)
{
    if (r<0.0f || r>1.0f) { fprintf(stderr,"error in tukey_s: r must be in [0.0 1.0] \n"); return 1; }
    if (L==0) { return 0; }
    if (L==1) { *X = 1.0f; return 0; }

    const float p = 2.0f*M_PIf/((L-1)*r);
    size_t l = 0;
    
    while (l<0.5f*r*L) { X[l] = 0.5f - 0.5f*cosf(p*l); l++; }
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


int tukey_d (double *X, const size_t L, const double r, const char normalize)
{
    if (r<0.0 || r>1.0) { fprintf(stderr,"error in tukey_d: r must be in [0.0 1.0] \n"); return 1; }
    if (L==0) { return 0; }
    if (L==1) { *X = 1.0; return 0; }

    const double p = 2.0*M_PI/((L-1)*r);
    size_t l = 0;

    while (l<0.5*r*L) { X[l] = 0.5 - 0.5*cos(p*l); l++; }
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


int tukey_c (float *X, const size_t L, const float r, const char normalize)
{
    if (r<0.0f || r>1.0f) { fprintf(stderr,"error in tukey_c: r must be in [0.0 1.0] \n"); return 1; }
    if (L==0) { return 0; }
    if (L==1) { *X++ = 1.0f; *X = 0.0f; return 0; }

    const float z = 0.0f, p = 2.0f*M_PIf/((L-1)*r);
    size_t l = 0;

    cblas_scopy((int)L,&z,0,&X[1],2);

    while (l<0.5f*r*L) { X[2*l] = 0.5f - 0.5f*cosf(p*l); l++; }
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


int tukey_z (double *X, const size_t L, const double r, const char normalize)
{
    if (r<0.0 || r>1.0) { fprintf(stderr,"error in tukey_z: r must be in [0.0 1.0] \n"); return 1; }
    if (L==0) { return 0; }
    if (L==1) { *X++ = 1.0; *X = 0.0; return 0; }

    const double z = 0.0, p = 2.0*M_PI/((L-1)*r);
    size_t l = 0;

    cblas_dcopy((int)L,&z,0,&X[1],2);

    while (l<0.5*r*L) { X[2*l] = 0.5 - 0.5*cos(p*l); l++; }
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
