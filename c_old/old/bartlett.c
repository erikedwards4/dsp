//For Barlett, I use Octave/Matlab convention,
//for which they cite Openheim & Shafer: Discrete-Time Signal Processing.

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int bartlett_s (float *X, const size_t L, const char normalize);
int bartlett_d (double *X, const size_t L, const char normalize);
int bartlett_c (float *X, const size_t L, const char normalize);
int bartlett_z (double *X, const size_t L, const char normalize);


int bartlett_s (float *X, const size_t L, const char normalize)
{
    const float p = 1.0f/(L-(L%2)*(L/2)-1);
    size_t l = 0;

    while (l<L/2) { X[l] = (l*(2-L%2))*p; l++; }
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

    if (normalize)
    {
        float nrm;
        if (norm==1) { nrm = cblas_sasum((int)L,X,1); }
        else if (norm==2) { nrm = cblas_snrm2((int)L,X,1); }
        else { fprintf(stderr,"error in _s: norm must be in {0,1,2}\n") }
        cblas_sscal((int)L,1.0f/nrm,X,1);
    }

    return 0;
}


int bartlett_d (double *X, const size_t L, const char normalize)
{
    const double p = 1.0/(L-(L%2)*(L/2)-1);
    size_t l = 0;

    while (l<L/2) { X[l] = (l*(2-L%2))*p; l++; }
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


int bartlett_c (float *X, const size_t L, const char normalize)
{
    const float z = 0.0f, p = 1.0f/(L-(L%2)*(L/2)-1);
    size_t l = 0;

    cblas_scopy((int)L,&z,0,&X[1],2);

    while (l<L/2) { X[2*l] = (l*(2-L%2))*p; l++; }
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


int bartlett_z (double *X, const size_t L, const char normalize)
{
    const double z = 0.0, p = 1.0/(L-(L%2)*(L/2)-1);
    size_t l = 0;

    cblas_dcopy((int)L,&z,0,&X[1],2);

    while (l<L/2) { X[2*l] = (l*(2-L%2))*p; l++; }
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
