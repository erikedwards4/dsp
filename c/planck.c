//Planck window with length L and parameter epsilon.

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int planck_s (float *Y, const size_t L, const float epsilon, const size_t norm);
int planck_d (double *Y, const size_t L, const double epsilon, const size_t norm);
int planck_c (float *Y, const size_t L, const float epsilon, const size_t norm);
int planck_z (double *Y, const size_t L, const double epsilon, const size_t norm);


int planck_s (float *Y, const size_t L, const float epsilon, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in planck_s: norm must be in {0,1,2,3}\n"); return 1; }
    if (epsilon<0.0f || epsilon>0.5f) { fprintf(stderr,"error in planck_s: epsilon must be in [0.0 0.5] \n"); return 1; }

    if (L==0) {}
    else if (L==1) { *Y++ = 1.0f; }
    else if (L==2) { *Y++ = 1.0f; *Y++ = 1.0f; }
    else
    {
        const float p = epsilon*L;
        *Y++ = 0.0f;
        for (size_t l=1; l<p && l<L/2; ++l, ++Y) { *Y = 1.0f/(1.0f+expf(p/l-p/(p-l))); }
        for (size_t l=(size_t)ceilf(p); l<L/2+L%2; ++l, ++Y) { *Y = 1.0f; }
        for (size_t l=0; l<L/2; ++l, ++Y) { *Y = *(Y-2*l-1-L%2); }
    }

    if (norm)
    {
        float nrm = 0.0f;
        if (norm==1) { for (size_t l=0; l<L; ++l) { --Y; nrm += *Y; } }
        else if (norm==2) { for (size_t l=0; l<L; ++l) { --Y; nrm += *Y**Y; } nrm = sqrtf(nrm); }
        else if (norm==3) { Y -= L-L/2; nrm = *Y; Y -= L/2; }
        for (size_t l=0; l<L; ++l, ++Y) { *Y /= nrm; }
    }

    return 0;
}


int planck_d (double *Y, const size_t L, const double epsilon, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in planck_s: norm must be in {0,1,2,3}\n"); return 1; }
    if (epsilon<0.0 || epsilon>0.5) { fprintf(stderr,"error in planck_d: epsilon must be in [0.0 0.5] \n"); return 1; }

    if (L==0) {}
    else if (L==1) { *Y++ = 1.0; }
    else if (L==2) { *Y++ = 1.0; *Y++ = 1.0; }
    else
    {
        const double p = epsilon*L;
        *Y++ = 0.0;
        for (size_t l=1; l<p && l<L/2; ++l, ++Y) { *Y = 1.0/(1.0+exp(p/l-p/(p-l))); }
        for (size_t l=(size_t)ceil(p); l<L/2+L%2; ++l, ++Y) { *Y = 1.0; }
        for (size_t l=0; l<L/2; ++l, ++Y) { *Y = *(Y-2*l-1-L%2); }
    }

    if (norm)
    {
        double nrm = 0.0;
        if (norm==1) { for (size_t l=0; l<L; ++l) { --Y; nrm += *Y; } }
        else if (norm==2) { for (size_t l=0; l<L; ++l) { --Y; nrm += *Y**Y; } nrm = sqrt(nrm); }
        else if (norm==3) { Y -= L-L/2; nrm = *Y; Y -= L/2; }
        for (size_t l=0; l<L; ++l, ++Y) { *Y /= nrm; }
    }

    return 0;
}


int planck_c (float *Y, const size_t L, const float epsilon, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in planck_c: norm must be in {0,1,2,3}\n"); return 1; }
    if (epsilon<0.0f || epsilon>0.5f) { fprintf(stderr,"error in planck_c: epsilon must be in [0.0 0.5] \n"); return 1; }

    if (L==0) {}
    else if (L==1) { *Y++ = 1.0f; *Y++ = 0.0f; }
    else if (L==2) { *Y++ = 1.0f; *Y++ = 0.0f; *Y++ = 1.0f; *Y++ = 0.0f; }
    else
    {
        const float p = epsilon*L;
        *Y++ = 0.0f; *Y++ = 0.0f;
        for (size_t l=1; l<p && l<L/2; ++l, ++Y) { *Y = 1.0f/(1.0f+expf(p/l-p/(p-l))); *++Y = 0.0f; }
        for (size_t l=(size_t)ceilf(p); l<L/2+L%2; ++l, ++Y) { *Y = 1.0f; *++Y = 0.0f; }
        for (size_t l=0; l<L/2; ++l, ++Y) { *Y = *(Y-4*l-2-2*(L%2)); *++Y = 0.0f; }
    }

    if (norm)
    {
        float nrm = 0.0f;
        if (norm==1) { for (size_t l=0; l<L; ++l) { Y-=2; nrm += *Y; } }
        else if (norm==2) { for (size_t l=0; l<L; ++l) { Y-=2; nrm += *Y**Y; } nrm = sqrtf(nrm); }
        else if (norm==3) { Y -= 2*(L-L/2); nrm = *Y; Y -= 2*(L/2); }
        for (size_t l=0; l<L; ++l, Y+=2) { *Y /= nrm; }
    }

    return 0;
}


int planck_z (double *Y, const size_t L, const double epsilon, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in planck_z: norm must be in {0,1,2,3}\n"); return 1; }
    if (epsilon<0.0 || epsilon>0.5) { fprintf(stderr,"error in planck_z: epsilon must be in [0.0 0.5] \n"); return 1; }

    if (L==0) {}
    else if (L==1) { *Y++ = 1.0; *Y++ = 0.0; }
    else if (L==2) { *Y++ = 1.0; *Y++ = 0.0; *Y++ = 1.0; *Y++ = 0.0; }
    else
    {
        const double p = epsilon*L;
        *Y++ = 0.0; *Y++ = 0.0;
        for (size_t l=1; l<p && l<L/2; ++l, ++Y) { *Y = 1.0/(1.0+exp(p/l-p/(p-l))); *++Y = 0.0; }
        for (size_t l=(size_t)ceil(p); l<L/2+L%2; ++l, ++Y) { *Y = 1.0; *++Y = 0.0; }
        for (size_t l=0; l<L/2; ++l, ++Y) { *Y = *(Y-4*l-2-2*(L%2)); *++Y = 0.0; }
    }

    if (norm)
    {
        double nrm = 0.0;
        if (norm==1) { for (size_t l=0; l<L; ++l) { Y-=2; nrm += *Y; } }
        else if (norm==2) { for (size_t l=0; l<L; ++l) { Y-=2; nrm += *Y**Y; } nrm = sqrt(nrm); }
        else if (norm==3) { Y -= 2*(L-L/2); nrm = *Y; Y -= 2*(L/2); }
        for (size_t l=0; l<L; ++l, Y+=2) { *Y /= nrm; }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
