//Gauss window with length L and stdev param a = (L-1u)/(2*stdev)

#include <stdio.h>
#include <math.h>
#include "codee_dsp.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int gauss_s (float *Y, const size_t L, const float a, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in gauss_s: norm must be in {0,1,2,3}\n"); return 1; }
    if (a<=0.0f) { fprintf(stderr,"error in gauss_s: a must be positive \n"); return 1; }

    const float p = -2.0f*(a*a/(float)((L-1u)*(L-1u)));
    const float m = 0.5f*(float)(L-1u);

    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = expf(((float)l-m)*((float)l-m)*p); }
    if (L%2u) { *Y++ = 1.0f; }
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = *(Y-(2u*l+1u+L%2u)); }

    if (norm)
    {
        float nrm = 0.0f;
        if (norm==1u) { for (size_t l=L; l>0u; --l) { --Y; nrm += *Y; } }
        else if (norm==2u) { for (size_t l=L; l>0u; --l) { --Y; nrm += *Y**Y; } nrm = sqrtf(nrm); }
        else if (norm==3u) { Y -= L-L/2u; nrm = *Y; Y -= L/2u; }
        for (size_t l=L; l>0u; --l, ++Y) { *Y /= nrm; }
    }

    return 0;
}


int gauss_d (double *Y, const size_t L, const double a, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in gauss_d: norm must be in {0,1,2,3}\n"); return 1; }
    if (a<=0.0) { fprintf(stderr,"error in gauss_d: a must be positive \n"); return 1; }

    const double p = -2.0*(a*a/(double)((L-1u)*(L-1u)));
    const double m = 0.5*(double)(L-1u);
    
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = exp(((double)l-m)*((double)l-m)*p); }
    if (L%2u) { *Y++ = 1.0; }
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = *(Y-(2u*l+1u+L%2u)); }

    if (norm)
    {
        double nrm = 0.0;
        if (norm==1u) { for (size_t l=L; l>0u; --l) { --Y; nrm += *Y; } }
        else if (norm==2u) { for (size_t l=L; l>0u; --l) { --Y; nrm += *Y**Y; } nrm = sqrt(nrm); }
        else if (norm==3u) { Y -= L-L/2u; nrm = *Y; Y -= L/2u; }
        for (size_t l=L; l>0u; --l, ++Y) { *Y /= nrm; }
    }

    return 0;
}


int gauss_c (float *Y, const size_t L, const float a, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in gauss_c: norm must be in {0,1,2,3}\n"); return 1; }
    if (a<=0.0f) { fprintf(stderr,"error in gauss_c: a must be positive \n"); return 1; }

    const float p = -2.0f*(a*a/(float)((L-1u)*(L-1u)));
    const float m = 0.5f*(float)(L-1u);
    
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = expf(((float)l-m)*((float)l-m)*p); *++Y = 0.0f; }
    if (L%2u) { *Y++ = 1.0f; *Y++ = 0.0f; }
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = *(Y-(2u*l+1u+L%2u)); *++Y = 0.0f; }

    if (norm)
    {
        float nrm = 0.0f;
        if (norm==1u) { for (size_t l=L; l>0u; --l) { Y-=2u; nrm += *Y; } }
        else if (norm==2u) { for (size_t l=L; l>0u; --l) { Y-=2u; nrm += *Y**Y; } nrm = sqrtf(nrm); }
        else if (norm==3u) { Y -= 2u*(L-L/2u); nrm = *Y; Y -= 2u*(L/2u); }
        for (size_t l=L; l>0u; --l, Y+=2u) { *Y /= nrm; }
    }

    return 0;
}


int gauss_z (double *Y, const size_t L, const double a, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in gauss_z: norm must be in {0,1,2,3}\n"); return 1; }
    if (a<=0.0) { fprintf(stderr,"error in gauss_z: a must be psoitive \n"); return 1; }
    
    const double p = -2.0*(a*a/(double)((L-1u)*(L-1u)));
    const double m = 0.5*(double)(L-1u);
    
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = exp(((double)l-m)*((double)l-m)*p); *++Y = 0.0; }
    if (L%2u) { *Y++ = 1.0; *Y++ = 0.0; }
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = *(Y-(2u*l+1u+L%2u)); *++Y = 0.0; }

    if (norm)
    {
        double nrm = 0.0;
        if (norm==1u) { for (size_t l=L; l>0u; --l) { Y-=2u; nrm += *Y; } }
        else if (norm==2u) { for (size_t l=L; l>0u; --l) { Y-=2u; nrm += *Y**Y; } nrm = sqrt(nrm); }
        else if (norm==3u) { Y -= 2u*(L-L/2u); nrm = *Y; Y -= 2u*(L/2u); }
        for (size_t l=L; l>0u; --l, Y+=2u) { *Y /= nrm; }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
