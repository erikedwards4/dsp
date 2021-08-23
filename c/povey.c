//From Kaldi: "like hamming but goes to zero at edges".

#include <stdio.h>
#include <math.h>

#ifndef M_PI
   #define M_PI 3.141592653589793238462643383279502884
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int povey_s (float *Y, const size_t L, const size_t norm);
int povey_d (double *Y, const size_t L, const size_t norm);
int povey_c (float *Y, const size_t L, const size_t norm);
int povey_z (double *Y, const size_t L, const size_t norm);


int povey_s (float *Y, const size_t L, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in povey_s: norm must be in {0,1,2,3}\n"); return 1; }

    const float p = (float)(2.0*M_PI/(double)(L-1u));
    
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = powf(0.5f-0.5f*cosf((float)l*p),0.85f); }
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


int povey_d (double *Y, const size_t L, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in povey_d: norm must be in {0,1,2,3}\n"); return 1; }

    const double p = 2.0*M_PI/(double)(L-1u);
    
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = pow(0.5-0.5*cos((double)l*p),0.85); }
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


int povey_c (float *Y, const size_t L, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in povey_c: norm must be in {0,1,2,3}\n"); return 1; }

    const float p = (float)(2.0*M_PI/(double)(L-1u));
    
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = powf(0.5f-0.5f*cosf((float)l*p),0.85f); *++Y = 0.0f; }
    if (L%2u) { *Y++ = 1.0f; *Y++ = 0.0f; }
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = *(Y-(4u*l+2u+2u*(L%2u))); *++Y = 0.0f; }

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


int povey_z (double *Y, const size_t L, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in povey_z: norm must be in {0,1,2,3}\n"); return 1; }
    
    const double p = 2.0*M_PI/(double)(L-1u);
    
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = pow(0.5-0.5*cos((double)l*p),0.85); *++Y = 0.0; }
    if (L%2u) { *Y++ = 1.0; *Y++ = 0.0; }
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = *(Y-(4u*l+2u+2u*(L%2u))); *++Y = 0.0; }

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
