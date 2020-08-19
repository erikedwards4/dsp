//Note that this is a continuation of a series of generalized symmetric cosine windows,
//with [Wikipedia]: a0 = 0.35875, a1 = 0.48829, a2 = 0.14128, a3 = 0.01168;

#include <stdio.h>
#include <math.h>

#ifndef M_PI
   #define M_PI 3.141592653589793238462643383279502884
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int blackmanharris_s (float *Y, const size_t L, const size_t norm);
int blackmanharris_d (double *Y, const size_t L, const size_t norm);
int blackmanharris_c (float *Y, const size_t L, const size_t norm);
int blackmanharris_z (double *Y, const size_t L, const size_t norm);


int blackmanharris_s (float *Y, const size_t L, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in blackmanharris_s: norm must be in {0,1,2,3}\n"); return 1; }

    const float p = (float)(2.0*M_PI/(L-1));

    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = 0.35875f - 0.48829f*cosf(l*p) + 0.14128f*cosf(2*l*p) - 0.01168f*cosf(3*l*p); }
    if (L%2) { *Y++ = 1.0f; }
    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = *(Y-2*l-1-L%2); }

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


int blackmanharris_d (double *Y, const size_t L, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in blackmanharris_d: norm must be in {0,1,2,3}\n"); return 1; }

    const double p = 2.0*M_PI/(L-1.0);
    
    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = 0.35875 - 0.48829*cos(l*p) + 0.14128*cos(2*l*p) - 0.01168*cos(3*l*p); }
    if (L%2) { *Y++ = 1.0; }
    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = *(Y-2*l-1-L%2); }

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


int blackmanharris_c (float *Y, const size_t L, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in blackmanharris_c: norm must be in {0,1,2,3}\n"); return 1; }

    const float p = (float)(2.0*M_PI/(L-1));
    
    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = 0.35875f - 0.48829f*cosf(l*p) + 0.14128f*cosf(2*l*p) - 0.01168f*cosf(3*l*p); *++Y = 0.0f; }
    if (L%2) { *Y++ = 1.0f; *Y++ = 0.0f; }
    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = *(Y-4*l-2-2*(L%2)); *++Y = 0.0f; }

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


int blackmanharris_z (double *Y, const size_t L, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in blackmanharris_z: norm must be in {0,1,2,3}\n"); return 1; }

    const double p = 2.0*M_PI/(L-1.0);
    
    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = 0.35875 - 0.48829*cos(l*p) + 0.14128*cos(2*l*p) - 0.01168*cos(3*l*p); *++Y = 0.0; }
    if (L%2) { *Y++ = 1.0; *Y++ = 0.0; }
    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = *(Y-4*l-2-2*(L%2)); *++Y = 0.0; }

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
