//Note that this is a continuation of a series of generalized symmetric cosine windows,
//with [Wikipedia]: a0 = 0.21557895, a1 = 0.41663158, a2 = 0.277263158, a3 = 0.083578947, a4 = 0.006947368.
//This window has negative values at the edges.

#include <stdio.h>
#include <math.h>

#ifndef M_PI
   #define M_PI 3.141592653589793238462643383279502884
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int flattop_s (float *Y, const size_t L, const size_t norm);
int flattop_d (double *Y, const size_t L, const size_t norm);
int flattop_c (float *Y, const size_t L, const size_t norm);
int flattop_z (double *Y, const size_t L, const size_t norm);


int flattop_s (float *Y, const size_t L, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in flattop_s: norm must be in {0,1,2,3}\n"); return 1; }

    const float p = (float)(2.0*M_PI/(double)(L-1u));
    
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = 0.21557895f - 0.41663158f*cosf((float)l*p) + 0.277263158f*cosf((float)(2u*l)*p) - 0.083578947f*cosf((float)(3u*l)*p) + 0.006947368f*cosf((float)(4u*l)*p); }
    if (L%2u) { *Y++ = 1.0f; }
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = *(Y-(2u*l+1u+L%2u)); }

    if (norm)
    {
        float nrm = 0.0f;
        if (norm==1u) { for (size_t l=0u; l<L; ++l) { --Y; nrm += fabsf(*Y); } }
        else if (norm==2u) { for (size_t l=0u; l<L; ++l) { --Y; nrm += *Y**Y; } nrm = sqrtf(nrm); }
        else if (norm==3u) { Y -= L-L/2u; nrm = *Y; Y -= L/2u; }
        for (size_t l=0u; l<L; ++l, ++Y) { *Y /= nrm; }
    }

    return 0;
}


int flattop_d (double *Y, const size_t L, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in flattop_d: norm must be in {0,1,2,3}\n"); return 1; }

    const double p = 2.0*M_PI/(double)(L-1u);
    
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = 0.21557895 - 0.41663158*cos((double)l*p) + 0.277263158*cos((double)(2u*l)*p) - 0.083578947*cos((double)(3u*l)*p) + 0.006947368*cos((double)(4u*l)*p); }
    if (L%2u) { *Y++ = 1.0; }
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = *(Y-(2u*l+1u+L%2u)); }

    if (norm)
    {
        double nrm = 0.0;
        if (norm==1u) { for (size_t l=0u; l<L; ++l) { --Y; nrm += fabs(*Y); } }
        else if (norm==2u) { for (size_t l=0u; l<L; ++l) { --Y; nrm += *Y**Y; } nrm = sqrt(nrm); }
        else if (norm==3u) { Y -= L-L/2u; nrm = *Y; Y -= L/2u; }
        for (size_t l=0u; l<L; ++l, ++Y) { *Y /= nrm; }
    }

    return 0;
}


int flattop_c (float *Y, const size_t L, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in flattop_c: norm must be in {0,1,2,3}\n"); return 1; }

    const float p = (float)(2.0*M_PI/(double)(L-1u));
    
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = 0.21557895f - 0.41663158f*cosf((float)l*p) + 0.277263158f*cosf((float)(2u*l)*p) - 0.083578947f*cosf((float)(3u*l)*p) + 0.006947368f*cosf((float)(4u*l)*p); *++Y = 0.0f; }
    if (L%2u) { *Y++ = 1.0f; *Y++ = 0.0f; }
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = *(Y-(4u*l+2u+2u*(L%2u))); *++Y = 0.0f; }

    if (norm)
    {
        float nrm = 0.0f;
        if (norm==1u) { for (size_t l=0u; l<L; ++l) { Y-=2u; nrm += fabsf(*Y); } }
        else if (norm==2u) { for (size_t l=0u; l<L; ++l) { Y-=2u; nrm += *Y**Y; } nrm = sqrtf(nrm); }
        else if (norm==3u) { Y -= 2u*(L-L/2u); nrm = *Y; Y -= 2u*(L/2u); }
        for (size_t l=0u; l<L; ++l, Y+=2u) { *Y /= nrm; }
    }

    return 0;
}


int flattop_z (double *Y, const size_t L, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in flattop_z: norm must be in {0,1,2,3}\n"); return 1; }

    const double p = 2.0*M_PI/(double)(L-1u);
    
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = 0.21557895 - 0.41663158*cos((double)l*p) + 0.277263158*cos((double)(2u*l)*p) - 0.083578947*cos((double)(3u*l)*p) + 0.006947368*cos((double)(4u*l)*p); *++Y = 0.0; }
    if (L%2u) { *Y++ = 1.0; *Y++ = 0.0; }
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = *(Y-(4u*l+2u+2u*(L%2u))); *++Y = 0.0; }

    if (norm)
    {
        double nrm = 0.0;
        if (norm==1u) { for (size_t l=0u; l<L; ++l) { Y-=2u; nrm += fabs(*Y); } }
        else if (norm==2u) { for (size_t l=0u; l<L; ++l) { Y-=2u; nrm += *Y**Y; } nrm = sqrt(nrm); }
        else if (norm==3u) { Y -= 2u*(L-L/2u); nrm = *Y; Y -= 2u*(L/2u); }
        for (size_t l=0u; l<L; ++l, Y+=2u) { *Y /= nrm; }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
