//There is an 'exact' blackman with:
//a0 = 7938/18608 ≈ 0.42659, a1 = 9240/18608 ≈ 0.49656, and a2 = 1430/18608 ≈ 0.076849
//However:
//a0 = 0.42, a1 = 0.5, a2 = 0.08
//are more-often used (they have improved 18 db/oct fall-off but do not null the sidelobes as well [Wikepedia]).
//However, this later window has negative values at the edges, so I include the exact version here as an option.

#include <stdio.h>
#include <math.h>

#ifndef M_PI
   #define M_PI 3.141592653589793238462643383279502884
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int blackman_s (float *Y, const size_t L, const int exact, const size_t norm);
int blackman_d (double *Y, const size_t L, const int exact, const size_t norm);
int blackman_c (float *Y, const size_t L, const int exact, const size_t norm);
int blackman_z (double *Y, const size_t L, const int exact, const size_t norm);


int blackman_s (float *Y, const size_t L, const int exact, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in blackman_s: norm must be in {0,1,2,3}\n"); return 1; }

    const float p = (float)(2.0*M_PI/(double)(L-1u));

    if (exact)
    {
        const float a0 = 7938.0f/18608.0f, a1 = 9240.0f/18608.0f, a2 = 1430.0f/18608.0f;
        for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = a0 - a1*cosf((float)l*p) + a2*cosf((float)(2u*l)*p); }
    }
    else
    {
        for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = 0.42f - 0.5f*cosf((float)l*p) + 0.08f*cosf((float)(2u*l)*p); }
    }
    if (L%2u) { *Y++ = 1.0f; }
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = *(Y-(2u*l+1u+L%2u)); }

    if (norm)
    {
        float nrm = 0.0f;
        if (norm==1u) { for (size_t l=0u; l<L; ++l) { --Y; nrm += *Y; } }
        else if (norm==2u) { for (size_t l=0u; l<L; ++l) { --Y; nrm += *Y**Y; } nrm = sqrtf(nrm); }
        else if (norm==3u) { Y -= L-L/2u; nrm = *Y; Y -= L/2u; }
        for (size_t l=0u; l<L; ++l, ++Y) { *Y /= nrm; }
    }

    return 0;
}


int blackman_d (double *Y, const size_t L, const int exact, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in blackman_d: norm must be in {0,1,2,3}\n"); return 1; }

    const double p = 2.0*M_PI/(double)(L-1u);
    
    if (exact)
    {
        const double a0 = 7938.0/18608.0, a1 = 9240.0/18608.0, a2 = 1430.0/18608.0;
        for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = a0 - a1*cos((double)l*p) + a2*cos((double)(2u*l)*p); }
    }
    else
    {
        for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = 0.42 - 0.5*cos((double)l*p) + 0.08*cos((double)(2u*l)*p); }
    }
    if (L%2u) { *Y++ = 1.0; }
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = *(Y-(2u*l+1u+L%2u)); }

    if (norm)
    {
        double nrm = 0.0;
        if (norm==1u) { for (size_t l=0u; l<L; ++l) { --Y; nrm += *Y; } }
        else if (norm==2u) { for (size_t l=0u; l<L; ++l) { --Y; nrm += *Y**Y; } nrm = sqrt(nrm); }
        else if (norm==3u) { Y -= L-L/2u; nrm = *Y; Y -= L/2u; }
        for (size_t l=0u; l<L; ++l, ++Y) { *Y /= nrm; }
    }

    return 0;
}


int blackman_c (float *Y, const size_t L, const int exact, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in blackman_c: norm must be in {0,1,2,3}\n"); return 1; }

    const float p = (float)(2.0*M_PI/(double)(L-1u));
    
    if (exact)
    {
        const float a0 = 7938.0f/18608.0f, a1 = 9240.0f/18608.0f, a2 = 1430.0f/18608.0f;
        for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = a0 - a1*cosf((float)l*p) + a2*cosf((float)(2u*l)*p); *++Y = 0.0f; }
    }
    else
    {
        for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = 0.42f - 0.5f*cosf((float)l*p) + 0.08f*cosf((float)(2u*l)*p); *++Y = 0.0f; }
    }
    if (L%2u) { *Y++ = 1.0f; *Y++ = 0.0f; }
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = *(Y-(4u*l+2u+2u*(L%2u))); *++Y = 0.0f; }

    if (norm)
    {
        float nrm = 0.0f;
        if (norm==1u) { for (size_t l=0u; l<L; ++l) { Y-=2u; nrm += *Y; } }
        else if (norm==2u) { for (size_t l=0u; l<L; ++l) { Y-=2u; nrm += *Y**Y; } nrm = sqrtf(nrm); }
        else if (norm==3u) { Y -= 2u*(L-L/2u); nrm = *Y; Y -= 2u*(L/2u); }
        for (size_t l=0u; l<L; ++l, Y+=2u) { *Y /= nrm; }
    }

    return 0;
}


int blackman_z (double *Y, const size_t L, const int exact, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in blackman_z: norm must be in {0,1,2,3}\n"); return 1; }

    const double p = 2.0*M_PI/(double)(L-1u);
    
    if (exact)
    {
        const double a0 = 7938.0/18608.0, a1 = 9240.0/18608.0, a2 = 1430.0/18608.0;
        for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = a0 - a1*cos((double)l*p) + a2*cos((double)(2u*l)*p); *++Y = 0.0; }
    }
    else
    {
        for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = 0.42 - 0.5*cos((double)l*p) + 0.08*cos((double)(2u*l)*p); *++Y = 0.0; }
    }
    if (L%2u) { *Y++ = 1.0; *Y++ = 0.0; }
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = *(Y-(4u*l+2u+2u*(L%2u))); *++Y = 0.0; }

    if (norm)
    {
        double nrm = 0.0;
        if (norm==1u) { for (size_t l=0u; l<L; ++l) { Y-=2u; nrm += *Y; } }
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
