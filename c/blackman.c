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

int blackman_s (float *Y, const size_t L, const char exact, const size_t norm);
int blackman_d (double *Y, const size_t L, const char exact, const size_t norm);
int blackman_c (float *Y, const size_t L, const char exact, const size_t norm);
int blackman_z (double *Y, const size_t L, const char exact, const size_t norm);


int blackman_s (float *Y, const size_t L, const char exact, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in blackman_s: norm must be in {0,1,2,3}\n"); return 1; }

    const float p = (float)(2.0*M_PI/(L-1));

    if (exact)
    {
        const float a0 = 7938.0f/18608.0f, a1 = 9240.0f/18608.0f, a2 = 1430.0f/18608.0f;
        for (size_t l=0; l<L/2; ++l, ++Y) { *Y = a0 - a1*cosf(l*p) + a2*cosf(2*l*p); }
    }
    else
    {
        for (size_t l=0; l<L/2; ++l, ++Y) { *Y = 0.42f - 0.5f*cosf(l*p) + 0.08f*cosf(2*l*p); }
    }
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


int blackman_d (double *Y, const size_t L, const char exact, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in blackman_d: norm must be in {0,1,2,3}\n"); return 1; }

    const double p = 2.0*M_PI/(L-1.0);
    
    if (exact)
    {
        const double a0 = 7938.0/18608.0, a1 = 9240.0/18608.0, a2 = 1430.0/18608.0;
        for (size_t l=0; l<L/2; ++l, ++Y) { *Y = a0 - a1*cos(l*p) + a2*cos(2*l*p); }
    }
    else
    {
        for (size_t l=0; l<L/2; ++l, ++Y) { *Y = 0.42 - 0.5*cos(l*p) + 0.08*cos(2*l*p); }
    }
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


int blackman_c (float *Y, const size_t L, const char exact, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in blackman_c: norm must be in {0,1,2,3}\n"); return 1; }

    const float p = (float)(2.0*M_PI/(L-1));
    
    if (exact)
    {
        const float a0 = 7938.0f/18608.0f, a1 = 9240.0f/18608.0f, a2 = 1430.0f/18608.0f;
        for (size_t l=0; l<L/2; ++l, ++Y) { *Y = a0 - a1*cosf(l*p) + a2*cosf(2*l*p); *++Y = 0.0f; }
    }
    else
    {
        for (size_t l=0; l<L/2; ++l, ++Y) { *Y = 0.42f - 0.5f*cosf(l*p) + 0.08f*cosf(2*l*p); *++Y = 0.0f; }
    }
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


int blackman_z (double *Y, const size_t L, const char exact, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in blackman_z: norm must be in {0,1,2,3}\n"); return 1; }

    const double p = 2.0*M_PI/(L-1.0);
    
    if (exact)
    {
        const double a0 = 7938.0/18608.0, a1 = 9240.0/18608.0, a2 = 1430.0/18608.0;
        for (size_t l=0; l<L/2; ++l, ++Y) { *Y = a0 - a1*cos(l*p) + a2*cos(2*l*p); *++Y = 0.0; }
    }
    else
    {
        for (size_t l=0; l<L/2; ++l, ++Y) { *Y = 0.42 - 0.5*cos(l*p) + 0.08*cos(2*l*p); *++Y = 0.0; }
    }
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
