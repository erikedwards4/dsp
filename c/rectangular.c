//Makes rectangular window of length L

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int rectangular_s (float *Y, const size_t L, const size_t norm);
int rectangular_d (double *Y, const size_t L, const size_t norm);
int rectangular_c (float *Y, const size_t L, const size_t norm);
int rectangular_z (double *Y, const size_t L, const size_t norm);


int rectangular_s (float *Y, const size_t L, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in rectangular_s: norm must be in {0,1,2,3}\n"); return 1; }

    if (norm==0u || norm==3u)
    {
        for (size_t l=0u; l<L; ++l, ++Y) { *Y = 1.0f; }
    }
    else if (norm==1u)
    {
        const float y = 1.0f/(float)L;
        for (size_t l=0u; l<L; ++l, ++Y) { *Y = y; }
    }
    else if (norm==2u)
    {
        const float y = 1.0f/sqrtf((float)L);
        for (size_t l=0u; l<L; ++l, ++Y) { *Y = y; }
    }

    return 0;
}


int rectangular_d (double *Y, const size_t L, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in rectangular_d: norm must be in {0,1,2,3}\n"); return 1; }

    if (norm==0u || norm==3u)
    {
        for (size_t l=0u; l<L; ++l, ++Y) { *Y = 1.0; }
    }
    else if (norm==1u)
    {
        const double y = 1.0/(double)L;
        for (size_t l=0u; l<L; ++l, ++Y) { *Y = y; }
    }
    else if (norm==2u)
    {
        const double y = 1.0/sqrt((double)L);
        for (size_t l=0u; l<L; ++l, ++Y) { *Y = y; }
    }

    return 0;
}


int rectangular_c (float *Y, const size_t L, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in rectangular_c: norm must be in {0,1,2,3}\n"); return 1; }

    if (norm==0u || norm==3u)
    {
        for (size_t l=0u; l<L; ++l, ++Y) { *Y = 1.0f; *++Y = 0.0f; }
    }
    else if (norm==1u)
    {
        const float y = 1.0f/(float)L;
        for (size_t l=0u; l<L; ++l, ++Y) { *Y = y; *++Y = 0.0f; }
    }
    else if (norm==2u)
    {
        const float y = 1.0f/sqrtf((float)L);
        for (size_t l=0u; l<L; ++l, ++Y) { *Y = y; *++Y = 0.0f; }
    }

    return 0;
}


int rectangular_z (double *Y, const size_t L, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in rectangular_z: norm must be in {0,1,2,3}\n"); return 1; }

    if (norm==0u || norm==3u)
    {
        for (size_t l=0u; l<L; ++l, ++Y) { *Y = 1.0; *++Y = 0.0; }
    }
    else if (norm==1u)
    {
        const double y = 1.0/(double)L;
        for (size_t l=0u; l<L; ++l, ++Y) { *Y = y; *++Y = 0.0; }
    }
    else if (norm==2u)
    {
        const double y = 1.0/sqrt((double)L);
        for (size_t l=0u; l<L; ++l, ++Y) { *Y = y; *++Y = 0.0; }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif

