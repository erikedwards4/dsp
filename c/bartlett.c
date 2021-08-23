//For Barlett, I use Octave/Matlab convention,
//for which they cite Openheim & Shafer: Discrete-Time Signal Processing.

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int bartlett_s (float *Y, const size_t L, const size_t norm);
int bartlett_d (double *Y, const size_t L, const size_t norm);
int bartlett_c (float *Y, const size_t L, const size_t norm);
int bartlett_z (double *Y, const size_t L, const size_t norm);


int bartlett_s (float *Y, const size_t L, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in bartlett_s: norm must be in {0,1,2,3}\n"); return 1; }

    const float p = 1.0f/(float)(L-(L%2u)*(L/2u)-1u);
    // struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = (float)(l*(2u-L%2u))*p; }
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

    // clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


int bartlett_d (double *Y, const size_t L, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in bartlett_d: norm must be in {0,1,2,3}\n"); return 1; }
    
    const double p = 1.0/(double)(L-(L%2u)*(L/2u)-1u);

    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = (double)(l*(2u-L%2u))*p; }
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


int bartlett_c (float *Y, const size_t L, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in bartlett_c: norm must be in {0,1,2,3}\n"); return 1; }

    const float p = 1.0f/(float)(L-(L%2u)*(L/2u)-1u);

    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = (float)(l*(2u-L%2u))*p; *++Y = 0.0f; }
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


int bartlett_z (double *Y, const size_t L, const size_t norm)
{
    if (norm>3u) { fprintf(stderr,"error in bartlett_z: norm must be in {0,1,2,3}\n"); return 1; }

    const double p = 1.0/(double)(L-(L%2u)*(L/2u)-1u);
    
    for (size_t l=0u; l<L/2u; ++l, ++Y) { *Y = (double)(l*(2u-L%2u))*p; *++Y = 0.0; }
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
