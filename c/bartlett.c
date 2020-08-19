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
    if (norm>3) { fprintf(stderr,"error in bartlett_s: norm must be in {0,1,2,3}\n"); return 1; }

    const float p = 1.0f/(L-(L%2)*(L/2)-1);
    // struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = (l*(2-L%2))*p; }
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

    // clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


int bartlett_d (double *Y, const size_t L, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in bartlett_d: norm must be in {0,1,2,3}\n"); return 1; }
    
    const double p = 1.0/(L-(L%2)*(L/2)-1);

    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = (l*(2-L%2))*p; }
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


int bartlett_c (float *Y, const size_t L, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in bartlett_c: norm must be in {0,1,2,3}\n"); return 1; }

    const float p = 1.0f/(L-(L%2)*(L/2)-1);

    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = (l*(2-L%2))*p; *++Y = 0.0f; }
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


int bartlett_z (double *Y, const size_t L, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in bartlett_z: norm must be in {0,1,2,3}\n"); return 1; }

    const double p = 1.0/(L-(L%2)*(L/2)-1);
    
    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = (l*(2-L%2))*p; *++Y = 0.0; }
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
