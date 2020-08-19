//Gauss window with length L and stdev param a = (L-1)/(2*stdev)

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int gauss_s (float *Y, const size_t L, const float r, const size_t norm);
int gauss_d (double *Y, const size_t L, const double r, const size_t norm);
int gauss_c (float *Y, const size_t L, const float r, const size_t norm);
int gauss_z (double *Y, const size_t L, const double r, const size_t norm);


int gauss_s (float *Y, const size_t L, const float a, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in gauss_s: norm must be in {0,1,2,3}\n"); return 1; }
    if (a<=0.0f) { fprintf(stderr,"error in gauss_s: a must be positive \n"); return 1; }

    const float p = -2.0f*(a*a/((L-1)*(L-1)));
    const float m = 0.5f*(L-1);

    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = expf((l-m)*(l-m)*p); }
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


int gauss_d (double *Y, const size_t L, const double a, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in gauss_d: norm must be in {0,1,2,3}\n"); return 1; }
    if (a<=0.0) { fprintf(stderr,"error in gauss_d: a must be positive \n"); return 1; }

    const double p = -2.0*(a*a/((L-1)*(L-1)));
    const double m = 0.5*(L-1);
    
    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = exp((l-m)*(l-m)*p); }
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


int gauss_c (float *Y, const size_t L, const float a, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in gauss_c: norm must be in {0,1,2,3}\n"); return 1; }
    if (a<=0.0f) { fprintf(stderr,"error in gauss_c: a must be positive \n"); return 1; }

    const float p = -2.0f*(a*a/((L-1)*(L-1)));
    const float m = 0.5f*(L-1);
    
    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = expf((l-m)*(l-m)*p); *++Y = 0.0f; }
    if (L%2) { *Y++ = 1.0f; *Y++ = 0.0f; }
    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = *(Y-2*l-1-L%2); *++Y = 0.0f; }

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


int gauss_z (double *Y, const size_t L, const double a, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in gauss_z: norm must be in {0,1,2,3}\n"); return 1; }
    if (a<=0.0) { fprintf(stderr,"error in gauss_z: a must be psoitive \n"); return 1; }
    
    const double p = -2.0*(a*a/((L-1)*(L-1)));
    const double m = 0.5*(L-1);
    
    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = exp((l-m)*(l-m)*p); *++Y = 0.0; }
    if (L%2) { *Y++ = 1.0; *Y++ = 0.0; }
    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = *(Y-2*l-1-L%2); *++Y = 0.0; }

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
