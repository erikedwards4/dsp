//For triangular, I used Wikipedia definition, but this leads a vanishing 1st sample.
//So, instead I just follow the Octave/Matlab convention.

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int triangular_s (float *Y, const size_t L, const size_t norm);
int triangular_d (double *Y, const size_t L, const size_t norm);
int triangular_c (float *Y, const size_t L, const size_t norm);
int triangular_z (double *Y, const size_t L, const size_t norm);


int triangular_s (float *Y, const size_t L, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in triangular_s: norm must be in {0,1,2,3}\n"); return 1; }

    const float p = 1.0f/(L-(L%2)*(L/2));

    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = (l*(2-L%2)+1)*p; }
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


int triangular_d (double *Y, const size_t L, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in triangular_d: norm must be in {0,1,2,3}\n"); return 1; }

    const double p = 1.0/(L-(L%2)*(L/2));
    
    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = (l*(2-L%2)+1)*p; }
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


int triangular_c (float *Y, const size_t L, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in triangular_c: norm must be in {0,1,2,3}\n"); return 1; }

    const float p = 1.0f/(L-(L%2)*(L/2));

    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = (l*(2-L%2)+1)*p; *++Y = 0.0f; }
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


int triangular_z (double *Y, const size_t L, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in triangular_z: norm must be in {0,1,2,3}\n"); return 1; }

    const double p = 1.0/(L-(L%2)*(L/2));
    
    for (size_t l=0; l<L/2; ++l, ++Y) { *Y = (l*(2-L%2)+1)*p; *++Y = 0.0; }
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
