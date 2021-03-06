//Tukey (tapered cosine) window with length L and cosine fraction r

#include <stdio.h>
#include <math.h>

#ifndef M_PI
   #define M_PI 3.141592653589793238462643383279502884
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int tukey_s (float *Y, const size_t L, const float r, const size_t norm);
int tukey_d (double *Y, const size_t L, const double r, const size_t norm);
int tukey_c (float *Y, const size_t L, const float r, const size_t norm);
int tukey_z (double *Y, const size_t L, const double r, const size_t norm);


int tukey_s (float *Y, const size_t L, const float r, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in tukey_s: norm must be in {0,1,2,3}\n"); return 1; }
    if (r<0.0f || r>1.0f) { fprintf(stderr,"error in tukey_s: r must be in [0.0 1.0] \n"); return 1; }

    if (L==0) {}
    else if (L==1) { *Y++ = 1.0f; }
    else if (L==2) { *Y++ = 1.0f; *Y++ = 1.0f; }
    else
    {
        const float p = (float)(2.0*M_PI/((L-1)*(double)r));
        for (size_t l=0; l<0.5f*r*L; ++l, ++Y) { *Y = 0.5f - 0.5f*cosf(l*p); }
        for (size_t l=(size_t)ceilf(0.5f*r*L); l<L/2+L%2; ++l, ++Y) { *Y = 1.0f; }
        for (size_t l=0; l<L/2; ++l, ++Y) { *Y = *(Y-2*l-1-L%2); }
    }

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


int tukey_d (double *Y, const size_t L, const double r, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in tukey_d: norm must be in {0,1,2,3}\n"); return 1; }
    if (r<0.0 || r>1.0) { fprintf(stderr,"error in tukey_d: r must be in [0.0 1.0] \n"); return 1; }

    if (L==0) {}
    else if (L==1) { *Y++ = 1.0; }
    else if (L==2) { *Y++ = 1.0; *Y++ = 1.0; }
    else
    {
        const double p = 2.0*M_PI/((L-1)*r);
        for (size_t l=0; l<0.5*r*L; ++l, ++Y) { *Y = 0.5 - 0.5*cos(l*p); }
        for (size_t l=(size_t)ceil(0.5*r*L); l<L/2+L%2; ++l, ++Y) { *Y = 1.0; }
        for (size_t l=0; l<L/2; ++l, ++Y) { *Y = *(Y-2*l-1-L%2); }
    }

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


int tukey_c (float *Y, const size_t L, const float r, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in tukey_c: norm must be in {0,1,2,3}\n"); return 1; }
    if (r<0.0f || r>1.0f) { fprintf(stderr,"error in tukey_c: r must be in [0.0 1.0] \n"); return 1; }

    if (L==0) {}
    else if (L==1) { *Y++ = 1.0; *Y++ = 0.0; }
    else if (L==2) { *Y++ = 1.0; *Y++ = 0.0; *Y++ = 1.0; *Y++ = 0.0; }
    else
    {
        const float p = (float)(2.0*M_PI/((L-1)*(double)r));
        for (size_t l=0; l<0.5f*r*L; ++l, ++Y) { *Y = 0.5f - 0.5f*cosf(l*p); *++Y = 0.0f; }
        for (size_t l=(size_t)ceilf(0.5f*r*L); l<L/2+L%2; ++l, ++Y) { *Y = 1.0f; *++Y = 0.0f; }
        for (size_t l=0; l<L/2; ++l, ++Y) { *Y = *(Y-4*l-2-2*(L%2)); *++Y = 0.0f; }
    }

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


int tukey_z (double *Y, const size_t L, const double r, const size_t norm)
{
    if (norm>3) { fprintf(stderr,"error in tukey_z: norm must be in {0,1,2,3}\n"); return 1; }
    if (r<0.0 || r>1.0) { fprintf(stderr,"error in tukey_z: r must be in [0.0 1.0] \n"); return 1; }

    if (L==0) {}
    else if (L==1) { *Y++ = 1.0; *Y++ = 0.0; }
    else if (L==2) { *Y++ = 1.0; *Y++ = 0.0; *Y++ = 1.0; *Y++ = 0.0; }
    else
    {
        const double p = 2.0*M_PI/((L-1)*r);
        for (size_t l=0; l<0.5*r*L; ++l, ++Y) { *Y = 0.5 - 0.5*cos(l*p); *++Y = 0.0; }
        for (size_t l=(size_t)ceil(0.5*r*L); l<L/2+L%2; ++l, ++Y) { *Y = 1.0; *++Y = 0.0; }
        for (size_t l=0; l<L/2; ++l, ++Y) { *Y = *(Y-4*l-2-2*(L%2)); *++Y = 0.0; }
    }

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
