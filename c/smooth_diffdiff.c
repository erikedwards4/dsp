//Generates filter coefficients for smooth differentiators as described in:
//Holoborodko P. 2008. Smooth noise robust differentiators. www.holoborodko.com.
//See also: de Matos MC. 2018. Seismic attributes from the complex Teager-Kaiser energy.

//These can be used with fir.
//See smooth_diff for more comments.

//N is the (odd) filter-length, and n is the degree of the monomial (1+x+x^2+...+x^n).

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int smooth_diffdiff_s (float *Y, const size_t N, const size_t n);
int smooth_diffdiff_d (double *Y, const size_t N, const size_t n);


int smooth_diffdiff_s (float *Y, const size_t N, const size_t n)
{
    if (N%2u==0u) { fprintf(stderr,"error in smooth_diffdiff_s: N (filter length) must be odd\n"); return 1; }
    if (n%2u==0u) { fprintf(stderr,"error in smooth_diffdiff_s: n (degree of polynomial exactness) must be odd\n"); return 1; }
    
    if (n==3u)
    {
        if (N==5u)
        {
            *Y++ = 1.0f/4.0f; *Y++ = 0.0f;
            *Y++ = -2.0f/4.0f;
            *Y++ = 0.0f; *Y++ = 1.0f/4.0f;
        }
        else if (N==7u)
        {
            *Y++ = 1.0f/16.0f; *Y++ = 2.0f/16.0f; *Y++ = -1.0f/16.0f;
            *Y++ = -4.0f/16.0f;
            *Y++ = -1.0f/16.0f; *Y++ = 2.0f/16.0f; *Y++ = 1.0f/16.0f;
        }
        else if (N==9u)
        {
            *Y++ = 1.0f/64.0f; *Y++ = 4.0f/64.0f; *Y++ = -4.0f/64.0f; *Y++ = -4.0f/64.0f;
            *Y++ = -10.0f/64.0f;
            *Y++ = -4.0f/64.0f; *Y++ = -4.0f/64.0f; *Y++ = 4.0f/64.0f; *Y++ = 1.0f/64.0f;
        }
        else
        {
            fprintf(stderr,"error in smooth_diffdiff_s: N (filter length) must be in {5,7,9} for n=3\n"); return 1;
        }
    }
    else if (n==5u)
    {
        if (N==7u)
        {
            *Y++ = -1.0f/12.0f; *Y++ = 5.0f/12.0f; *Y++ = 1.0f/12.0f;
            *Y++ = -10.0f/12.0f;
            *Y++ = 1.0f/12.0f; *Y++ = 5.0f/12.0f; *Y++ = -1.0f/12.0f;
        }
        else if (N==9u)
        {
            *Y++ = -7.0f/192.0f; *Y++ = 12.0f/192.0f; *Y++ = 52.0f/192.0f; *Y++ = -12.0f/192.0f;
            *Y++ = -90.0f/192.0f;
            *Y++ = -12.0f/192.0f; *Y++ = 52.0f/192.0f; *Y++ = 12.0f/192.0f; *Y++ = -7.0f/192.0f;
        }
        else
        {
            fprintf(stderr,"error in smooth_diffdiff_s: N (filter length) must be in {7,9} for n=5\n"); return 1;
        }
    }
    else
    {
        fprintf(stderr,"error in smooth_diffdiff_s: n (degree of polynomial exactness) must be 3 or 5\n"); return 1;
    }

    return 0;
}


int smooth_diffdiff_d (double *Y, const size_t N, const size_t n)
{
    if (N%2u==0u) { fprintf(stderr,"error in smooth_diffdiff_d: N (filter length) must be odd\n"); return 1; }
    if (n%2u==0u) { fprintf(stderr,"error in smooth_diffdiff_d: n (degree of polynomial exactness) must be odd\n"); return 1; }
    
    if (n==3u)
    {
        if (N==5u)
        {
            *Y++ = 1.0/4.0; *Y++ = 0.0;
            *Y++ = -2.0/4.0;
            *Y++ = 0.0; *Y++ = 1.0/4.0;
        }
        else if (N==7u)
        {
            *Y++ = 1.0/16.0; *Y++ = 2.0/16.0; *Y++ = -1.0/16.0;
            *Y++ = -4.0/16.0;
            *Y++ = -1.0/16.0; *Y++ = 2.0/16.0; *Y++ = 1.0/16.0;
        }
        else if (N==9u)
        {
            *Y++ = 1.0/64.0; *Y++ = 4.0/64.0; *Y++ = -4.0/64.0; *Y++ = -4.0/64.0;
            *Y++ = -10.0/64.0;
            *Y++ = -4.0/64.0; *Y++ = -4.0/64.0; *Y++ = 4.0/64.0; *Y++ = 1.0/64.0;
        }
        else
        {
            fprintf(stderr,"error in smooth_diffdiff_d: N (filter length) must be in {5,7,9} for n=3\n"); return 1;
        }
    }
    else if (n==5u)
    {
        if (N==7u)
        {
            *Y++ = -1.0/12.0; *Y++ = 5.0/12.0; *Y++ = 1.0/12.0;
            *Y++ = -10.0/12.0;
            *Y++ = 1.0/12.0; *Y++ = 5.0/12.0; *Y++ = -1.0/12.0;
        }
        else if (N==9u)
        {
            *Y++ = -7.0/192.0; *Y++ = 12.0/192.0; *Y++ = 52.0/192.0; *Y++ = -12.0/192.0;
            *Y++ = -90.0/192.0;
            *Y++ = -12.0/192.0; *Y++ = 52.0/192.0; *Y++ = 12.0/192.0; *Y++ = -7.0/192.0;
        }
        else
        {
            fprintf(stderr,"error in smooth_diffdiff_d: N (filter length) must be in {7,9} for n=5\n"); return 1;
        }
    }
    else
    {
        fprintf(stderr,"error in smooth_diffdiff_d: n (degree of polynomial exactness) must be 3 or 5\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
