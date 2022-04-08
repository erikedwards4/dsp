//Generates filter coefficients for smooth differentiators as described in:
//Holoborodko P. 2008. Smooth noise robust differentiators. www.holoborodko.com.
//See also: de Matos MC. 2018. Seismic attributes from the complex Teager-Kaiser energy.

//These can be used with fir.
//The same papers also give smooth_diffdiff filters.
//It also mentions that causal, 2D and bandpass differentiators are possible.

//Using the smooth diff and diffdiff allow improved TKEO, as used in de Matos [2018].
//[The later paper also describes the variational version of the TK (TKV).]

//These are derived for the following properties:
//exactness on polynomials
//preciseness on low frequencies
//noise-robust (suppression of high freqs)

//N is the (odd) filter-length, and n is the degree of the monomial (1+x+x^2+...+x^n).

#include <stdio.h>
#include "codee_dsp.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int smooth_diff_s (float *Y, const size_t N, const size_t n)
{
    if (N%2u==0u) { fprintf(stderr,"error in smooth_diff_s: N (filter length) must be odd\n"); return 1; }
    if (n%2u==1u) { fprintf(stderr,"error in smooth_diff_s: n (degree of polynomial exactness) must be even\n"); return 1; }
    
    if (n==2u)
    {
        if (N==5u)
        {
            *Y++ = -1.0f/8.0f; *Y++ = -2.0f/8.0f;
            *Y++ = 0.0f;
            *Y++ = 2.0f/8.0f; *Y++ = 1.0f/8.0f;
        }
        else if (N==7u)
        {
            *Y++ = -1.0f/32.0f; *Y++ = -4.0f/32.0f; *Y++ = -5.0f/32.0f;
            *Y++ = 0.0f;
            *Y++ = 5.0f/32.0f; *Y++ = 4.0f/32.0f; *Y++ = 1.0f/32.0f;
        }
        else if (N==9u)
        {
            *Y++ = -1.0f/128.0f; *Y++ = -6.0f/128.0f; *Y++ = -14.0f/128.0f; *Y++ = -14.0f/128.0f;
            *Y++ = 0.0f;
            *Y++ = 14.0f/128.0f; *Y++ = 14.0f/128.0f; *Y++ = 6.0/128.0f; *Y++ = 1.0f/128.0f;
        }
        else if (N==11u)
        {
            *Y++ = -1.0f/512.0f; *Y++ = -8.0f/512.0f; *Y++ = -27.0f/512.0f; *Y++ = -48.0f/512.0f; *Y++ = -42.0f/512.0f;
            *Y++ = 0.0f;
            *Y++ = 42.0f/512.0f; *Y++ = 48.0f/512.0f; *Y++ = 27.0/512.0f; *Y++ = 8.0f/512.0f; *Y++ = 1.0f/512.0f;
        }
        else
        {
            fprintf(stderr,"error in smooth_diff_s: N (filter length) must be in {5,7,9,11} for n=2\n"); return 1;
        }
    }
    else if (n==4u)
    {
        if (N==7u)
        {
            *Y++ = 5.0f/96.0f; *Y++ = -12.0f/96.0f; *Y++ = -39.0f/96.0f;
            *Y++ = 0.0f;
            *Y++ = 39.0f/96.0f; *Y++ = 12.0f/96.0f; *Y++ = -5.0f/96.0f;
        }
        else if (N==9u)
        {
            *Y++ = 2.0f/96.0f; *Y++ = 1.0f/96.0f; *Y++ = -16.0f/96.0f; *Y++ = -27.0f/96.0f;
            *Y++ = 0.0f;
            *Y++ = 27.0f/96.0f; *Y++ = 16.0f/96.0f; *Y++ = -1.0f/96.0f; *Y++ = -2.0f/96.0f;
        }
        else if (N==11u)
        {
            *Y++ = 11.0f/1536.0f; *Y++ = 32.0f/1536.0f; *Y++ = -39.0f/1536.0f; *Y++ = -256.0f/1536.0f; *Y++ = -322.0f/1536.0f;
            *Y++ = 0.0f;
            *Y++ = 322.0f/1536.0f; *Y++ = 256.0f/1536.0f; *Y++ = 39.0f/1536.0f; *Y++ = -32.0f/1536.0f; *Y++ = -11.0f/1536.0f;
        }
        else
        {
            fprintf(stderr,"error in smooth_diff_s: N (filter length) must be in {7,9,11} for n=4\n"); return 1;
        }
    }
    else
    {
        fprintf(stderr,"error in smooth_diff_s: n (degree of polynomial exactness) must be 2 or 4\n"); return 1;
    }

    return 0;
}


int smooth_diff_d (double *Y, const size_t N, const size_t n)
{
    if (N%2u==0u) { fprintf(stderr,"error in smooth_diff_d: N (filter length) must be odd\n"); return 1; }
    if (n%2u==1u) { fprintf(stderr,"error in smooth_diff_d: n (degree of polynomial exactness) must be even\n"); return 1; }
    
    if (n==2u)
    {
        if (N==5u)
        {
            *Y++ = -1.0/8.0; *Y++ = -2.0/8.0;
            *Y++ = 0.0;
            *Y++ = 2.0/8.0; *Y++ = 1.0/8.0;
        }
        else if (N==7u)
        {
            *Y++ = -1.0/32.0; *Y++ = -4.0/32.0; *Y++ = -5.0/32.0;
            *Y++ = 0.0;
            *Y++ = 5.0/32.0; *Y++ = 4.0/32.0; *Y++ = 1.0/32.0;
        }
        else if (N==9u)
        {
            *Y++ = -1.0/128.0; *Y++ = -6.0/128.0; *Y++ = -14.0/128.0; *Y++ = -14.0/128.0;
            *Y++ = 0.0;
            *Y++ = 14.0/128.0; *Y++ = 14.0/128.0; *Y++ = 6.0/128.0; *Y++ = 1.0/128.0;
        }
        else if (N==11u)
        {
            *Y++ = -1.0/512.0; *Y++ = -8.0/512.0; *Y++ = -27.0/512.0; *Y++ = -48.0/512.0; *Y++ = -42.0/512.0;
            *Y++ = 0.0;
            *Y++ = 42.0/512.0; *Y++ = 48.0/512.0; *Y++ = 27.0/512.0; *Y++ = 8.0/512.0; *Y++ = 1.0/512.0;
        }
        else
        {
            fprintf(stderr,"error in smooth_diff_d: N (filter length) must be in {5,7,9,11} for n=2\n"); return 1;
        }
    }
    else if (n==4u)
    {
        if (N==7u)
        {
            *Y++ = 5.0/96.0; *Y++ = -12.0/96.0; *Y++ = -39.0/96.0;
            *Y++ = 0.0;
            *Y++ = 39.0/96.0; *Y++ = 12.0/96.0; *Y++ = -5.0/96.0;
        }
        else if (N==9u)
        {
            *Y++ = 2.0/96.0; *Y++ = 1.0/96.0; *Y++ = -16.0/96.0; *Y++ = -27.0/96.0;
            *Y++ = 0.0;
            *Y++ = 27.0/96.0; *Y++ = 16.0/96.0; *Y++ = -1.0/96.0; *Y++ = -2.0/96.0;
        }
        else if (N==11u)
        {
            *Y++ = 11.0/1536.0; *Y++ = 32.0/1536.0; *Y++ = -39.0/1536.0; *Y++ = -256.0/1536.0; *Y++ = -322.0/1536.0;
            *Y++ = 0.0;
            *Y++ = 322.0/1536.0; *Y++ = 256.0/1536.0; *Y++ = 39.0/1536.0; *Y++ = -32.0/1536.0; *Y++ = -11.0/1536.0;
        }
        else
        {
            fprintf(stderr,"error in smooth_diff_d: N (filter length) must be in {7,9,11} for n=4\n"); return 1;
        }
    }
    else
    {
        fprintf(stderr,"error in smooth_diff_d: n (degree of polynomial exactness) must be 2 or 4\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
