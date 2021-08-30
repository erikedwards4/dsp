//Generates filter coefficients for smooth differentiators as described in:
//Holoborodko P. 2008. Smooth noise robust differentiators. www.holoborodko.com.
//See also: de Matos MC. 2018. Seismic attributes from the complex Teager-Kaiser energy.

//These can be used with fir.
//The same papers also give smooth_diffdiff filters.
//Together, these allow improved TKEO.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int smooth_diff_s (float *Y, const size_t N);
int smooth_diff_d (double *Y, const size_t N);


int smooth_diff_s (float *Y, const size_t N)
{

    return 0;
}


int smooth_diff_d (double *Y, const size_t N)
{

    return 0;
}


#ifdef __cplusplus
}
}
#endif
