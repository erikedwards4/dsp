//Generates filter coefficients for Apply Spencer's 15-point MA filter.
//B = [-3, -6, -5, 3, 21, 46, 67, 74, 67, 46, 21, 3, -5, -6, -3] / 320;
//These can be used with fir.

//This is a symmetric filter that is designed to allow all
//linear, quadratic and cubic functions to pass through unaltered.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int spencer_s (float *Y);
int spencer_d (double *Y);


int spencer_s (float *Y)
{
    *Y++ = -3.0f / 320.0f;
    *Y++ = -6.0f / 320.0f;
    *Y++ = -5.0f / 320.0f;
    *Y++ =  3.0f / 320.0f;
    *Y++ = 21.0f / 320.0f;
    *Y++ = 46.0f / 320.0f;
    *Y++ = 67.0f / 320.0f;
    *Y++ = 74.0f / 320.0f;
    *Y++ = 67.0f / 320.0f;
    *Y++ = 46.0f / 320.0f;
    *Y++ = 21.0f / 320.0f;
    *Y++ =  3.0f / 320.0f;
    *Y++ = -5.0f / 320.0f;
    *Y++ = -6.0f / 320.0f;
    *Y   = -3.0f / 320.0f;

    return 0;
}


int spencer_d (double *Y)
{
    *Y++ = -3.0 / 320.0;
    *Y++ = -6.0 / 320.0;
    *Y++ = -5.0 / 320.0;
    *Y++ =  3.0 / 320.0;
    *Y++ = 21.0 / 320.0;
    *Y++ = 46.0 / 320.0;
    *Y++ = 67.0 / 320.0;
    *Y++ = 74.0 / 320.0;
    *Y++ = 67.0 / 320.0;
    *Y++ = 46.0 / 320.0;
    *Y++ = 21.0 / 320.0;
    *Y++ =  3.0 / 320.0;
    *Y++ = -5.0 / 320.0;
    *Y++ = -6.0 / 320.0;
    *Y   = -3.0 / 320.0;

    return 0;
}


#ifdef __cplusplus
}
}
#endif
