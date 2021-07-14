//Gets real-valued power for complex-valued FFT.
//This just squares input X element-wise, that is:
// Y = |X|.^2 = Xr.*Xr + Xi.*Xi.

//This has in-place and not-in-place versions.

#include <stdio.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int fft_power_c (float *Y, const float *X, const size_t N);
int fft_power_z (double *Y, const double *X, const size_t N);


int fft_power_c (float *Y, const float *X, const size_t N)
{
    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    
    for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X**X; ++X; *Y += *X**X; }
    //for (size_t n=0u; n<N; ++n, X+=2, ++Y) { *Y = *X**X + *(X+1)**(X+1); }

    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    return 0;
}


int fft_power_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X**X; ++X; *Y += *X**X; }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
