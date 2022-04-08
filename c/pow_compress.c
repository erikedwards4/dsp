//Applies power compression to each element of X.
//This is intended for use with nonnegative X only (such as power).
//A small power regularization (preg) is also added.
//If pow is 1, then Y = X + preg.
//If pow is 0, then Y = log(X+preg)
//If 0<pow<1, Y = (X+preg).^pow;

#include <stdio.h>
#include <float.h>
#include <math.h>
#include "codee_dsp.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int pow_compress_s (float *Y, const float *X, const size_t N, const float p, const float preg)
{
    if (p<0.0f || p>1.0f) { fprintf(stderr,"error in pow_compress_s: pow must be in [0.0 1.0]\n"); return 1; }
    if (preg<0.0f) { fprintf(stderr,"error in pow_compress_s: preg must be nonnegative\n"); return 1; }

    if (p==1.0f)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X + preg; }
    }
    else if (p==0.0f)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = logf(*X+preg); }
    }
    else if (p==0.5f)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = sqrtf(*X+preg); }
    }
    else if (fabsf(p-1.0f/3.0f)<=FLT_EPSILON)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = cbrtf(*X+preg); }
    }
    else
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = powf(*X+preg,p); }
    }
    
    return 0;
}


int pow_compress_d (double *Y, const double *X, const size_t N, const double p, const double preg)
{
    if (p<0.0 || p>1.0) { fprintf(stderr,"error in pow_compress_d: pow must be in [0.0 1.0]\n"); return 1; }
    if (preg<0.0) { fprintf(stderr,"error in pow_compress_d: preg must be nonnegative\n"); return 1; }

    if (p==1.0)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X + preg; }
    }
    else if (p==0.0)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = log(*X+preg); }
    }
    else if (p==0.5)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = sqrt(*X+preg); }
    }
    else if (fabs(p-1.0/3.0)<=(double)FLT_EPSILON)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = cbrt(*X+preg); }
    }
    else
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = pow(*X+preg,p); }
    }
    
    return 0;
}


int pow_compress_inplace_s (float *X, const size_t N, const float p, const float preg)
{
    if (p<0.0f || p>1.0f) { fprintf(stderr,"error in pow_compress_inplace_s: pow must be in [0.0 1.0]\n"); return 1; }
    if (preg<0.0f) { fprintf(stderr,"error in pow_compress_inplace_s: preg must be nonnegative\n"); return 1; }

    if (p==1.0f)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X += preg; }
    }
    else if (p==0.0f)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = logf(*X+preg); }
    }
    else if (p==0.5f)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = sqrtf(*X+preg); }
    }
    else if (fabsf(p-1.0f/3.0f)<=FLT_EPSILON)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = cbrtf(*X+preg); }
    }
    else
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = powf(*X+preg,p); }
    }
    
    return 0;
}


int pow_compress_inplace_d (double*X, const size_t N, const double p, const double preg)
{
    if (p<0.0 || p>1.0) { fprintf(stderr,"error in pow_compress_inplace_d: pow must be in [0.0 1.0]\n"); return 1; }
    if (preg<0.0) { fprintf(stderr,"error in pow_compress_inplace_d: preg must be nonnegative\n"); return 1; }

    if (p==1.0)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X += preg; }
    }
    else if (p==0.0)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = log(*X+preg); }
    }
    else if (p==0.5)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = sqrt(*X+preg); }
    }
    else if (fabs(p-1.0/3.0)<=(double)FLT_EPSILON)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = cbrt(*X+preg); }
    }
    else
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = pow(*X+preg,p); }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
