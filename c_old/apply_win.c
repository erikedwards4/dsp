//Applies window (X2) to input tensor (X2) along dim
//This is a 2-input elementwise function.
//Does elementwise multiplication: Y = X1*X2,
//such that Y has the same size as X1.
//X2 is a vector of length L, and X1 must have length L along dim.
//This has in-place and not-in-place versions.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int apply_win_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const size_t dim, const char iscolmajor);
int apply_win_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const size_t dim, const char iscolmajor);
int apply_win_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const size_t dim, const char iscolmajor);
int apply_win_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const size_t dim, const char iscolmajor);

int apply_win_inplace_s (float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const size_t dim, const char iscolmajor);
int apply_win_inplace_d (double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const size_t dim, const char iscolmajor);
int apply_win_inplace_c (float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const size_t dim, const char iscolmajor);
int apply_win_inplace_z (double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const size_t dim, const char iscolmajor);


int apply_win_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const size_t dim, const char iscolmajor)
{
    if (dim>3u) { fprintf(stderr,"error in apply_win_s: dim must be in [0 3]\n"); return 1; }
    if (dim==0u && R!=L) { fprintf(stderr,"error in apply_win_s: L (winlength) must equal R (nrows X) for dim=0\n"); return 1; }
    if (dim==1u && C!=L) { fprintf(stderr,"error in apply_win_s: L (winlength) must equal C (ncols X) for dim=1\n"); return 1; }
    if (dim==2u && S!=L) { fprintf(stderr,"error in apply_win_s: L (winlength) must equal S (nslices X) for dim=2\n"); return 1; }
    if (dim==3u && H!=L) { fprintf(stderr,"error in apply_win_s: L (winlength) must equal H (nhyperslices X) for dim=3\n"); return 1; }

    const size_t N = R*C*S*H;

    if (N==0u) {}
    else if (L==1u)
    {
        const float x2 = *X2;
        for (size_t n=0u; n<N; ++n, ++X1, ++Y) { *Y = *X1 * x2; }
    }
    else if (N==L)
    {
        for (size_t n=0u; n<N; ++n, ++X1, ++X2, ++Y) { *Y = *X1 * *X2; }
    }
    else if (iscolmajor)
    {
        const int r1i = (int)(R>1u), r2i = (int)(dim==0u);
        const int c1i = (int)R*((int)(C>1u)-(int)(R>1u)), c2i = (int)(dim==1u);
        const int s1i = (int)(R*C)*((int)(S>1u)-(int)(C>1u)), s2i = (int)(dim==2u);
        const int h1i = (int)(R*C*S)*((int)(H>1u)-(int)(S>1u)), h2i = (int)(dim==3u);
        for (size_t h=0u; h<H; ++h, X1+=h1i, X2+=h2i)
        {
            for (size_t s=0u; s<S; ++s, X1+=s1i, X2+=s2i)
            {
                for (size_t c=0u; c<C; ++c, X1+=c1i, X2+=c2i)
                {
                    for (size_t r=0u; r<R; ++r, X1+=r1i, X2+=r2i, ++Y)
                    {
                        *Y = *X1 * *X2;
                    }
                }
            }
        }
    }
    else
    {
        const int h1i = (int)(H>1u), h2i = (int)(dim==3u);
        const int s1i = (int)H*((int)(S>1u)-(int)(H>1u)), s2i = (int)(dim==2u);
        const int c1i = (int)(H*S)*((int)(C>1u)-(int)(S>1u)), c2i = (int)(dim==1u);
        const int r1i = (int)(H*S*C)*((int)(R>1u)-(int)(C>1u)), r2i = (int)(dim==0u);
        for (size_t r=0u; r<R; ++r, X1+=r1i, X2+=r2i)
        {
            for (size_t c=0u; c<C; ++c, X1+=c1i, X2+=c2i)
            {
                for (size_t s=0u; s<S; ++s, X1+=s1i, X2+=s2i)
                {
                    for (size_t h=0u; h<H; ++h, X1+=h1i, X2+=h2i, ++Y)
                    {
                        *Y = *X1 * *X2;
                    }
                }
            }
        }
    }

    return 0;
}


int apply_win_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const size_t dim, const char iscolmajor)
{
    if (dim>3u) { fprintf(stderr,"error in apply_win_d: dim must be in [0 3]\n"); return 1; }
    if (dim==0u && R!=L) { fprintf(stderr,"error in apply_win_d: L (winlength) must equal R (nrows X) for dim=0\n"); return 1; }
    if (dim==1u && C!=L) { fprintf(stderr,"error in apply_win_d: L (winlength) must equal C (ncols X) for dim=1\n"); return 1; }
    if (dim==2u && S!=L) { fprintf(stderr,"error in apply_win_d: L (winlength) must equal S (nslices X) for dim=2\n"); return 1; }
    if (dim==3u && H!=L) { fprintf(stderr,"error in apply_win_d: L (winlength) must equal H (nhyperslices X) for dim=3\n"); return 1; }

    const size_t N = R*C*S*H;

    if (N==0u) {}
    else if (L==1u)
    {
        const float x2 = *X2;
        for (size_t n=0u; n<N; ++n, ++X1, ++Y) { *Y = *X1 * x2; }
    }
    else if (N==L)
    {
        for (size_t n=0u; n<N; ++n, ++X1, ++X2, ++Y) { *Y = *X1 * *X2; }
    }
    else if (iscolmajor)
    {
        const int r1i = (int)(R>1u), r2i = (int)(dim==0u);
        const int c1i = (int)R*((int)(C>1u)-(int)(R>1u)), c2i = (int)(dim==1u);
        const int s1i = (int)(R*C)*((int)(S>1u)-(int)(C>1u)), s2i = (int)(dim==2u);
        const int h1i = (int)(R*C*S)*((int)(H>1u)-(int)(S>1u)), h2i = (int)(dim==3u);
        for (size_t h=0u; h<H; ++h, X1+=h1i, X2+=h2i)
        {
            for (size_t s=0u; s<S; ++s, X1+=s1i, X2+=s2i)
            {
                for (size_t c=0u; c<C; ++c, X1+=c1i, X2+=c2i)
                {
                    for (size_t r=0u; r<R; ++r, X1+=r1i, X2+=r2i, ++Y)
                    {
                        *Y = *X1 * *X2;
                    }
                }
            }
        }
    }
    else
    {
        const int h1i = (int)(H>1u), h2i = (int)(dim==3u);
        const int s1i = (int)H*((int)(S>1u)-(int)(H>1u)), s2i = (int)(dim==2u);
        const int c1i = (int)(H*S)*((int)(C>1u)-(int)(S>1u)), c2i = (int)(dim==1u);
        const int r1i = (int)(H*S*C)*((int)(R>1u)-(int)(C>1u)), r2i = (int)(dim==0u);
        for (size_t r=0u; r<R; ++r, X1+=r1i, X2+=r2i)
        {
            for (size_t c=0u; c<C; ++c, X1+=c1i, X2+=c2i)
            {
                for (size_t s=0u; s<S; ++s, X1+=s1i, X2+=s2i)
                {
                    for (size_t h=0u; h<H; ++h, X1+=h1i, X2+=h2i, ++Y)
                    {
                        *Y = *X1 * *X2;
                    }
                }
            }
        }
    }

    return 0;
}


int apply_win_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const size_t dim, const char iscolmajor)
{
    if (dim>3u) { fprintf(stderr,"error in apply_win_c: dim must be in [0 3]\n"); return 1; }
    if (dim==0u && R!=L) { fprintf(stderr,"error in apply_win_c: L (winlength) must equal R (nrows X) for dim=0\n"); return 1; }
    if (dim==1u && C!=L) { fprintf(stderr,"error in apply_win_c: L (winlength) must equal C (ncols X) for dim=1\n"); return 1; }
    if (dim==2u && S!=L) { fprintf(stderr,"error in apply_win_c: L (winlength) must equal S (nslices X) for dim=2\n"); return 1; }
    if (dim==3u && H!=L) { fprintf(stderr,"error in apply_win_c: L (winlength) must equal H (nhyperslices X) for dim=3\n"); return 1; }

    const size_t N = R*C*S*H;

    if (N==0u) {}
    else if (L==1u)
    {
        const float x2r = *X2, x2i = *(X2+1);
        for (size_t n=0u; n<N; ++n, X1+=2, ++Y)
        {
            *Y = *X1*x2r - *(X1+1)*x2i;
            *++Y = *X1*x2i + *(X1+1)*x2r;
        }
    }
    else if (N==L)
    {
        for (size_t n=0u; n<N; ++n, X1+=2, X2+=2, ++Y)
        {
            *Y = *X1**X2 - *(X1+1)**(X2+1);
            *++Y = *X1**(X2+1) + *(X1+1)**X2;
        }
    }
    else if (iscolmajor)
    {
        const int r1i = 2*(int)(R>1u), r2i = 2*(int)(dim==0u);
        const int c1i = 2*(int)R*((int)(C>1u)-(int)(R>1u)), c2i = 2*(int)(dim==1u);
        const int s1i = 2*(int)(R*C)*((int)(S>1u)-(int)(C>1u)), s2i = 2*(int)(dim==2u);
        const int h1i = 2*(int)(R*C*S)*((int)(H>1u)-(int)(S>1u)), h2i = 2*(int)(dim==3u);
        for (size_t h=0u; h<H; ++h, X1+=h1i, X2+=h2i)
        {
            for (size_t s=0u; s<S; ++s, X1+=s1i, X2+=s2i)
            {
                for (size_t c=0u; c<C; ++c, X1+=c1i, X2+=c2i)
                {
                    for (size_t r=0u; r<R; ++r, X1+=r1i, X2+=r2i, ++Y)
                    {
                        *Y = *X1**X2 - *(X1+1)**(X2+1);
                        *++Y = *X1**(X2+1) + *(X1+1)**X2;
                    }
                }
            }
        }
    }
    else
    {
        const int h1i = (int)(H>1u), h2i = (int)(dim==3u);
        const int s1i = (int)H*((int)(S>1u)-(int)(H>1u)), s2i = (int)(dim==2u);
        const int c1i = (int)(H*S)*((int)(C>1u)-(int)(S>1u)), c2i = (int)(dim==1u);
        const int r1i = (int)(H*S*C)*((int)(R>1u)-(int)(C>1u)), r2i = (int)(dim==0u);
        for (size_t r=0u; r<R; ++r, X1+=r1i, X2+=r2i)
        {
            for (size_t c=0u; c<C; ++c, X1+=c1i, X2+=c2i)
            {
                for (size_t s=0u; s<S; ++s, X1+=s1i, X2+=s2i)
                {
                    for (size_t h=0u; h<H; ++h, X1+=h1i, X2+=h2i, ++Y)
                    {
                        *Y = *X1**X2 - *(X1+1)**(X2+1);
                        *++Y = *X1**(X2+1) + *(X1+1)**X2;
                    }
                }
            }
        }
    }

    return 0;
}


int apply_win_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const size_t dim, const char iscolmajor)
{
    if (dim>3u) { fprintf(stderr,"error in apply_win_z: dim must be in [0 3]\n"); return 1; }
    if (dim==0u && R!=L) { fprintf(stderr,"error in apply_win_z: L (winlength) must equal R (nrows X) for dim=0\n"); return 1; }
    if (dim==1u && C!=L) { fprintf(stderr,"error in apply_win_z: L (winlength) must equal C (ncols X) for dim=1\n"); return 1; }
    if (dim==2u && S!=L) { fprintf(stderr,"error in apply_win_z: L (winlength) must equal S (nslices X) for dim=2\n"); return 1; }
    if (dim==3u && H!=L) { fprintf(stderr,"error in apply_win_z: L (winlength) must equal H (nhyperslices X) for dim=3\n"); return 1; }

    const size_t N = R*C*S*H;

    if (N==0u) {}
    else if (L==1u)
    {
        const float x2r = *X2, x2i = *(X2+1);
        for (size_t n=0u; n<N; ++n, X1+=2, ++Y)
        {
            *Y = *X1*x2r - *(X1+1)*x2i;
            *++Y = *X1*x2i + *(X1+1)*x2r;
        }
    }
    else if (N==L)
    {
        for (size_t n=0u; n<N; ++n, X1+=2, X2+=2, ++Y)
        {
            *Y = *X1**X2 - *(X1+1)**(X2+1);
            *++Y = *X1**(X2+1) + *(X1+1)**X2;
        }
    }
    else if (iscolmajor)
    {
        const int r1i = 2*(int)(R>1u), r2i = 2*(int)(dim==0u);
        const int c1i = 2*(int)R*((int)(C>1u)-(int)(R>1u)), c2i = 2*(int)(dim==1u);
        const int s1i = 2*(int)(R*C)*((int)(S>1u)-(int)(C>1u)), s2i = 2*(int)(dim==2u);
        const int h1i = 2*(int)(R*C*S)*((int)(H>1u)-(int)(S>1u)), h2i = 2*(int)(dim==3u);
        for (size_t h=0u; h<H; ++h, X1+=h1i, X2+=h2i)
        {
            for (size_t s=0u; s<S; ++s, X1+=s1i, X2+=s2i)
            {
                for (size_t c=0u; c<C; ++c, X1+=c1i, X2+=c2i)
                {
                    for (size_t r=0u; r<R; ++r, X1+=r1i, X2+=r2i, ++Y)
                    {
                        *Y = *X1**X2 - *(X1+1)**(X2+1);
                        *++Y = *X1**(X2+1) + *(X1+1)**X2;
                    }
                }
            }
        }
    }
    else
    {
        const int h1i = (int)(H>1u), h2i = (int)(dim==3u);
        const int s1i = (int)H*((int)(S>1u)-(int)(H>1u)), s2i = (int)(dim==2u);
        const int c1i = (int)(H*S)*((int)(C>1u)-(int)(S>1u)), c2i = (int)(dim==1u);
        const int r1i = (int)(H*S*C)*((int)(R>1u)-(int)(C>1u)), r2i = (int)(dim==0u);
        for (size_t r=0u; r<R; ++r, X1+=r1i, X2+=r2i)
        {
            for (size_t c=0u; c<C; ++c, X1+=c1i, X2+=c2i)
            {
                for (size_t s=0u; s<S; ++s, X1+=s1i, X2+=s2i)
                {
                    for (size_t h=0u; h<H; ++h, X1+=h1i, X2+=h2i, ++Y)
                    {
                        *Y = *X1**X2 - *(X1+1)**(X2+1);
                        *++Y = *X1**(X2+1) + *(X1+1)**X2;
                    }
                }
            }
        }
    }

    return 0;
}


int apply_win_inplace_s (float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const size_t dim, const char iscolmajor)
{
    if (dim>3u) { fprintf(stderr,"error in apply_win_inplace_s: dim must be in [0 3]\n"); return 1; }
    if (dim==0u && R!=L) { fprintf(stderr,"error in apply_win_inplace_s: L (winlength) must equal R (nrows X) for dim=0\n"); return 1; }
    if (dim==1u && C!=L) { fprintf(stderr,"error in apply_win_inplace_s: L (winlength) must equal C (ncols X) for dim=1\n"); return 1; }
    if (dim==2u && S!=L) { fprintf(stderr,"error in apply_win_inplace_s: L (winlength) must equal S (nslices X) for dim=2\n"); return 1; }
    if (dim==3u && H!=L) { fprintf(stderr,"error in apply_win_inplace_s: L (winlength) must equal H (nhyperslices X) for dim=3\n"); return 1; }

    const size_t N = R*C*S*H;

    if (N==0u) {}
    else if (L==1u)
    {
        const float x2 = *X2;
        for (size_t n=0u; n<N; ++n, ++X1) { *X1 *= x2; }
    }
    else if (N==L)
    {
        for (size_t n=0u; n<N; ++n, ++X1, ++X2) { *X1 *= *X2; }
    }
    else if (iscolmajor)
    {
        const int r1i = (int)(R>1u), r2i = (int)(dim==0u);
        const int c1i = (int)R*((int)(C>1u)-(int)(R>1u)), c2i = (int)(dim==1u);
        const int s1i = (int)(R*C)*((int)(S>1u)-(int)(C>1u)), s2i = (int)(dim==2u);
        const int h1i = (int)(R*C*S)*((int)(H>1u)-(int)(S>1u)), h2i = (int)(dim==3u);
        for (size_t h=0u; h<H; ++h, X1+=h1i, X2+=h2i)
        {
            for (size_t s=0u; s<S; ++s, X1+=s1i, X2+=s2i)
            {
                for (size_t c=0u; c<C; ++c, X1+=c1i, X2+=c2i)
                {
                    for (size_t r=0u; r<R; ++r, X1+=r1i, X2+=r2i)
                    {
                        *X1 *= *X2;
                    }
                }
            }
        }
    }
    else
    {
        const int h1i = (int)(H>1u), h2i = (int)(dim==3u);
        const int s1i = (int)H*((int)(S>1u)-(int)(H>1u)), s2i = (int)(dim==2u);
        const int c1i = (int)(H*S)*((int)(C>1u)-(int)(S>1u)), c2i = (int)(dim==1u);
        const int r1i = (int)(H*S*C)*((int)(R>1u)-(int)(C>1u)), r2i = (int)(dim==0u);
        for (size_t r=0u; r<R; ++r, X1+=r1i, X2+=r2i)
        {
            for (size_t c=0u; c<C; ++c, X1+=c1i, X2+=c2i)
            {
                for (size_t s=0u; s<S; ++s, X1+=s1i, X2+=s2i)
                {
                    for (size_t h=0u; h<H; ++h, X1+=h1i, X2+=h2i)
                    {
                        *X1 *= *X2;
                    }
                }
            }
        }
    }

    return 0;
}


int apply_win_inplace_d (double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const size_t dim, const char iscolmajor)
{
    if (dim>3u) { fprintf(stderr,"error in apply_win_inplace_d: dim must be in [0 3]\n"); return 1; }
    if (dim==0u && R!=L) { fprintf(stderr,"error in apply_win_inplace_d: L (winlength) must equal R (nrows X) for dim=0\n"); return 1; }
    if (dim==1u && C!=L) { fprintf(stderr,"error in apply_win_inplace_d: L (winlength) must equal C (ncols X) for dim=1\n"); return 1; }
    if (dim==2u && S!=L) { fprintf(stderr,"error in apply_win_inplace_d: L (winlength) must equal S (nslices X) for dim=2\n"); return 1; }
    if (dim==3u && H!=L) { fprintf(stderr,"error in apply_win_inplace_d: L (winlength) must equal H (nhyperslices X) for dim=3\n"); return 1; }

    const size_t N = R*C*S*H;

    if (N==0u) {}
    else if (L==1u)
    {
        const float x2 = *X2;
        for (size_t n=0u; n<N; ++n, ++X1) { *X1 *= x2; }
    }
    else if (N==L)
    {
        for (size_t n=0u; n<N; ++n, ++X1, ++X2) { *X1 *= *X2; }
    }
    else if (iscolmajor)
    {
        const int r1i = (int)(R>1u), r2i = (int)(dim==0u);
        const int c1i = (int)R*((int)(C>1u)-(int)(R>1u)), c2i = (int)(dim==1u);
        const int s1i = (int)(R*C)*((int)(S>1u)-(int)(C>1u)), s2i = (int)(dim==2u);
        const int h1i = (int)(R*C*S)*((int)(H>1u)-(int)(S>1u)), h2i = (int)(dim==3u);
        for (size_t h=0u; h<H; ++h, X1+=h1i, X2+=h2i)
        {
            for (size_t s=0u; s<S; ++s, X1+=s1i, X2+=s2i)
            {
                for (size_t c=0u; c<C; ++c, X1+=c1i, X2+=c2i)
                {
                    for (size_t r=0u; r<R; ++r, X1+=r1i, X2+=r2i)
                    {
                        *X1 *= *X2;
                    }
                }
            }
        }
    }
    else
    {
        const int h1i = (int)(H>1u), h2i = (int)(dim==3u);
        const int s1i = (int)H*((int)(S>1u)-(int)(H>1u)), s2i = (int)(dim==2u);
        const int c1i = (int)(H*S)*((int)(C>1u)-(int)(S>1u)), c2i = (int)(dim==1u);
        const int r1i = (int)(H*S*C)*((int)(R>1u)-(int)(C>1u)), r2i = (int)(dim==0u);
        for (size_t r=0u; r<R; ++r, X1+=r1i, X2+=r2i)
        {
            for (size_t c=0u; c<C; ++c, X1+=c1i, X2+=c2i)
            {
                for (size_t s=0u; s<S; ++s, X1+=s1i, X2+=s2i)
                {
                    for (size_t h=0u; h<H; ++h, X1+=h1i, X2+=h2i)
                    {
                        *X1 *= *X2;
                    }
                }
            }
        }
    }

    return 0;
}


int apply_win_inplace_c (float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const size_t dim, const char iscolmajor)
{
    if (dim>3u) { fprintf(stderr,"error in apply_win_inplace_c: dim must be in [0 3]\n"); return 1; }
    if (dim==0u && R!=L) { fprintf(stderr,"error in apply_win_inplace_c: L (winlength) must equal R (nrows X) for dim=0\n"); return 1; }
    if (dim==1u && C!=L) { fprintf(stderr,"error in apply_win_inplace_c: L (winlength) must equal C (ncols X) for dim=1\n"); return 1; }
    if (dim==2u && S!=L) { fprintf(stderr,"error in apply_win_inplace_c: L (winlength) must equal S (nslices X) for dim=2\n"); return 1; }
    if (dim==3u && H!=L) { fprintf(stderr,"error in apply_win_inplace_c: L (winlength) must equal H (nhyperslices X) for dim=3\n"); return 1; }

    const size_t N = R*C*S*H;
    float xr, xi;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X1)
        {
            xr = *X1**X2 - *(X1+1)**(X2+1);
            xi = *X1**(X2+1) + *(X1+1)**X2;
            *X1 = xr; *++X1 = xi;
        }
    }
    else if (N==L)
    {
        for (size_t n=0u; n<N; ++n, ++X1, X2+=2)
        {
            xr = *X1**X2 - *(X1+1)**(X2+1);
            xi = *X1**(X2+1) + *(X1+1)**X2;
            *X1 = xr; *++X1 = xi;
        }
    }
    else if (iscolmajor)
    {
        const int r1i = 2*(int)(R>1u), r2i = 2*(int)(dim==0u);
        const int c1i = 2*(int)R*((int)(C>1u)-(int)(R>1u)), c2i = 2*(int)(dim==1u);
        const int s1i = 2*(int)(R*C)*((int)(S>1u)-(int)(C>1u)), s2i = 2*(int)(dim==2u);
        const int h1i = 2*(int)(R*C*S)*((int)(H>1u)-(int)(S>1u)), h2i = 2*(int)(dim==3u);
        for (size_t h=0u; h<H; ++h, X1+=h1i, X2+=h2i)
        {
            for (size_t s=0u; s<S; ++s, X1+=s1i, X2+=s2i)
            {
                for (size_t c=0u; c<C; ++c, X1+=c1i, X2+=c2i)
                {
                    for (size_t r=0u; r<R; ++r, X1+=r1i, X2+=r2i)
                    {
                        xr = *X1**X2 - *(X1+1)**(X2+1);
                        xi = *X1**(X2+1) + *(X1+1)**X2;
                        *X1 = xr; *++X1 = xi;
                    }
                }
            }
        }
    }
    else
    {
        const int h1i = 2*(int)(H>1u), h2i = 2*(int)(dim==3u);
        const int s1i = 2*(int)H*((int)(S>1u)-(int)(H>1u)), s2i = 2*(int)(dim==2u);
        const int c1i = 2*(int)(H*S)*((int)(C>1u)-(int)(S>1u)), c2i = 2*(int)(dim==1u);
        const int r1i = 2*(int)(H*S*C)*((int)(R>1u)-(int)(C>1u)), r2i = 2*(int)(dim==0u);
        for (size_t r=0u; r<R; ++r, X1+=r1i, X2+=r2i)
        {
            for (size_t c=0u; c<C; ++c, X1+=c1i, X2+=c2i)
            {
                for (size_t s=0u; s<S; ++s, X1+=s1i, X2+=s2i)
                {
                    for (size_t h=0u; h<H; ++h, X1+=h1i, X2+=h2i)
                    {
                        xr = *X1**X2 - *(X1+1)**(X2+1);
                        xi = *X1**(X2+1) + *(X1+1)**X2;
                        *X1 = xr; *++X1 = xi;
                    }
                }
            }
        }
    }

    return 0;
}


int apply_win_inplace_z (double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const size_t dim, const char iscolmajor)
{
    if (dim>3u) { fprintf(stderr,"error in apply_win_inplace_z: dim must be in [0 3]\n"); return 1; }
    if (dim==0u && R!=L) { fprintf(stderr,"error in apply_win_inplace_z: L (winlength) must equal R (nrows X) for dim=0\n"); return 1; }
    if (dim==1u && C!=L) { fprintf(stderr,"error in apply_win_inplace_z: L (winlength) must equal C (ncols X) for dim=1\n"); return 1; }
    if (dim==2u && S!=L) { fprintf(stderr,"error in apply_win_inplace_z: L (winlength) must equal S (nslices X) for dim=2\n"); return 1; }
    if (dim==3u && H!=L) { fprintf(stderr,"error in apply_win_inplace_z: L (winlength) must equal H (nhyperslices X) for dim=3\n"); return 1; }

    const size_t N = R*C*S*H;
    double xr, xi;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X1)
        {
            xr = *X1**X2 - *(X1+1)**(X2+1);
            xi = *X1**(X2+1) + *(X1+1)**X2;
            *X1 = xr; *++X1 = xi;
        }
    }
    else if (N==L)
    {
        for (size_t n=0u; n<N; ++n, ++X1, X2+=2)
        {
            xr = *X1**X2 - *(X1+1)**(X2+1);
            xi = *X1**(X2+1) + *(X1+1)**X2;
            *X1 = xr; *++X1 = xi;
        }
    }
    else if (iscolmajor)
    {
        const int r1i = 2*(int)(R>1u), r2i = 2*(int)(dim==0u);
        const int c1i = 2*(int)R*((int)(C>1u)-(int)(R>1u)), c2i = 2*(int)(dim==1u);
        const int s1i = 2*(int)(R*C)*((int)(S>1u)-(int)(C>1u)), s2i = 2*(int)(dim==2u);
        const int h1i = 2*(int)(R*C*S)*((int)(H>1u)-(int)(S>1u)), h2i = 2*(int)(dim==3u);
        for (size_t h=0u; h<H; ++h, X1+=h1i, X2+=h2i)
        {
            for (size_t s=0u; s<S; ++s, X1+=s1i, X2+=s2i)
            {
                for (size_t c=0u; c<C; ++c, X1+=c1i, X2+=c2i)
                {
                    for (size_t r=0u; r<R; ++r, X1+=r1i, X2+=r2i)
                    {
                        xr = *X1**X2 - *(X1+1)**(X2+1);
                        xi = *X1**(X2+1) + *(X1+1)**X2;
                        *X1 = xr; *++X1 = xi;
                    }
                }
            }
        }
    }
    else
    {
        const int h1i = 2*(int)(H>1u), h2i = 2*(int)(dim==3u);
        const int s1i = 2*(int)H*((int)(S>1u)-(int)(H>1u)), s2i = 2*(int)(dim==2u);
        const int c1i = 2*(int)(H*S)*((int)(C>1u)-(int)(S>1u)), c2i = 2*(int)(dim==1u);
        const int r1i = 2*(int)(H*S*C)*((int)(R>1u)-(int)(C>1u)), r2i = 2*(int)(dim==0u);
        for (size_t r=0u; r<R; ++r, X1+=r1i, X2+=r2i)
        {
            for (size_t c=0u; c<C; ++c, X1+=c1i, X2+=c2i)
            {
                for (size_t s=0u; s<S; ++s, X1+=s1i, X2+=s2i)
                {
                    for (size_t h=0u; h<H; ++h, X1+=h1i, X2+=h2i)
                    {
                        xr = *X1**X2 - *(X1+1)**(X2+1);
                        xi = *X1**(X2+1) + *(X1+1)**X2;
                        *X1 = xr; *++X1 = xi;
                    }
                }
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
