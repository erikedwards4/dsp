//Polynomial from roots (like Matlab/Octave with vector input)

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int poly_s (float *AS, const float *roots, const int R);
int poly_d (double *AS, const double *roots, const int R);
int poly_c (float *AS, const float *roots, const int R);
int poly_z (double *AS, const double *roots, const int R);


int poly_s (float *AS, const float *roots, const int R)
{
    int r, p;

    //Checks
    if (R<1) { fprintf(stderr,"error in poly_s: R (length of roots) must be positive\n"); return 1; }

    AS[0] = 1.0f;
    for (r=0; r<R; r++)
    {
        AS[r+1] = 0.0f;
        for (p=r+1; p>0; p--) { AS[p] = fmaf(-roots[r],AS[p-1],AS[p]); }
        //for (p=r+1; p>0; p--) { AS[p] -= roots[r] * AS[p-1]; }
    }

    return 0;
}


int poly_d (double *AS, const double *roots, const int R)
{
    int r, p;

    //Checks
    if (R<1) { fprintf(stderr,"error in poly_d: R (length of roots) must be positive\n"); return 1; }

    AS[0] = 1.0;
    for (r=0; r<R; r++)
    {
        AS[r+1] = 0.0;
        for (p=r+1; p>0; p--) { AS[p] = fma(-roots[r],AS[p-1],AS[p]); }
        //for (p=r+1; p>0; p--) { AS[p] -= roots[r] * AS[p-1]; }
    }

    return 0;
}


int poly_c (float *AS, const float *roots, const int R)
{
    int r, p;

    //Checks
    if (R<1) { fprintf(stderr,"error in poly_c: R (length of roots) must be positive\n"); return 1; }

    AS[0] = 1.0f; AS[1] = 0.0f;
    for (r=0; r<R; r++)
    {
        AS[2*r+2] = AS[2*r+3] = 0.0f;
        for (p=r+1; p>0; p--)
        {
            AS[2*p] = fmaf(-roots[2*r],AS[2*p-2],fmaf(roots[2*r+1],AS[2*p-1],AS[2*p]));
            AS[2*p+1] = fmaf(-roots[2*r],AS[2*p-1],fmaf(-roots[2*r+1],AS[2*p-2],AS[2*p+1]));
            //AS[2*p] -= roots[2*r]*AS[2*p-2] - roots[2*r+1]*AS[2*p-1];
            //AS[2*p+1] -= roots[2*r]*AS[2*p-1] + roots[2*r+1]*AS[2*p-2];
        }
    }

    return 0;
}


int poly_z (double *AS, const double *roots, const int R)
{
    int r, p;

    //Checks
    if (R<1) { fprintf(stderr,"error in poly_z: R (length of roots) must be positive\n"); return 1; }

    AS[0] = 1.0; AS[1] = 0.0;
    for (r=0; r<R; r++)
    {
        AS[2*r+2] = AS[2*r+3] = 0.0;
        for (p=r+1; p>0; p--)
        {
            AS[2*p] = fma(-roots[2*r],AS[2*p-2],fma(roots[2*r+1],AS[2*p-1],AS[2*p]));
            AS[2*p+1] = fma(-roots[2*r],AS[2*p-1],fma(-roots[2*r+1],AS[2*p-2],AS[2*p+1]));
            //AS[2*p] -= roots[2*r]*AS[2*p-2] - roots[2*r+1]*AS[2*p-1];
            //AS[2*p+1] -= roots[2*r]*AS[2*p-1] + roots[2*r+1]*AS[2*p-2];
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif

