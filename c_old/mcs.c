//Mean-crossings

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int mcs_s (char *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int going);
int mcs_d (char *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int going);


int mcs_s (char *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int going)
{
    const float o = 1.0f;
    const int N = R*C;
    int n = -1, ns;
    float m;

    //Checks
    if (R<1) { fprintf(stderr,"error in mcs_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in mcs_s: C (ncols X) must be positive\n"); return 1; }

    //Get mean
    if (dim==0)
    {
        if (iscolmajor)
        {
            for (ns=0; ns<N; ns+=R)
            {
                m = cblas_sdot(R,&X[ns],1,&o,0) / R;
                n = ns - 1;
                if (going>0) { while (++n<ns+R) { Y[n] = (X[n]>=m); } }
                else { while (++n<ns+R) { Y[n] = (X[n]<m); } }
                if (going) { while (--n>ns) { Y[n] *= (Y[n]!=Y[n-1]); } }
                else { while (--n>ns) { Y[n] = (Y[n]!=Y[n-1]); } }
            }
            while (n<N) { Y[n] = 0; n += R; }
        }
        else
        {
            for (ns=0; ns<C; ns++)
            {
                m = cblas_sdot(R,&X[ns],C,&o,0) / R;
                n = ns;
                if (going>0) { while (n<N) { Y[n] = (X[n]>=m); n += C; } }
                else { while (n<N) { Y[n] = (X[n]<m); n += C; } }
                n -= C;
                if (going) { while (n>ns) { Y[n] *= (Y[n]!=Y[n-C]); n -= C; } }
                else { while (n>ns) { Y[n] = (Y[n]!=Y[n-C]); n -= C; } }
            }
            while (n>=0) { Y[n--] = 0; }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (ns=0; ns<R; ns++)
            {
                m = cblas_sdot(C,&X[ns],R,&o,0) / C;
                n = ns;
                if (going>0) { while (n<N) { Y[n] = (X[n]>=m); n += R; } }
                else { while (n<N) { Y[n] = (X[n]<m); n += R; } }
                n -= R;
                if (going) { while (n>ns) { Y[n] *= (Y[n]!=Y[n-R]); n -= R; } }
                else { while (n>ns) { Y[n] = (Y[n]!=Y[n-R]); n -= R; }  }
            }
            while (n>=0) { Y[n--] = 0; }
        }
        else
        {
            for (ns=0; ns<N; ns+=C)
            {
                m = cblas_sdot(C,&X[ns],1,&o,0) / C;
                n = ns - 1;
                if (going>0) { while (++n<ns+C) { Y[n] = (X[n]>=m); } }
                else { while (++n<ns+C) { Y[n] = (X[n]<m); } }
                if (going) { while (--n>ns) { Y[n] *= (Y[n]!=Y[n-1]); } }
                else { while (--n>ns) { Y[n] = (Y[n]!=Y[n-1]); } }
            }
            while (n<N) { Y[n] = 0; n += C; }
        }
    }
    else
    {
        fprintf(stderr,"error in mcs_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int mcs_d (char *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int going)
{
    const double o = 1.0;
    const int N = R*C;
    int n = -1, ns;
    double m;

    //Checks
    if (R<1) { fprintf(stderr,"error in mcs_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in mcs_d: C (ncols X) must be positive\n"); return 1; }

    //Get mean
    if (dim==0)
    {
        if (iscolmajor)
        {
            for (ns=0; ns<N; ns+=R)
            {
                m = cblas_ddot(R,&X[ns],1,&o,0) / R;
                n = ns - 1;
                if (going>0) { while (++n<ns+R) { Y[n] = (X[n]>=m); } }
                else { while (++n<ns+R) { Y[n] = (X[n]<m); } }
                if (going) { while (--n>ns) { Y[n] *= (Y[n]!=Y[n-1]); } }
                else { while (--n>ns) { Y[n] = (Y[n]!=Y[n-1]); } }
            }
            while (n<N) { Y[n] = 0; n += R; }
        }
        else
        {
            for (ns=0; ns<C; ns++)
            {
                m = cblas_ddot(R,&X[ns],C,&o,0) / R;
                n = ns;
                if (going>0) { while (n<N) { Y[n] = (X[n]>=m); n += C; } }
                else { while (n<N) { Y[n] = (X[n]<m); n += C; } }
                n -= C;
                if (going) { while (n>ns) { Y[n] *= (Y[n]!=Y[n-C]); n -= C; } }
                else { while (n>ns) { Y[n] = (Y[n]!=Y[n-C]); n -= C; } }
            }
            while (n>=0) { Y[n--] = 0; }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (ns=0; ns<R; ns++)
            {
                m = cblas_ddot(C,&X[ns],R,&o,0) / C;
                n = ns;
                if (going>0) { while (n<N) { Y[n] = (X[n]>=m); n += R; } }
                else { while (n<N) { Y[n] = (X[n]<m); n += R; } }
                n -= R;
                if (going) { while (n>ns) { Y[n] *= (Y[n]!=Y[n-R]); n -= R; } }
                else { while (n>ns) { Y[n] = (Y[n]!=Y[n-R]); n -= R; } }
            }
            while (n>=0) { Y[n--] = 0; }
        }
        else
        {
            for (ns=0; ns<N; ns+=C)
            {
                m = cblas_ddot(C,&X[ns],1,&o,0) / C;
                n = ns - 1;
                if (going>0) { while (++n<ns+C) { Y[n] = (X[n]>=m); } }
                else { while (++n<ns+C) { Y[n] = (X[n]<m); } }
                if (going) { while (--n>ns) { Y[n] *= (Y[n]!=Y[n-1]); } }
                else { while (--n>ns) { Y[n] = (Y[n]!=Y[n-1]); } }
            }
            while (n<N) { Y[n] = 0; n += C; }
        }
    }
    else
    {
        fprintf(stderr,"error in mcs_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif

