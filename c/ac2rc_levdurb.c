//Gets reflection coefficients (RCs) from the autocorrelation (AC) function for each row or col of X.
//Starts with a Levinson-Durbin recursion from the AC values.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ac2rc_levdurb_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2rc_levdurb_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2rc_levdurb_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2rc_levdurb_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);


int ac2rc_levdurb_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ac2rc_levdurb_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else
    {
        if (L==N)
        {
        }
        else
        {}
    }

    const int P = (dim==0) ? R-1 : C-1;
    int r, c, p, q;
	float g, v;
    float *AS, *AStmp;

    //Checks
    if (R<1) { fprintf(stderr,"error in ac2rcs_s: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ac2rcs_s: ncols X must be positive\n"); return 1; }
    if (P<1) { fprintf(stderr,"error in ac2rcs_s: P (length of poly coeffs including a0) must be positive\n"); return 1; }
	
    //Initialize
    if (!(AS=(float *)malloc((size_t)P*sizeof(float)))) { fprintf(stderr,"error in ac2rcs_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(AStmp=(float *)malloc((size_t)(P-1)*sizeof(float)))) { fprintf(stderr,"error in ac2rcs_s: problem with malloc. "); perror("malloc"); return 1; }
	AS[0] = 1.0f;

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                Y[c*R] = AS[1] = AStmp[0] = g = -X[1+c*R]/X[c*R];
                v = fmaf(X[1+c*R],g,X[c*R]);
                for (p=2; p<R; p++)
                {
                    g = X[p+c*R];
                    for (q=1; q<p; q++) { g = fmaf(X[q+c*R],AS[p-q],g); }
                    Y[p-1+c*R] = AS[p] = g = -g/v;
                    for (q=1; q<p; q++) { AS[q] = fmaf(g,AStmp[p-q-1],AS[q]); }
                    cblas_scopy(p,&AS[1],1,&AStmp[0],1);
                    v *= fmaf(g,-g,1.0f);
                }
                g = X[p+c*R];
                for (q=1; q<p; q++) { g = fmaf(X[q+c*R],AS[p-q],g); }
                Y[p-1+c*R] = -g/v;
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                Y[c] = AS[1] = AStmp[0] = g = -X[c+C]/X[c];
                v = fmaf(X[c+C],g,X[c]);
                for (p=2; p<P; p++)
                {
                    g = X[c+p*C];
                    for (q=1; q<p; q++) { g = fmaf(X[c+q*C],AS[p-q],g); }
                    Y[c+(p-1)*C] = AS[p] = g = -g/v;
                    for (q=1; q<p; q++) { AS[q] = fmaf(g,AStmp[p-q-1],AS[q]); }
                    cblas_scopy(p,&AS[1],C,&AStmp[0],1);
                    v *= fmaf(g,-g,1.0f);
                }
                g = X[c+p*C];
                for (q=1; q<p; q++) { g = fmaf(X[c+q*C],AS[p-q],g); }
                Y[c+(p-1)*C] = -g/v;
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                Y[r] = AS[1] = AStmp[0] = g = -X[r+R]/X[r];
                v = fmaf(X[r+R],g,X[r]);
                for (p=2; p<P; p++)
                {
                    g = X[r+p*R];
                    for (q=1; q<p; q++) { g = fmaf(X[r+q*R],AS[p-q],g); }
                    Y[r+(p-1)*R] = AS[p] = g = -g/v;
                    for (q=1; q<p; q++) { AS[q] = fmaf(g,AStmp[p-q-1],AS[q]); }
                    cblas_scopy(p,&AS[1],R,&AStmp[0],1);
                    v *= fmaf(g,-g,1.0f);
                }
                g = X[r+p*R];
                for (q=1; q<p; q++) { g = fmaf(X[r+q*R],AS[p-q],g); }
                Y[r+(p-1)*R] = -g/v;
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                Y[r*C] = AS[1] = AStmp[0] = g = -X[1+r*C]/X[r*C];
                v = fmaf(X[1+r*C],g,X[r*C]);
                for (p=2; p<P; p++)
                {
                    g = X[p+r*C];
                    for (q=1; q<p; q++) { g = fmaf(X[q+r*C],AS[p-q],g); }
                    Y[p-1+r*C] = AS[p] = g = -g/v;
                    for (q=1; q<p; q++) { AS[q] = fmaf(g,AStmp[p-q-1],AS[q]); }
                    cblas_scopy(p,&AS[1],1,&AStmp[0],1);
                    v *= fmaf(g,-g,1.0f);
                }
                g = X[p+r*C];
                for (q=1; q<p; q++) { g = fmaf(X[q+r*C],AS[p-q],g); }
                Y[p-1+r*C] = -g/v;
            }
        }
    }
    else
    {
        fprintf(stderr,"error in ac2rcs_s: dim must be 0 or 1.\n"); return 1;
    }
	
	return 0;
}


int ac2rcs_d (double *Y, const double *X, const int iscolmajor, const int R, const int C, const int dim)
{
	const int P = (dim==0) ? R-1 : C-1;
    int r, c, p, q;
	double g, v;
    double *AS, *AStmp;

    //Checks
    if (R<1) { fprintf(stderr,"error in ac2rcs_d: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ac2rcs_d: ncols X must be positive\n"); return 1; }
    if (P<1) { fprintf(stderr,"error in ac2rcs_d: P (length of poly coeffs including a0) must be positive\n"); return 1; }
	
    //Initialize
    if (!(AS=(double *)malloc((size_t)P*sizeof(double)))) { fprintf(stderr,"error in ac2rcs_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(AStmp=(double *)malloc((size_t)(P-1)*sizeof(double)))) { fprintf(stderr,"error in ac2rcs_d: problem with malloc. "); perror("malloc"); return 1; }
    AS[0] = 1.0;

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                Y[c*R] = AS[1] = AStmp[0] = g = -X[1+c*R]/X[c*R];
                v = fma(X[1+c*R],g,X[c*R]);
                for (p=2; p<P; p++)
                {
                    g = X[p+c*R];
                    for (q=1; q<p; q++) { g = fma(X[q+c*R],AS[p-q],g); }
                    Y[p-1+c*R] = AS[p] = g = -g/v;
                    for (q=1; q<p; q++) { AS[q] = fma(g,AStmp[p-q-1],AS[q]); }
                    cblas_dcopy(p,&AS[1],1,&AStmp[0],1);
                    v *= fma(g,-g,1.0);
                }
                g = X[p+c*R];
                for (q=1; q<p; q++) { g = fma(X[q+c*R],AS[p-q],g); }
                Y[p-1+c*R] = -g/v;
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                Y[c] = AS[1] = AStmp[0] = g = -X[c+C]/X[c];
                v = fma(X[c+C],g,X[c]);
                for (p=2; p<P; p++)
                {
                    g = X[c+p*C];
                    for (q=1; q<p; q++) { g = fma(X[c+q*C],AS[p-q],g); }
                    Y[c+(p-1)*C] = AS[p] = g = -g/v;
                    for (q=1; q<p; q++) { AS[q] = fma(g,AStmp[p-q-1],AS[q]); }
                    cblas_dcopy(p,&AS[1],C,&AStmp[0],1);
                    v *= fma(g,-g,1.0);
                }
                g = X[c+p*C];
                for (q=1; q<p; q++) { g = fma(X[c+q*C],AS[p-q],g); }
                Y[c+(p-1)*C] = -g/v;
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                Y[r] = AS[1] = AStmp[0] = g = -X[r+R]/X[r];
                v = fma(X[r+R],g,X[r]);
                for (p=2; p<P; p++)
                {
                    g = X[r+p*R];
                    for (q=1; q<p; q++) { g = fma(X[r+q*R],AS[p-q],g); }
                    Y[r+(p-1)*R] = AS[p] = g = -g/v;
                    for (q=1; q<p; q++) { AS[q] = fma(g,AStmp[p-q-1],AS[q]); }
                    cblas_dcopy(p,&AS[1],R,&AStmp[0],1);
                    v *= fma(g,-g,1.0);
                }
                g = X[r+p*R];
                for (q=1; q<p; q++) { g = fma(X[r+q*R],AS[p-q],g); }
                Y[r+(p-1)*R] = -g/v;
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                Y[r*C] = AS[1] = AStmp[0] = g = -X[1+r*C]/X[r*C];
                v = fma(X[1+r*C],g,X[r*C]);
                for (p=2; p<P; p++)
                {
                    g = X[p+r*C];
                    for (q=1; q<p; q++) { g = fma(X[q+r*C],AS[p-q],g); }
                    Y[p-1+r*C] = AS[p] = g = -g/v;
                    for (q=1; q<p; q++) { AS[q] = fma(g,AStmp[p-q-1],AS[q]); }
                    cblas_dcopy(p,&AS[1],1,&AStmp[0],1);
                    v *= fma(g,-g,1.0);
                }
                g = X[p+r*C];
                for (q=1; q<p; q++) { g = fma(X[q+r*C],AS[p-q],g); }
                Y[p-1+r*C] = -g/v;
            }
        }
    }
    else
    {
        fprintf(stderr,"error in ac2rcs_d: dim must be 0 or 1.\n"); return 1;
    }
	
	return 0;
}


#ifdef __cplusplus
}
#endif

