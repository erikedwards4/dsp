//This gets LCs as usual, and then applies window W to get rate.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int lcr_windowed_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const float *W, const int L, const int dim, const int c0, const float stp, const int going, const float level);
int lcr_windowed_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const double *W, const int L, const int dim, const int c0, const double stp, const int going, const double level);


int lcr_windowed_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const float *W, const int L, const int dim, const int c0, const float stp, const int going, const float level)
{
    const float z = 0.0f;
    const int N = R*C;
    const int T = (dim==0) ? 1 + (int)(floorf(((int)R-1-c0)/stp)) : 1 + (int)(floorf(((int)C-1-c0)/stp));
    const int Lpre = L/2;          //nsamps before center samp
    const int Lpost = L - L/2 - 1; //nsamps after center samp
    int t = 0, cs = c0;
    int c, r, l, n = -1;
    char *Z;
    float *W2;

    //Checks
    if (R<1) { fprintf(stderr,"error in lcr_windowed_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in lcr_windowed_s: C (ncols X) must be positive\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in lcr_windowed_s: L (winlength) must be positive\n"); return 1; }

    //Initialize Z
    if (!(Z=(char *)malloc((size_t)(N)*sizeof(char)))) { fprintf(stderr,"error in lcr_windowed_s: problem with malloc. "); perror("malloc"); return 1; }

    //Initialize W2 (same as W if L even, else 0.5*W added to itself with lag 1)
    if (!(W2=(float *)malloc((size_t)(L)*sizeof(float)))) { fprintf(stderr,"error in lcr_windowed_s: problem with malloc. "); perror("malloc"); return 1; }
    if (L%2) { W2[0] = 0.5f*W[0]; for (l=1; l<L; l++) { W2[l] = 0.5f*(W[l]+W[l-1]); } }
    else { for (l=0; l<L; l++) { W2[l] = W[l]; } }

    //Raw material
    if (going>0) { while (++n<N) { Z[n] = (X[n]>=level); } }
    else { while (++n<N) { Z[n] = (X[n]<level); } }

    if (dim==0)
    {
        cblas_scopy(T*C,&z,0,&Y[0],1);
        if (iscolmajor)
        {
            if (going) { while (--n>0) { Z[n] *= (Z[n]!=Z[n-1]); } }
            else { while (--n>0) { Z[n] = (Z[n]!=Z[n-1]); } }
            while (n<N) { Z[n] = 0; n += R; }
            while (t<T && cs<=Lpre)
            {
                for (c=0; c<C; c++) { for (l=Lpre-cs; l<L; l++) { if (Z[c*R+cs+l-Lpre]) { Y[t+c*T] += W2[l]; } } }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
            while (t<T && cs<R-Lpost)
            {
                for (c=0; c<C; c++) { for (l=0; l<L; l++) { if (Z[c*R+cs+l-Lpre]) { Y[t+c*T] += W2[l]; } } }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
            while (t<T)
            {
                for (c=0; c<C; c++) { for (l=0; l<Lpre+R-cs; l++) { if (Z[c*R+cs+l-Lpre]) { Y[t+c*T] += W2[l]; } } }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
        }
        else
        {
            if (going) { while (--n>=C) { Z[n] *= (Z[n]!=Z[n-C]); } }
            else { while (--n>=C) { Z[n] = (Z[n]!=Z[n-C]); } }
            while (n>=0) { Z[n--] = 0; }
            while (t<T && cs<=Lpre)
            {
                for (c=0; c<C; c++) { for (l=Lpre-cs; l<L; l++) { if (Z[c+(cs+l-Lpre)*C]) { Y[c+t*C] += W2[l]; } } }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
            while (t<T && cs<R-Lpost)
            {
                for (c=0; c<C; c++) { for (l=0; l<L; l++) { if (Z[c+(cs+l-Lpre)*C]) { Y[c+t*C] += W2[l]; } } }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
            while (t<T)
            {
                for (c=0; c<C; c++) { for (l=0; l<Lpre+R-cs; l++) { if (Z[c+(cs+l-Lpre)*C]) { Y[c+t*C] += W2[l]; } } }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
        }
    }
    else if (dim==1)
    {
        cblas_scopy(R*T,&z,0,&Y[0],1);
        if (iscolmajor)
        {
            if (going) { while (--n>=R) { Z[n] *= (Z[n]!=Z[n-R]); } }
            else { while (--n>=R) { Z[n] = (Z[n]!=Z[n-R]); } }
            while (n>=0) { Z[n--] = 0; }
            while (t<T && cs<=Lpre)
            {
                for (r=0; r<R; r++) { for (l=Lpre-cs; l<L; l++) { if (Z[r+(cs+l-Lpre)*R]) { Y[r+t*R] += W2[l]; } } }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
            while (t<T && cs<R-Lpost)
            {
                for (r=0; r<R; r++) { for (l=0; l<L; l++) { if (Z[r+(cs+l-Lpre)*R]) { Y[r+t*R] += W2[l]; } } }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
            while (t<T)
            {
                for (r=0; r<R; r++) { for (l=0; l<Lpre+C-cs; l++) { if (Z[r+(cs+l-Lpre)*R]) { Y[r+t*R] += W2[l]; } } }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
        }
        else
        {
            if (going) { while (--n>0) { Z[n] *= (Z[n]!=Z[n-1]); } }
            else { while (--n>0) { Z[n] = (Z[n]!=Z[n-1]); } }
            while (n<N) { Z[n] = 0; n += C; }
            while (t<T && cs<=Lpre)
            {
                for (r=0; r<R; r++) { for (l=Lpre-cs; l<L; l++) { if (Z[r*C+cs+l-Lpre]) { Y[t+r*T] += W2[l]; } } }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
            while (t<T && cs<R-Lpost)
            {
                for (r=0; r<R; r++) { for (l=0; l<L; l++) { if (Z[r*C+cs+l-Lpre]) { Y[t+r*T] += W2[l]; } } }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
            while (t<T)
            {
                for (r=0; r<R; r++) { for (l=0; l<Lpre+C-cs; l++) { if (Z[r*C+cs+l-Lpre]) { Y[t+r*T] += W2[l]; } } }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
        }
    }
    else
    {
        fprintf(stderr,"error in lcr_windowed_s: dim must be 0 or 1.\n"); return 1;
    }

    free(Z); free(W2);
    return 0;
}


int lcr_windowed_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const double *W, const int L, const int dim, const int c0, const double stp, const int going, const double level)
{
    const double z = 0.0;
    const int N = R*C;
    const int T = (dim==0) ? 1 + (int)(floor(((int)R-1-c0)/stp)) : 1 + (int)(floor(((int)C-1-c0)/stp));
    const int Lpre = L/2;          //nsamps before center samp
    const int Lpost = L - L/2 - 1; //nsamps after center samp
    int cs = c0, t = 0;
    int c, r, l, n = -1;
    char *Z;
    double *W2;

    //Checks
    if (R<1) { fprintf(stderr,"error in lcr_windowed_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in lcr_windowed_d: C (ncols X) must be positive\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in lcr_windowed_d: L (winlength) must be positive\n"); return 1; }

    //Initialize Z
    if (!(Z=(char *)malloc((size_t)(N)*sizeof(char)))) { fprintf(stderr,"error in lcr_windowed_d: problem with malloc. "); perror("malloc"); return 1; }

    //Initialize W2 (same as W if L even, else 0.5*W added to itself with lag 1)
    if (!(W2=(double *)malloc((size_t)(L)*sizeof(double)))) { fprintf(stderr,"error in lcr_windowed_d: problem with malloc. "); perror("malloc"); return 1; }
    if (L%2) { W2[0] = 0.5*W[0]; for (l=1; l<L; l++) { W2[l] = 0.5*(W[l]+W[l-1]); } }
    else { for (l=0; l<L; l++) { W2[l] = W[l]; } }

    //Raw material
    if (going>0) { while (++n<N) { Z[n] = (X[n]>=level); } }
    else { while (++n<N) { Z[n] = (X[n]<level); } }

    if (dim==0)
    {
        cblas_dcopy(T*C,&z,0,&Y[0],1);
        if (iscolmajor)
        {
            if (going) { while (--n>0) { Z[n] *= (Z[n]!=Z[n-1]); } }
            else { while (--n>0) { Z[n] = (Z[n]!=Z[n-1]); } }
            while (n<N) { Z[n] = 0; n += R; }
            while (t<T && cs<=Lpre)
            {
                for (c=0; c<C; c++) { for (l=Lpre-cs; l<L; l++) { if (Z[c*R+cs+l-Lpre]) { Y[t+c*T] += W2[l]; } } }
                t++; cs = (int)(round(t*stp)) + c0;
            }
            while (t<T && cs<R-Lpost)
            {
                for (c=0; c<C; c++) { for (l=0; l<L; l++) { if (Z[c*R+cs+l-Lpre]) { Y[t+c*T] += W2[l]; } } }
                t++; cs = (int)(round(t*stp)) + c0;
            }
            while (t<T)
            {
                for (c=0; c<C; c++) { for (l=0; l<Lpre+R-cs; l++) { if (Z[c*R+cs+l-Lpre]) { Y[t+c*T] += W2[l]; } } }
                t++; cs = (int)(round(t*stp)) + c0;
            }
        }
        else
        {
            if (going) { while (--n>=C) { Z[n] *= (Z[n]!=Z[n-C]); } }
            else { while (--n>=C) { Z[n] = (Z[n]!=Z[n-C]); } }
            while (n>=0) { Z[n--] = 0; }
            while (t<T && cs<=Lpre)
            {
                for (c=0; c<C; c++) { for (l=Lpre-cs; l<L; l++) { if (Z[c+(cs+l-Lpre)*C]) { Y[c+t*C] += W2[l]; } } }
                t++; cs = (int)(round(t*stp)) + c0;
            }
            while (t<T && cs<R-Lpost)
            {
                for (c=0; c<C; c++) { for (l=0; l<L; l++) { if (Z[c+(cs+l-Lpre)*C]) { Y[c+t*C] += W2[l]; } } }
                t++; cs = (int)(round(t*stp)) + c0;
            }
            while (t<T)
            {
                for (c=0; c<C; c++) { for (l=0; l<Lpre+R-cs; l++) { if (Z[c+(cs+l-Lpre)*C]) { Y[c+t*C] += W2[l]; } } }
                t++; cs = (int)(round(t*stp)) + c0;
            }
        }
    }
    else if (dim==1)
    {
        cblas_dcopy(R*T,&z,0,&Y[0],1);
        if (iscolmajor)
        {
            if (going) { while (--n>=R) { Z[n] *= (Z[n]!=Z[n-R]); } }
            else { while (--n>=R) { Z[n] = (Z[n]!=Z[n-R]); } }
            while (n>=0) { Z[n--] = 0; }
            while (t<T && cs<=Lpre)
            {
                for (r=0; r<R; r++) { for (l=Lpre-cs; l<L; l++) { if (Z[r+(cs+l-Lpre)*R]) { Y[r+t*R] += W2[l]; } } }
                t++; cs = (int)(round(t*stp)) + c0;
            }
            while (t<T && cs<R-Lpost)
            {
                for (r=0; r<R; r++) { for (l=0; l<L; l++) { if (Z[r+(cs+l-Lpre)*R]) { Y[r+t*R] += W2[l]; } } }
                t++; cs = (int)(round(t*stp)) + c0;
            }
            while (t<T)
            {
                for (r=0; r<R; r++) { for (l=0; l<Lpre+C-cs; l++) { if (Z[r+(cs+l-Lpre)*R]) { Y[r+t*R] += W2[l]; } } }
                t++; cs = (int)(round(t*stp)) + c0;
            }
        }
        else
        {
            if (going) { while (--n>0) { Z[n] *= (Z[n]!=Z[n-1]); } }
            else { while (--n>0) { Z[n] = (Z[n]!=Z[n-1]); } }
            while (n<N) { Z[n] = 0; n += C; }
            while (t<T && cs<=Lpre)
            {
                for (r=0; r<R; r++) { for (l=Lpre-cs; l<L; l++) { if (Z[r*C+cs+l-Lpre]) { Y[t+r*T] += W2[l]; } } }
                t++; cs = (int)(round(t*stp)) + c0;
            }
            while (t<T && cs<R-Lpost)
            {
                for (r=0; r<R; r++) { for (l=0; l<L; l++) { if (Z[r*C+cs+l-Lpre]) { Y[t+r*T] += W2[l]; } } }
                t++; cs = (int)(round(t*stp)) + c0;
            }
            while (t<T)
            {
                for (r=0; r<R; r++) { for (l=0; l<Lpre+C-cs; l++) { if (Z[r*C+cs+l-Lpre]) { Y[t+r*T] += W2[l]; } } }
                t++; cs = (int)(round(t*stp)) + c0;
            }
        }
    }
    else
    {
        fprintf(stderr,"error in lcr_windowed_d: dim must be 0 or 1.\n"); return 1;
    }

    free(Z); free(W2);
    return 0;
}


#ifdef __cplusplus
}
}
#endif

