//This gets MCs as usual, and then sums over (rectangular) window of length L.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int mcr_s (float *Y, const float *X, const int iscolmajor, const int R, const int C, const int L, const int dim, const int c0, const float stp, const int going);
int mcr_d (double *Y, const double *X, const int iscolmajor, const int R, const int C, const int L, const int dim, const int c0, const double stp, const int going);


int mcr_s (float *Y, const float *X, const int iscolmajor, const int R, const int C, const int L, const int dim, const int c0, const float stp, const int going)
{
    const float o = 1.0f;
    const int N = R*C;
    const int T = (dim==0) ? 1 + (int)(floorf(((int)R-1-c0)/stp)) : 1 + (int)(floorf(((int)C-1-c0)/stp));
    const int Lpre = L/2;          //nsamps before center samp
    const int Lpost = L - L/2 - 1; //nsamps after center samp
    int cs = c0, t = 0;
    int c, r, ns, n = -1;
    int *Z;
    float sc, m;

    //Checks
    if (R<1) { fprintf(stderr,"error in mcr_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in mcr_s: C (ncols X) must be positive\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in mcr_s: L (winlength) must be positive\n"); return 1; }

    //Initialize Z
    if (!(Z=(int *)malloc((size_t)(N)*sizeof(int)))) { fprintf(stderr,"error in mcr_s: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (ns=0; ns<N; ns+=R)
            {
                m = cblas_sdot(R,&X[ns],1,&o,0) / R;
                n = ns - 1;
                if (going>0) { while (++n<ns+R) { Z[n] = (X[n]>=m); } }
                else { while (++n<ns+R) { Z[n] = (X[n]<m); } }
                if (going) { while (--n>ns) { Z[n] *= (Z[n]!=Z[n-1]); } }
                else { while (--n>ns) { Z[n] = (Z[n]!=Z[n-1]); } }
            }
            while (n<N) { Z[n] = 0; n += R; }
            for (c=0; c<C; c++) { for (n=c*R+1; n<c*R+R; n++) { Z[n] += Z[n-1]; } } //cumsum
            while (t<T && cs<=Lpre)
            {
                sc = (float)L/(cs+1+Lpost);
                for (c=0; c<C; c++) { Y[t+c*T] = sc * Z[c*R+cs+Lpost]; }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
            while (t<T && cs<R-Lpost)
            {
                for (c=0; c<C; c++) { Y[t+c*T] = Z[c*R+cs+Lpost] - Z[c*R+cs-Lpre-1]; }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
            while (t<T)
            {
                sc = (float)L/(R-cs+Lpre);
                for (c=0; c<C; c++) { Y[t+c*T] = sc * (Z[c*R+R-1] - Z[c*R+cs-Lpre-1]); }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
        }
        else
        {
            for (ns=0; ns<C; ns++)
            {
                m = cblas_sdot(R,&X[ns],C,&o,0) / R;
                n = ns;
                if (going>0) { while (n<N) { Z[n] = (X[n]>=m); n += C; } }
                else { while (n<N) { Z[n] = (X[n]<m); n += C; } }
                n -= C;
                if (going) { while (n>ns) { Z[n] *= (Z[n]!=Z[n-C]); n -= C; } }
                else { while (n>ns) { Z[n] = (Z[n]!=Z[n-C]); n -= C; } }
            }
            while (n>=0) { Z[n--] = 0; }
            n += C; while (++n<N) { Z[n] += Z[n-C]; }  //cumsum
            while (t<T && cs<=Lpre)
            {
                sc = (float)L/(cs+1+Lpost);
                for (c=0; c<C; c++) { Y[c+t*C] = sc * Z[c+(cs+Lpost)*C]; }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
            while (t<T && cs<R-Lpost)
            {
                for (c=0; c<C; c++) { Y[c+t*C] = Z[c+(cs+Lpost)*C] - Z[c+(cs-Lpre-1)*C]; }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
            while (t<T)
            {
                sc = (float)L/(R-cs+Lpre);
                for (c=0; c<C; c++) { Y[c+t*C] = sc * (Z[N-C+c] - Z[c+(cs-Lpre-1)*C]); }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
        }
        cblas_sscal(T*C,1.0f/L,&Y[0],1);
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (ns=0; ns<R; ns++)
            {
                m = cblas_sdot(C,&X[ns],R,&o,0) / C;
                n = ns;
                if (going>0) { while (n<N) { Z[n] = (X[n]>=m); n += R; } }
                else { while (n<N) { Z[n] = (X[n]<m); n += R; } }
                n -= R;
                if (going) { while (n>ns) { Z[n] *= (Z[n]!=Z[n-R]); n -= R; } }
                else { while (n>ns) { Z[n] = (Z[n]!=Z[n-R]); n -= R; }  }
            }
            while (n>=0) { Z[n--] = 0; }
            n += R; while (++n<N) { Z[n] += Z[n-R]; }  //cumsum
            while (t<T && cs<=Lpre)
            {
                sc = (float)L/(cs+1+Lpost);
                for (r=0; r<R; r++) { Y[r+t*R] = sc * Z[r+(cs+Lpost)*R]; }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
            while (t<T && cs<C-Lpost)
            {
                for (r=0; r<R; r++) { Y[r+t*R] = Z[r+(cs+Lpost)*R] - Z[r+(cs-Lpre-1)*R]; }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
            while (t<T)
            {
                sc = (float)L/(C-cs+Lpre);
                for (r=0; r<R; r++) { Y[r+t*R] = sc * (Z[N-R+r] - Z[r+(cs-Lpre-1)*R]); }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
        }
        else
        {
            for (ns=0; ns<N; ns+=C)
            {
                m = cblas_sdot(C,&X[ns],1,&o,0) / C;
                n = ns - 1;
                if (going>0) { while (++n<ns+C) { Z[n] = (X[n]>=m); } }
                else { while (++n<ns+C) { Z[n] = (X[n]<m); } }
                if (going) { while (--n>ns) { Z[n] *= (Z[n]!=Z[n-1]); } }
                else { while (--n>ns) { Z[n] = (Z[n]!=Z[n-1]); } }
            }
            while (n<N) { Z[n] = 0; n += C; }
            for (r=0; r<R; r++) { for (n=r*C+1; n<r*C+C; n++) { Z[n] += Z[n-1]; } }  //cumsum
            while (t<T && cs<=Lpre)
            {
                sc = (float)L/(cs+1+Lpost);
                for (r=0; r<R; r++) { Y[t+r*T] = sc * Z[r*C+cs+Lpost]; }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
            while (t<T && cs<C-Lpost)
            {
                for (r=0; r<R; r++) { Y[t+r*T] = Z[r*C+cs+Lpost] - Z[r*C+cs-Lpre-1]; }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
            while (t<T)
            {
                sc = (float)L/(C-cs+Lpre);
                for (r=0; r<R; r++) { Y[t+r*T] = sc * (Z[r*C+C-1] - Z[r*C+cs-Lpre-1]); }
                t++; cs = (int)(roundf(t*stp)) + c0;
            }
        }
        cblas_sscal(R*T,1.0f/L,&Y[0],1);
    }
    else
    {
        fprintf(stderr,"error in mcr_s: dim must be 0 or 1.\n"); return 1;
    }

    free(Z);
    return 0;
}


int mcr_d (double *Y, const double *X, const int iscolmajor, const int R, const int C, const int L, const int dim, const int c0, const double stp, const int going)
{
    const double o = 1.0;
    const int N = R*C;
    const int T = (dim==0) ? 1 + (int)(floor(((int)R-1-c0)/stp)) : 1 + (int)(floor(((int)C-1-c0)/stp));
    const int Lpre = L/2;          //nsamps before center samp
    const int Lpost = L - L/2 - 1; //nsamps after center samp
    int cs = c0, t = 0;
    int c, r, ns, n = -1;
    int *Z;
    double sc, m;

    //Checks
    if (R<1) { fprintf(stderr,"error in mcr_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in mcr_d: C (ncols X) must be positive\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in mcr_d: L (winlength) must be positive\n"); return 1; }

    //Initialize Z
    if (!(Z=(int *)malloc((size_t)(N)*sizeof(int)))) { fprintf(stderr,"error in mcr_d: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (ns=0; ns<N; ns+=R)
            {
                m = cblas_ddot(R,&X[ns],1,&o,0) / R;
                n = ns - 1;
                if (going>0) { while (++n<ns+R) { Z[n] = (X[n]>=m); } }
                else { while (++n<ns+R) { Z[n] = (X[n]<m); } }
                if (going) { while (--n>ns) { Z[n] *= (Z[n]!=Z[n-1]); } }
                else { while (--n>ns) { Z[n] = (Z[n]!=Z[n-1]); } }
            }
            while (n<N) { Z[n] = 0; n += R; }
            for (c=0; c<C; c++) { for (n=c*R+1; n<c*R+R; n++) { Z[n] += Z[n-1]; } } //cumsum
            while (t<T && cs<=Lpre)
            {
                sc = (double)L/(cs+1+Lpost);
                for (c=0; c<C; c++) { Y[t+c*T] = sc * Z[c*R+cs+Lpost]; }
                t++; cs = (int)(round(t*stp)) + c0;
            }
            while (t<T && cs<R-Lpost)
            {
                for (c=0; c<C; c++) { Y[t+c*T] = Z[c*R+cs+Lpost] - Z[c*R+cs-Lpre-1]; }
                t++; cs = (int)(round(t*stp)) + c0;
            }
            while (t<T)
            {
                sc = (double)L/(R-cs+Lpre);
                for (c=0; c<C; c++) { Y[t+c*T] = sc * (Z[c*R+R-1] - Z[c*R+cs-Lpre-1]); }
                t++; cs = (int)(round(t*stp)) + c0;
            }
        }
        else
        {
            for (ns=0; ns<C; ns++)
            {
                m = cblas_ddot(R,&X[ns],C,&o,0) / R;
                n = ns;
                if (going>0) { while (n<N) { Z[n] = (X[n]>=m); n += C; } }
                else { while (n<N) { Z[n] = (X[n]<m); n += C; } }
                n -= C;
                if (going) { while (n>ns) { Z[n] *= (Z[n]!=Z[n-C]); n -= C; } }
                else { while (n>ns) { Z[n] = (Z[n]!=Z[n-C]); n -= C; } }
            }
            while (n>=0) { Z[n--] = 0; }
            n += C; while (++n<N) { Z[n] += Z[n-C]; }  //cumsum
            while (t<T && cs<=Lpre)
            {
                sc = (double)L/(cs+1+Lpost);
                for (c=0; c<C; c++) { Y[c+t*C] = sc * Z[c+(cs+Lpost)*C]; }
                t++; cs = (int)(round(t*stp)) + c0;
            }
            while (t<T && cs<R-Lpost)
            {
                for (c=0; c<C; c++) { Y[c+t*C] = Z[c+(cs+Lpost)*C] - Z[c+(cs-Lpre-1)*C]; }
                t++; cs = (int)(round(t*stp)) + c0;
            }
            while (t<T)
            {
                sc = (double)L/(R-cs+Lpre);
                for (c=0; c<C; c++) { Y[c+t*C] = sc * (Z[N-C+c] - Z[c+(cs-Lpre-1)*C]); }
                t++; cs = (int)(round(t*stp)) + c0;
            }
        }
        cblas_dscal(T*C,1.0/L,&Y[0],1);
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (ns=0; ns<R; ns++)
            {
                m = cblas_ddot(C,&X[ns],R,&o,0) / C;
                n = ns;
                if (going>0) { while (n<N) { Z[n] = (X[n]>=m); n += R; } }
                else { while (n<N) { Z[n] = (X[n]<m); n += R; } }
                n -= R;
                if (going) { while (n>ns) { Z[n] *= (Z[n]!=Z[n-R]); n -= R; } }
                else { while (n>ns) { Z[n] = (Z[n]!=Z[n-R]); n -= R; } }
            }
            while (n>=0) { Z[n--] = 0; }
            n += R; while (++n<N) { Z[n] += Z[n-R]; }  //cumsum
            while (t<T && cs<=Lpre)
            {
                sc = (double)L/(cs+1+Lpost);
                for (r=0; r<R; r++) { Y[r+t*R] = sc * Z[r+(cs+Lpost)*R]; }
                t++; cs = (int)(round(t*stp)) + c0;
            }
            while (t<T && cs<C-Lpost)
            {
                for (r=0; r<R; r++) { Y[r+t*R] = Z[r+(cs+Lpost)*R] - Z[r+(cs-Lpre-1)*R]; }
                t++; cs = (int)(round(t*stp)) + c0;
            }
            while (t<T)
            {
                sc = (double)L/(C-cs+Lpre);
                for (r=0; r<R; r++) { Y[r+t*R] = sc * (Z[N-R+r] - Z[r+(cs-Lpre-1)*R]); }
                t++; cs = (int)(round(t*stp)) + c0;
            }
        }
        else
        {
            for (ns=0; ns<N; ns+=C)
            {
                m = cblas_ddot(C,&X[ns],1,&o,0) / C;
                n = ns - 1;
                if (going>0) { while (++n<ns+C) { Z[n] = (X[n]>=m); } }
                else { while (++n<ns+C) { Z[n] = (X[n]<m); } }
                if (going) { while (--n>ns) { Z[n] *= (Z[n]!=Z[n-1]); } }
                else { while (--n>ns) { Z[n] = (Z[n]!=Z[n-1]); } }
            }
            while (n<N) { Z[n] = 0; n += C; }
            for (r=0; r<R; r++) { for (n=r*C+1; n<r*C+C; n++) { Z[n] += Z[n-1]; } }  //cumsum
            while (t<T && cs<=Lpre)
            {
                sc = (double)L/(cs+1+Lpost);
                for (r=0; r<R; r++) { Y[t+r*T] = sc * Z[r*C+cs+Lpost]; }
                t++; cs = (int)(round(t*stp)) + c0;
            }
            while (t<T && cs<C-Lpost)
            {
                for (r=0; r<R; r++) { Y[t+r*T] = Z[r*C+cs+Lpost] - Z[r*C+cs-Lpre-1]; }
                t++; cs = (int)(round(t*stp)) + c0;
            }
            while (t<T)
            {
                sc = (double)L/(C-cs+Lpre);
                for (r=0; r<R; r++) { Y[t+r*T] = sc * (Z[r*C+C-1] - Z[r*C+cs-Lpre-1]); }
                t++; cs = (int)(round(t*stp)) + c0;
            }
        }
        cblas_dscal(R*T,1.0/L,&Y[0],1);
    }
    else
    {
        fprintf(stderr,"error in mcr_d: dim must be 0 or 1.\n"); return 1;
    }

    free(Z);
    return 0;
}


#ifdef __cplusplus
}
}
#endif

