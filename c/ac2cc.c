//Gets cepstral coefficients (CCs) from the autocorrelation (AC) function for each row or col of AC.
//Starts with a Levinson-Durbin recursion from the AC values.

//To test compile:
//gcc -c ac2ccs.c -O2 -std=c99 -Wall -Wextra
//clang -c ac2ccs.c -O2 -std=c99 -Weverything
//g++ -c ac2ccs.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c ac2ccs.c -O2 -std=c++11 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int ac2ccs_s (float *CC, const float *AC, const char iscolmajor, const int R, const int C, const int dim, const int K, const float preg)
{
    const float o = 1.0f;
    const int P = K - 1;
    int r, c, p, q, j, k;
	float g, v, sm;
    float *AS, *AStmp;

    //Checks
    if (dim==0 && P>R) { fprintf(stderr,"error in ac2ccs_s: P must be < nrows AC for dim==0\n"); return 1; }
    if (dim==1 && P>C) { fprintf(stderr,"error in ac2ccs_s: P must be < ncols AC for dim==1\n"); return 1; }
	
    //Initialize
    if (!(AS=(float *)malloc((size_t)P*sizeof(float)))) { fprintf(stderr,"error in ac2ccs_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(AStmp=(float *)malloc((size_t)(P-1)*sizeof(float)))) { fprintf(stderr,"error in ac2ccs_s: problem with malloc. "); perror("malloc"); return 1; }
	
    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_scopy(C,&o,0,&AS[0],P);
            for (c=0; c<C; c++)
            {
                AS[1] = AStmp[0] = g = -AC[1+c*R]/AC[c*R];
                v = fmaf(AC[1+c*R],g,AC[c*R]);
                for (p=2; p<P; p++)
                {
                    g = AC[p+c*R];
                    for (q=1; q<p; q++) { g = fmaf(AC[q+c*R],AS[p-q],g); }
                    AS[p] = g = -g/v;
                    for (q=1; q<p; q++) { AS[q] = fmaf(g,AStmp[p-q-1],AS[q]); }
                    cblas_scopy(p,&AS[1],1,&AStmp[0],1);
                    v *= fmaf(g,-g,1.0f);
                }
	            CC[c*K] = logf(v+preg);
                for (k=1; k<K; k++)
                {
                    sm = 0.0f;
                    for (j=1; j<=k; j++) { sm = fmaf(AS[j]*(j-k),CC[k-j+c*K],sm); }
                    CC[k+c*K] = sm/k - AS[k];
                }
            }
        }
        else
        {
            cblas_scopy(C,&o,0,&AS[0],1);
            for (c=0; c<C; c++)
            {
                AS[1] = AStmp[0] = g = -AC[c+C]/AC[c];
                v = fmaf(AC[c+C],g,AC[c]);
                for (p=2; p<P; p++)
                {
                    g = AC[c+p*C];
                    for (q=1; q<p; q++) { g = fmaf(AC[c+q*C],AS[p-q],g); }
                    AS[p] = g = -g/v;
                    for (q=1; q<p; q++) { AS[q] = fmaf(g,AStmp[p-q-1],AS[q]); }
                    cblas_scopy(p,&AS[1],C,&AStmp[0],1);
                    v *= fmaf(g,-g,1.0f);
                }
                CC[c] = logf(v+preg);
                for (k=1; k<K; k++)
                {
                    sm = 0.0f;
                    for (j=1; j<=k; j++) { sm = fmaf(AS[j]*(j-k),CC[c+(k-j)*C],sm); }
                    CC[c+k*C] = sm/k - AS[k];
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_scopy(R,&o,0,&AS[0],1);
            for (r=0; r<R; r++)
            {
                AS[1] = AStmp[0] = g = -AC[r+R]/AC[r];
                v = fmaf(AC[r+R],g,AC[r]);
                for (p=2; p<P; p++)
                {
                    g = AC[r+p*R];
                    for (q=1; q<p; q++) { g = fmaf(AC[r+q*R],AS[p-q],g); }
                    AS[p] = g = -g/v;
                    for (q=1; q<p; q++) { AS[q] = fmaf(g,AStmp[p-q-1],AS[q]); }
                    cblas_scopy(p,&AS[1],R,&AStmp[0],1);
                    v *= fmaf(g,-g,1.0f);
                }
                CC[r] = logf(v+preg);
                for (k=1; k<K; k++)
                {
                    sm = 0.0f;
                    for (j=1; j<=k; j++) { sm = fmaf(AS[j]*(j-k),CC[r+(k-j)*R],sm); }
                    CC[r+k*R] = sm/k - AS[k];
                }
            }
        }
        else
        {
            cblas_scopy(R,&o,0,&AS[0],P);
            for (r=0; r<R; r++)
            {
                AS[1] = AStmp[0] = g = -AC[1+r*C]/AC[r*C];
                v = fmaf(AC[1+r*C],g,AC[r*C]);
                for (p=2; p<P; p++)
                {
                    g = AC[p+r*C];
                    for (q=1; q<p; q++) { g = fmaf(AC[q+r*C],AS[p-q],g); }
                    AS[p] = g = -g/v;
                    for (q=1; q<p; q++) { AS[q] = fmaf(g,AStmp[p-q-1],AS[q]); }
                    cblas_scopy(p,&AS[1],1,&AStmp[0],1);
                    v *= fmaf(g,-g,1.0f);
                }
                CC[r*K] = logf(v+preg);
                for (k=1; k<K; k++)
                {
                    sm = 0.0f;
                    for (j=1; j<=k; j++) { sm = fmaf(AS[j]*(j-k),CC[k-j+r*K],sm); }
                    CC[k+r*K] = sm/k - AS[k];
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in ac2ccs_s: dim must be 0 or 1.\n"); return 1;
    }
	
	return 0;
}


int ac2ccs_d (double *CC, const double *AC, const char iscolmajor, const int R, const int C, const int dim, const int K, const double preg)
{
	const double o = 1.0;
    const int P = K - 1;
    int r, c, p, q, j, k;
	double g, v, sm;
    double *AS, *AStmp;

    //Checks
    if (dim==0 && P>R) { fprintf(stderr,"error in ac2ccs_d: P must be < nrows AC for dim==0\n"); return 1; }
    if (dim==1 && P>C) { fprintf(stderr,"error in ac2ccs_d: P must be < ncols AC for dim==1\n"); return 1; }
	
    //Initialize
    if (!(AS=(double *)malloc((size_t)P*sizeof(double)))) { fprintf(stderr,"error in ac2ccs_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(AStmp=(double *)malloc((size_t)(P-1)*sizeof(double)))) { fprintf(stderr,"error in ac2ccs_d: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_dcopy(C,&o,0,&AS[0],P);
            for (c=0; c<C; c++)
            {
                AS[1] = AStmp[0] = g = -AC[1+c*R]/AC[c*R];
                v = fma(AC[1+c*R],g,AC[c*R]);
                for (p=2; p<P; p++)
                {
                    g = AC[p+c*R];
                    for (q=1; q<p; q++) { g = fma(AC[q+c*R],AS[p-q],g); }
                    AS[p] = g = -g/v;
                    for (q=1; q<p; q++) { AS[q] = fma(g,AStmp[p-q-1],AS[q]); }
                    cblas_dcopy(p,&AS[1],1,&AStmp[0],1);
                    v *= fma(g,-g,1.0);
                }
                CC[c*K] = log(v+preg);
                for (k=1; k<K; k++)
                {
                    sm = 0.0;
                    for (j=1; j<=k; j++) { sm = fma(AS[j]*(j-k),CC[k-j+c*K],sm); }
                    CC[k+c*K] = sm/k - AS[k];
                }
            }
        }
        else
        {
            cblas_dcopy(C,&o,0,&AS[0],1);
            for (c=0; c<C; c++)
            {
                AS[1] = AStmp[0] = g = -AC[c+C]/AC[c];
                v = fma(AC[c+C],g,AC[c]);
                for (p=2; p<P; p++)
                {
                    g = AC[c+p*C];
                    for (q=1; q<p; q++) { g = fma(AC[c+q*C],AS[p-q],g); }
                    AS[p] = g = -g/v;
                    for (q=1; q<p; q++) { AS[q] = fma(g,AStmp[p-q-1],AS[q]); }
                    cblas_dcopy(p,&AS[1],C,&AStmp[0],1);
                    v *= fma(g,-g,1.0);
                }
                CC[c] = log(v+preg);
                for (k=1; k<K; k++)
                {
                    sm = 0.0;
                    for (j=1; j<=k; j++) { sm = fma(AS[j]*(j-k),CC[c+(k-j)*C],sm); }
                    CC[c+k*C] = sm/k - AS[k];
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_dcopy(R,&o,0,&AS[0],1);
            for (r=0; r<R; r++)
            {
                AS[1] = AStmp[0] = g = -AC[r+R]/AC[r];
                v = fma(AC[r+R],g,AC[r]);
                for (p=2; p<P; p++)
                {
                    g = AC[r+p*R];
                    for (q=1; q<p; q++) { g = fma(AC[r+q*R],AS[p-q],g); }
                    AS[p] = g = -g/v;
                    for (q=1; q<p; q++) { AS[q] = fma(g,AStmp[p-q-1],AS[q]); }
                    cblas_dcopy(p,&AS[1],R,&AStmp[0],1);
                    v *= fma(g,-g,1.0);
                }
                CC[r] = log(v+preg);
                for (k=1; k<K; k++)
                {
                    sm = 0.0;
                    for (j=1; j<=k; j++) { sm = fma(AS[j]*(j-k),CC[r+(k-j)*R],sm); }
                    CC[r+k*R] = sm/k - AS[k];
                }
            }
        }
        else
        {
            cblas_dcopy(R,&o,0,&AS[0],P);
            for (r=0; r<R; r++)
            {
                AS[1] = AStmp[0] = g = -AC[1+r*C]/AC[r*C];
                v = fma(AC[1+r*C],g,AC[r*C]);
                for (p=2; p<P; p++)
                {
                    g = AC[p+r*C];
                    for (q=1; q<p; q++) { g = fma(AC[q+r*C],AS[p-q],g); }
                    AS[p] = g = -g/v;
                    for (q=1; q<p; q++) { AS[q] = fma(g,AStmp[p-q-1],AS[q]); }
                    cblas_dcopy(p,&AS[1],1,&AStmp[0],1);
                    v *= fma(g,-g,1.0);
                }
                CC[r*K] = log(v+preg);
                for (k=1; k<K; k++)
                {
                    sm = 0.0;
                    for (j=1; j<=k; j++) { sm = fma(AS[j]*(j-k),CC[k-j+r*K],sm); }
                    CC[k+r*K] = sm/k - AS[k];
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in ac2ccs_d: dim must be 0 or 1.\n"); return 1;
    }
	
	return 0;
}


#ifdef __cplusplus
}
#endif

