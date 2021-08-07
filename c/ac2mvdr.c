//Gets multivariate distortionless response (MVDR) starting from the autocorrelation (AC) function.
//Starts with a Levinson-Durbin recursion from the AC values.
//Then does FFT according to Marple (use only real part of output, so I use FFTW_R2HC).
//Also note that Marple used the Burg method (here I use Levinson-Durbin recursion).

//Input F is the usual nfft/2 + 1, and is the size of output MVDR along dimension dim.
//Marple [p. 359] uses a large nfft of fixed size 4096, and requires nfft is a power of 2.
//I believe that nfft must be > 2*P, where P = nlags -1 is the AR model order.
//Thus, set nfft to nextpow2(2*R) for dim==0 or nextpow2(2*C) for dim==1.
//Then take usual F = nfft/2 + 1, and set MVDR size to be FxC or RxF.
//The output is then on a linear frequency scale starting at 0 Hz.

//For the PMVDR, I had used
//nfft = 4*B; //this ensures that nfft>2*M for all practical choices of B, but keeps correct freq scale at output

//To test compile:
//gcc -c ac2mvdr.c -O2 -std=c99 -Wall -Wextra
//clang -c ac2mvdr.c -O2 -std=c99 -Weverything
//g++ -c ac2mvdr.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c ac2mvdr.c -O2 -std=c++11 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int ac2mvdr_s (float *MVDR, const float *AC, const char iscolmajor, const int R, const int C, const int dim, const int F, const float preg)
{
    const float z = 0.0f, o = 1.0f;
    const int P = (dim==0) ? R-1 : C-1;
    int r, c, p, q, f, nfft;
	float g, v;
    float *AS, *AStmp, *X1, *Y1;
    fftwf_plan fplan;

    //Checks
    if (dim==0 && F<R) { fprintf(stderr,"error in ac2mvdr_s: F must be >= nlags AC\n"); return 1; }
    if (dim==1 && F<C) { fprintf(stderr,"error in ac2mvdr_s: F must be >= nlags AC\n"); return 1; }
	
    //Initialize LP
    if (!(AS=(float *)malloc((size_t)P*sizeof(float)))) { fprintf(stderr,"error in ac2mvdr_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(AStmp=(float *)malloc((size_t)(P-1)*sizeof(float)))) { fprintf(stderr,"error in ac2mvdr_s: problem with malloc. "); perror("malloc"); return 1; }

    //Initialize FFT
    nfft = 1; while (nfft<F) { nfft *= 2; }  //assumes F=nfft/2+1, where nfft is a power of 2
    X1 = fftwf_alloc_real((size_t)nfft);
    Y1 = fftwf_alloc_real((size_t)nfft);
    fplan = fftwf_plan_r2r_1d(nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
    if (!fplan) { fprintf(stderr,"error in ac2mvdr_s: problem creating fftw plan"); return 1; }

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
                cblas_scopy(nfft-P,&z,0,&X1[0],1);
                for (p=0; p<=P; p++)   //enter mus directly into X1
                {
                    for (q=0; q<=P-p; q++) { X1[p] += (P+1-p-2*q)*AS[q]*AS[q+p]; }
                    if (p>0) { X1[nfft-p] = X1[p]; }
                }
                fftwf_execute(fplan);  //Marple FFT part (but may need to scale by fs and/or nfft)
                for (f=0; f<F; f++) { MVDR[f+c*F] = v/(Y1[f]+preg); }
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
                cblas_scopy(nfft-P,&z,0,&X1[0],1);
                for (p=0; p<=P; p++)   //enter mus directly into X1
                {
                    for (q=0; q<=P-p; q++) { X1[p] += (P+1-p-2*q)*AS[q]*AS[q+p]; }
                    if (p>0) { X1[nfft-p] = X1[p]; }
                }
                fftwf_execute(fplan);  //Marple FFT part (but may need to scale by fs and/or nfft)
                for (f=0; f<F; f++) { MVDR[c+f*C] = v/(Y1[f]+preg); }
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
                cblas_scopy(nfft-P,&z,0,&X1[0],1);
                for (p=0; p<=P; p++)   //enter mus directly into X1
                {
                    for (q=0; q<=P-p; q++) { X1[p] += (P+1-p-2*q)*AS[q]*AS[q+p]; }
                    if (p>0) { X1[nfft-p] = X1[p]; }
                }
                fftwf_execute(fplan);  //Marple FFT part (but may need to scale by fs and/or nfft)
                for (f=0; f<F; f++) { MVDR[r+f*R] = v/(Y1[f]+preg); }
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
                cblas_scopy(nfft-P,&z,0,&X1[0],1);
                for (p=0; p<=P; p++)   //enter mus directly into X1
                {
                    for (q=0; q<=P-p; q++) { X1[p] += (P+1-p-2*q)*AS[q]*AS[q+p]; }
                    if (p>0) { X1[nfft-p] = X1[p]; }
                }
                fftwf_execute(fplan);  //Marple FFT part (but may need to scale by fs and/or nfft)
                for (f=0; f<F; f++) { MVDR[f+r*F] = v/(Y1[f]+preg); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in ac2mvdr_s: dim must be 0 or 1.\n"); return 1;
    }
	
	return 0;
}


int ac2mvdr_d (double *MVDR, const double *AC, const char iscolmajor, const int R, const int C, const int dim, const int F, const double preg)
{
	const double z = 0.0, o = 1.0;
    const int P = (dim==0) ? R-1 : C-1;
    int r, c, p, q, f, nfft;
	double g, v;
    double *AS, *AStmp, *X1, *Y1;
    fftw_plan fplan;

    //Checks
    if (dim==0 && F<R) { fprintf(stderr,"error in ac2mvdr_d: F must be >= nlags AC\n"); return 1; }
    if (dim==1 && F<C) { fprintf(stderr,"error in ac2mvdr_d: F must be >= nlags AC\n"); return 1; }
	
    //Initialize LP
    if (!(AS=(double *)malloc((size_t)P*sizeof(double)))) { fprintf(stderr,"error in ac2mvdr_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(AStmp=(double *)malloc((size_t)(P-1)*sizeof(double)))) { fprintf(stderr,"error in ac2mvdr_d: problem with malloc. "); perror("malloc"); return 1; }
	
    //Initialize FFT
    nfft = 1; while (nfft<F) { nfft *= 2; }  //assumes F=nfft/2+1, where nfft is a power of 2
    X1 = fftw_alloc_real((size_t)nfft);
    Y1 = fftw_alloc_real((size_t)nfft);
    fplan = fftw_plan_r2r_1d(nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
    if (!fplan) { fprintf(stderr,"error in ac2mvdr_d: problem creating fftw plan"); return 1; }

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
                cblas_dcopy(nfft-P,&z,0,&X1[0],1);
                for (p=0; p<=P; p++)   //enter mus directly into X1
                {
                    for (q=0; q<=P-p; q++) { X1[p] += (P+1-p-2*q)*AS[q]*AS[q+p]; }
                    if (p>0) { X1[nfft-p] = X1[p]; }
                }
                fftw_execute(fplan);  //Marple FFT part (but may need to scale by fs and/or nfft)
                for (f=0; f<F; f++) { MVDR[f+c*F] = v/(Y1[f]+preg); }
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
                cblas_dcopy(nfft-P,&z,0,&X1[0],1);
                for (p=0; p<=P; p++)   //enter mus directly into X1
                {
                    for (q=0; q<=P-p; q++) { X1[p] += (P+1-p-2*q)*AS[q]*AS[q+p]; }
                    if (p>0) { X1[nfft-p] = X1[p]; }
                }
                fftw_execute(fplan);  //Marple FFT part (but may need to scale by fs and/or nfft)
                for (f=0; f<F; f++) { MVDR[c+f*C] = v/(Y1[f]+preg); }
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
                cblas_dcopy(nfft-P,&z,0,&X1[0],1);
                for (p=0; p<=P; p++)   //enter mus directly into X1
                {
                    for (q=0; q<=P-p; q++) { X1[p] += (P+1-p-2*q)*AS[q]*AS[q+p]; }
                    if (p>0) { X1[nfft-p] = X1[p]; }
                }
                fftw_execute(fplan);  //Marple FFT part (but may need to scale by fs and/or nfft)
                for (f=0; f<F; f++) { MVDR[r+f*R] = v/(Y1[f]+preg); }
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
                cblas_dcopy(nfft-P,&z,0,&X1[0],1);
                for (p=0; p<=P; p++)   //enter mus directly into X1
                {
                    for (q=0; q<=P-p; q++) { X1[p] += (P+1-p-2*q)*AS[q]*AS[q+p]; }
                    if (p>0) { X1[nfft-p] = X1[p]; }
                }
                fftw_execute(fplan);  //Marple FFT part (but may need to scale by fs and/or nfft)
                for (f=0; f<F; f++) { MVDR[f+r*F] = v/(Y1[f]+preg); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in ac2mvdr_d: dim must be 0 or 1.\n"); return 1;
    }
	
	return 0;
}


#ifdef __cplusplus
}
#endif

