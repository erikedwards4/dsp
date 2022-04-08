//Roots of polynomial using compan matrix, as in Matlab/Octave.
//Luckily, the compan matrix is upper Hessenberg,
//and only eigenvalues are required, so can easily use LAPACKE directly.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int roots_s (float *roots, const float *AS, const int P);
int roots_d (double *roots, const double *AS, const int P);
int roots_c (float *roots, const float *AS, const int P);
int roots_z (double *roots, const double *AS, const int P);


int roots_s (float *roots, const float *AS, const int P)
{
    const float o = 1.0f;
    const char job = 'E', compz = 'N';  //eigenvalues only
    const lapack_int ldh = P-1, n = P-1, ldz = 1;
    const lapack_int ilo = 1, ihi = n;  //avoids balancing
    lapack_int info;
    float *compan, *wr, *wi, z[1];
    int p;

    //Checks
    if (P<1) { fprintf(stderr,"error in roots_s: P (length of polynomial coeffs including a0=1) must be positive\n"); return 1; }

    //Make compan matrix
    if (!(wr=(float *)malloc((size_t)(n)*sizeof(float)))) { fprintf(stderr,"error in roots_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(wi=(float *)malloc((size_t)(n)*sizeof(float)))) { fprintf(stderr,"error in roots_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(compan=(float *)calloc((size_t)(n*n),sizeof(float)))) { fprintf(stderr,"error in roots_s: problem with calloc. "); perror("calloc"); return 1; }
    cblas_scopy(P-2,&o,0,&compan[1],P);
    for (p=0; p<P-1; p++) { compan[p*n] = -AS[p+1]/AS[0]; }
    //fprintf(stderr,"compan = \n"); for (p=0; p<n; p++) { for (int p2=0; p2<n; p2++) { fprintf(stderr,"%f ",compan[p+p2*n]); } fprintf(stderr,"\n"); }

    //Do eig
    info = LAPACKE_shseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,z,ldz);
    if (info) { fprintf(stderr,"error in roots_s: lapacke decomposition failed\n"); return 1; }

    //Copy to (complex-valued) roots
    cblas_scopy(P-1,wr,1,&roots[0],2);
    cblas_scopy(P-1,wi,1,&roots[1],2);

    //Exit
    free(compan); free(wr); free(wi);
    return 0;
}


int roots_d (double *roots, const double *AS, const int P)
{
    const double o = 1.0;
    const char job = 'E', compz = 'N';  //eigenvalues only
    const lapack_int ldh = P-1, n = P-1, ldz = 1;
    const lapack_int ilo = 1, ihi = n;  //avoids balancing
    lapack_int info;
    double *compan, *wr, *wi, z[1];
    int p;

    //Checks
    if (P<1) { fprintf(stderr,"error in roots_d: P (length of polynomial coeffs including a0=1) must be positive\n"); return 1; }

    //Make compan matrix
    if (!(wr=(double *)malloc((size_t)(n)*sizeof(double)))) { fprintf(stderr,"error in roots_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(wi=(double *)malloc((size_t)(n)*sizeof(double)))) { fprintf(stderr,"error in roots_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(compan=(double *)calloc((size_t)(n*n),sizeof(double)))) { fprintf(stderr,"error in roots_d: problem with calloc. "); perror("calloc"); return 1; }
    cblas_dcopy(P-2,&o,0,&compan[1],P);
    for (p=0; p<P-1; p++) { compan[p*n] = -AS[p+1]/AS[0]; }

    //Do eig
    info = LAPACKE_dhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,z,ldz);
    if (info) { fprintf(stderr,"error in roots_d: eigendecomposition failed\n"); return 1; }

    //Copy to (complex-valued) roots
    cblas_dcopy(P-1,wr,1,&roots[0],2);
    cblas_dcopy(P-1,wi,1,&roots[1],2);

    //Exit
    free(compan); free(wr); free(wi);
    return 0;
}


int roots_c (float *roots, const float *AS, const int P)
{
    const float o[2] =  {1.0f,0.0f};
    char job = 'E', compz = 'N';  //eigenvalues only
    const lapack_int ldh = P-1, n = P-1, ldz = 1;
    const lapack_int ilo = 1, ihi = n;  //avoids balancing
    lapack_int info;
    float *compan, z[2];
    int p;

    //Checks
    if (P<1) { fprintf(stderr,"error in roots_c: P (length of polynomial coeffs including a0=1) must be positive\n"); return 1; }

    //Make compan matrix
    if (!(compan=(float *)calloc((size_t)(4*n*n),sizeof(float)))) { fprintf(stderr,"error in roots_c: problem with calloc. "); perror("calloc"); return 1; }
    cblas_ccopy(P-2,o,0,&compan[2],P);
    for (p=0; p<P-1; p++) { compan[2*p*n] = -AS[2*p+2]/AS[0]; compan[2*p*n+1] = -AS[2*p+3]/AS[0]; }
    //fprintf(stderr,"compan = \n"); for (p=0; p<P-1; p++) { for (int p2=0; p2<P-1; p2++) { fprintf(stderr,"%f ",compan[p+p2*(P-1)]); } fprintf(stderr,"\n"); }

    //Do eig
    info = LAPACKE_chseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_float *)compan,ldh,(lapack_complex_float *)roots,(lapack_complex_float *)z,ldz);
    if (info) { fprintf(stderr,"error in roots_c: lapacke decomposition failed\n"); return 1; }
    
    //Exit
    free(compan);
    return 0;
}


int roots_z (double *roots, const double *AS, const int P)
{
    const double o[2] =  {1.0,0.0};
    char job = 'E', compz = 'N';  //eigenvalues only
    const lapack_int ldh = P-1, n = P-1, ldz = 1;
    const lapack_int ilo = 1, ihi = n;  //avoids balancing
    lapack_int info;
    double *compan, z[2];
    int p;

    //Checks
    if (P<1) { fprintf(stderr,"error in roots_z: P (length of polynomial coeffs including a0=1) must be positive\n"); return 1; }

    //Make compan matrix
    if (!(compan=(double *)calloc((size_t)(4*n*n),sizeof(double)))) { fprintf(stderr,"error in roots_z: problem with calloc. "); perror("calloc"); return 1; }
    cblas_zcopy(P-2,o,0,&compan[2],P);
    for (p=0; p<P-1; p++) { compan[2*p*n] = -AS[2*p+2]/AS[0]; compan[2*p*n+1] = -AS[2*p+3]/AS[0]; }
    //fprintf(stderr,"compan = \n"); for (p=0; p<P-1; p++) { for (int p2=0; p2<P-1; p2++) { fprintf(stderr,"%f ",compan[p+p2*(P-1)]); } fprintf(stderr,"\n"); }

    //Do eig
    info = LAPACKE_zhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_double *)compan,ldh,(lapack_complex_double *)roots,(lapack_complex_double *)z,ldz);
    if (info) { fprintf(stderr,"error in roots_z: lapacke decomposition failed\n"); return 1; }

    //Exit
    free(compan);
    return 0;
}


#ifdef __cplusplus
}
}
#endif

