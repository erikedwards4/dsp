//Roots of polynomial using compan matrix, as in Matlab/Octave.
//This does roots for each row or col of X

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int poly2roots_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int poly2roots_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int poly2roots_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int poly2roots_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int poly2roots_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in poly2roots_s: dim must be in [0 3]\n"); return 1; }
    if (P<1) { fprintf(stderr,"error in poly2roots_s: P (length of polynomial coeffs including a0=1) must be positive\n"); return 1; }

    const float z = 0.0f, o = 1.0f;
    const size_t P = (dim==0) ? R : C;
    const char job = 'E', compz = 'N';  //eigenvalues only
    const lapack_int ldh = P-1, n = P-1, ldz = 1;
    const lapack_int ilo = 1, ihi = n;  //avoids balancing
    lapack_int info;
    float *compan, *wr, *wi, zz[1];
    if (!(wr=(float *)malloc((size_t)(n)*sizeof(float)))) { fprintf(stderr,"error in poly2roots_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(wi=(float *)malloc((size_t)(n)*sizeof(float)))) { fprintf(stderr,"error in poly2roots_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(compan=(float *)malloc((size_t)(n*n)*sizeof(float)))) { fprintf(stderr,"error in poly2roots_s: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_scopy(n*n,&z,0,compan,1); cblas_scopy(n-1,&o,0,&compan[1],P);
                cblas_scopy(n,&X[c*R+1],1,compan,n); cblas_sscal(n,-1.0f/X[c*R],compan,n);
                //fprintf(stderr,"compan = \n"); for (size_t r=0; r<n; r++) { for (c2=0; c2<n; c2++) { fprintf(stderr,"%f ",compan[r+c2*n]); } fprintf(stderr,"\n"); }
                info = LAPACKE_shseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
                if (info) { fprintf(stderr,"error in poly2roots_s: lapacke decomposition failed\n"); return 1; }
                cblas_scopy(n,wr,1,&Y[2*c*n],2); cblas_scopy(n,wi,1,&Y[2*c*n+1],2);  //copy to complex-valued roots
            }
        }
        else
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_scopy(n*n,&z,0,compan,1); cblas_scopy(n-1,&o,0,&compan[1],P);
                cblas_scopy(n,&X[c+C],(int)C,compan,n); cblas_sscal(n,-1.0f/X[c],compan,n);
                //cblas_scopy(n*n,&z,0,compan,1); cblas_scopy(n-1,&o,0,&compan[n],P);
                //cblas_scopy(n,&X[c+C],(int)C,compan,1); cblas_sscal(n,-1.0f/X[c],compan,1);
                //info = LAPACKE_shseqr(LAPACK_ROW_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
                info = LAPACKE_shseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
                if (info) { fprintf(stderr,"error in poly2roots_s: lapacke decomposition failed\n"); return 1; }
                cblas_scopy(n,wr,1,&Y[2*c],2*(int)C); cblas_scopy(n,wi,1,&Y[2*c+1],2*(int)C);  //copy to complex-valued roots
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_scopy(n*n,&z,0,compan,1); cblas_scopy(n-1,&o,0,&compan[1],P);
                cblas_scopy(n,&X[r+R],(int)R,compan,n); cblas_sscal(n,-1.0f/X[r],compan,n);
                info = LAPACKE_shseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
                if (info) { fprintf(stderr,"error in poly2roots_s: lapacke decomposition failed\n"); return 1; }
                cblas_scopy(n,wr,1,&Y[2*r],2*(int)R); cblas_scopy(n,wi,1,&Y[2*r+1],2*(int)R);  //copy to complex-valued roots
            }
        }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_scopy(n*n,&z,0,compan,1); cblas_scopy(n-1,&o,0,&compan[1],P);
                cblas_scopy(n,&X[r*C+1],1,compan,n); cblas_sscal(n,-1.0f/X[r*C],compan,n);
                info = LAPACKE_shseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
                if (info) { fprintf(stderr,"error in poly2roots_s: lapacke decomposition failed\n"); return 1; }
                cblas_scopy(n,wr,1,&Y[2*r*n],2); cblas_scopy(n,wi,1,&Y[2*r*n+1],2);  //copy to complex-valued roots
            }
        }
    }
    else
    {
        fprintf(stderr,"error in poly2roots_s: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(compan); free(wr); free(wi);
    return 0;
}


int poly2roots_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in poly2roots_d: dim must be in [0 3]\n"); return 1; }
    if (P<1) { fprintf(stderr,"error in poly2roots_d: P (length of polynomial coeffs including a0=1) must be positive\n"); return 1; }

    const double z = 0.0, o = 1.0;
    const size_t P = (dim==0) ? R : C;
    const char job = 'E', compz = 'N';  //eigenvalues only
    const lapack_int ldh = P-1, n = P-1, ldz = 1;
    const lapack_int ilo = 1, ihi = n;  //avoids balancing
    lapack_int info;
    double *compan, *wr, *wi, zz[1];
    if (!(wr=(double *)malloc((size_t)(n)*sizeof(double)))) { fprintf(stderr,"error in poly2roots_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(wi=(double *)malloc((size_t)(n)*sizeof(double)))) { fprintf(stderr,"error in poly2roots_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(compan=(double *)malloc((size_t)(n*n)*sizeof(double)))) { fprintf(stderr,"error in poly2roots_d: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_dcopy(n*n,&z,0,compan,1); cblas_dcopy(n-1,&o,0,&compan[1],P);
                cblas_dcopy(n,&X[c*R+1],1,compan,n); cblas_dscal(n,-1.0/X[c*R],compan,n);
                info = LAPACKE_dhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
                if (info) { fprintf(stderr,"error in poly2roots_d: lapacke decomposition failed\n"); return 1; }
                cblas_dcopy(n,wr,1,&Y[2*c*n],2); cblas_dcopy(n,wi,1,&Y[2*c*n+1],2);  //copy to complex-valued roots
            }
        }
        else
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_dcopy(n*n,&z,0,compan,1); cblas_dcopy(n-1,&o,0,&compan[1],P);
                cblas_dcopy(n,&X[c+C],(int)C,compan,n); cblas_dscal(n,-1.0/X[c],compan,n);
                info = LAPACKE_dhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
                if (info) { fprintf(stderr,"error in poly2roots_d: lapacke decomposition failed\n"); return 1; }
                cblas_dcopy(n,wr,1,&Y[2*c],2*(int)C); cblas_dcopy(n,wi,1,&Y[2*c+1],2*(int)C);  //copy to complex-valued roots
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_dcopy(n*n,&z,0,compan,1); cblas_dcopy(n-1,&o,0,&compan[1],P);
                cblas_dcopy(n,&X[r+R],(int)R,compan,n); cblas_dscal(n,-1.0/X[r],compan,n);
                info = LAPACKE_dhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
                if (info) { fprintf(stderr,"error in poly2roots_d: lapacke decomposition failed\n"); return 1; }
                cblas_dcopy(n,wr,1,&Y[2*r],2*(int)R); cblas_dcopy(n,wi,1,&Y[2*r+1],2*(int)R);  //copy to complex-valued roots
            }
        }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_dcopy(n*n,&z,0,compan,1); cblas_dcopy(n-1,&o,0,&compan[1],P);
                cblas_dcopy(n,&X[r*C+1],1,compan,n); cblas_dscal(n,-1.0/X[r*C],compan,n);
                info = LAPACKE_dhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
                if (info) { fprintf(stderr,"error in poly2roots_d: lapacke decomposition failed\n"); return 1; }
                cblas_dcopy(n,wr,1,&Y[2*r*n],2); cblas_dcopy(n,wi,1,&Y[2*r*n+1],2);  //copy to complex-valued roots
            }
        }
    }
    else
    {
        fprintf(stderr,"error in poly2roots_d: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(compan); free(wr); free(wi);
    return 0;
}


int poly2roots_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in poly2roots_c: dim must be in [0 3]\n"); return 1; }
    if (P<1) { fprintf(stderr,"error in poly2roots_c: P (length of polynomial coeffs including a0=1) must be positive\n"); return 1; }

    const float z[2] =  {0.0f,0.0f}, o[2] =  {1.0f,0.0f};
    const size_t P = (dim==0) ? R : C;
    const char job = 'E', compz = 'N';  //eigenvalues only
    const lapack_int ldh = P-1, n = P-1, ldz = 1;
    const lapack_int ilo = 1, ihi = n;  //avoids balancing
    lapack_int info;
    float *compan, *roots, zz[2], sc[2];
    if (!(roots=(float *)malloc((size_t)(2*n)*sizeof(float)))) { fprintf(stderr,"error in poly2roots_c: problem with malloc. "); perror("malloc"); return 1; }
    if (!(compan=(float *)malloc((size_t)(4*n*n)*sizeof(float)))) { fprintf(stderr,"error in poly2roots_c: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_ccopy(n*n,z,0,compan,1); cblas_ccopy(n-1,o,0,&compan[2],P);
                cblas_ccopy(n,&X[2*c*R+2],1,compan,n);
                sc[0] = 1.0f/(X[2*c*R]*X[2*c*R]+X[2*c*R+1]*X[2*c*R+1]); sc[1] = X[2*c*R+1]*sc[0]; sc[0] *= -X[2*c*R];
                cblas_cscal(n,sc,compan,n);
                info = LAPACKE_chseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_float *)compan,ldh,(lapack_complex_float *)&Y[2*c*n],(lapack_complex_float *)zz,ldz);
                if (info) { fprintf(stderr,"error in poly2roots_c: lapacke decomposition failed\n"); return 1; }
                //cblas_ccopy(n,roots,1,&Y[2*c*n],1);
            }
        }
        else
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_ccopy(n*n,z,0,compan,1); cblas_ccopy(n-1,o,0,&compan[2],P);
                cblas_ccopy(n,&X[2*(c+C)],(int)C,compan,n);
                sc[0] = 1.0f/(X[2*c]*X[2*c]+X[2*c+1]*X[2*c+1]); sc[1] = X[2*c+1]*sc[0]; sc[0] *= -X[2*c];
                cblas_cscal(n,sc,compan,n);
                info = LAPACKE_chseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_float *)compan,ldh,(lapack_complex_float *)roots,(lapack_complex_float *)zz,ldz);
                if (info) { fprintf(stderr,"error in poly2roots_c: lapacke decomposition failed\n"); return 1; }
                cblas_ccopy(n,roots,1,&Y[2*c],(int)C);
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_ccopy(n*n,z,0,compan,1); cblas_ccopy(n-1,o,0,&compan[2],P);
                cblas_ccopy(n,&X[2*(r+R)],(int)R,compan,n);
                sc[0] = 1.0f/(X[2*r]*X[2*r]+X[2*r+1]*X[2*r+1]); sc[1] = X[2*r+1]*sc[0]; sc[0] *= -X[2*r];
                cblas_cscal(n,sc,compan,n);
                info = LAPACKE_chseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_float *)compan,ldh,(lapack_complex_float *)roots,(lapack_complex_float *)zz,ldz);
                if (info) { fprintf(stderr,"error in poly2roots_c: lapacke decomposition failed\n"); return 1; }
                cblas_ccopy(n,roots,1,&Y[2*r],(int)R);
            }
        }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_ccopy(n*n,z,0,compan,1); cblas_ccopy(n-1,o,0,&compan[2],P);
                cblas_ccopy(n,&X[2*r*C+2],1,compan,n);
                sc[0] = 1.0f/(X[2*r*C]*X[2*r*C]+X[2*r*C+1]*X[2*r*C+1]); sc[1] = X[2*r*C+1]*sc[0]; sc[0] *= -X[2*r*C];
                cblas_cscal(n,sc,compan,n);
                info = LAPACKE_chseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_float *)compan,ldh,(lapack_complex_float *)&Y[2*r*n],(lapack_complex_float *)zz,ldz);
                if (info) { fprintf(stderr,"error in poly2roots_c: lapacke decomposition failed\n"); return 1; }
                //cblas_ccopy(n,roots,1,&Y[2*r*n],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in poly2roots_c: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(compan); free(roots);
    return 0;
}


int poly2roots_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in poly2roots_z: dim must be in [0 3]\n"); return 1; }
    if (P<1) { fprintf(stderr,"error in poly2roots_z: P (length of polynomial coeffs including a0=1) must be positive\n"); return 1; }
    
    const double z[2] =  {0.0,0.0}, o[2] =  {1.0,0.0};
    const size_t P = (dim==0) ? R : C;
    const char job = 'E', compz = 'N';  //eigenvalues only
    const lapack_int ldh = P-1, n = P-1, ldz = 1;
    const lapack_int ilo = 1, ihi = n;  //avoids balancing
    lapack_int info;
    double *compan, *roots, zz[2], sc[2];
    if (!(roots=(double *)malloc((size_t)(2*n)*sizeof(double)))) { fprintf(stderr,"error in poly2roots_z: problem with malloc. "); perror("malloc"); return 1; }
    if (!(compan=(double *)malloc((size_t)(4*n*n)*sizeof(double)))) { fprintf(stderr,"error in poly2roots_z: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_zcopy(n*n,z,0,compan,1); cblas_zcopy(n-1,o,0,&compan[2],P);
                cblas_zcopy(n,&X[2*c*R+2],1,compan,n);
                sc[0] = 1.0/(X[2*c*R]*X[2*c*R]+X[2*c*R+1]*X[2*c*R+1]); sc[1] = X[2*c*R+1]*sc[0]; sc[0] *= -X[2*c*R];
                cblas_zscal(n,sc,compan,n);
                info = LAPACKE_zhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_double *)compan,ldh,(lapack_complex_double *)&Y[2*c*n],(lapack_complex_double *)zz,ldz);
                if (info) { fprintf(stderr,"error in poly2roots_z: lapacke decomposition failed\n"); return 1; }
                //cblas_zcopy(n,roots,1,&Y[2*c*n],1);
            }
        }
        else
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_zcopy(n*n,z,0,compan,1); cblas_zcopy(n-1,o,0,&compan[2],P);
                cblas_zcopy(n,&X[2*(c+C)],(int)C,compan,n);
                sc[0] = 1.0/(X[2*c]*X[2*c]+X[2*c+1]*X[2*c+1]); sc[1] = X[2*c+1]*sc[0]; sc[0] *= -X[2*c];
                cblas_zscal(n,sc,compan,n);
                info = LAPACKE_zhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_double *)compan,ldh,(lapack_complex_double *)roots,(lapack_complex_double *)zz,ldz);
                if (info) { fprintf(stderr,"error in poly2roots_z: lapacke decomposition failed\n"); return 1; }
                cblas_zcopy(n,roots,1,&Y[2*c],(int)C);
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_zcopy(n*n,z,0,compan,1); cblas_zcopy(n-1,o,0,&compan[2],P);
                cblas_zcopy(n,&X[2*(r+R)],(int)R,compan,n);
                sc[0] = 1.0/(X[2*r]*X[2*r]+X[2*r+1]*X[2*r+1]); sc[1] = X[2*r+1]*sc[0]; sc[0] *= -X[2*r];
                cblas_zscal(n,sc,compan,n);
                info = LAPACKE_zhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_double *)compan,ldh,(lapack_complex_double *)roots,(lapack_complex_double *)zz,ldz);
                if (info) { fprintf(stderr,"error in poly2roots_z: lapacke decomposition failed\n"); return 1; }
                cblas_zcopy(n,roots,1,&Y[2*r],(int)R);
            }
        }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_zcopy(n*n,z,0,compan,1); cblas_zcopy(n-1,o,0,&compan[2],P);
                cblas_zcopy(n,&X[2*r*C+2],1,compan,n);
                sc[0] = 1.0/(X[2*r*C]*X[2*r*C]+X[2*r*C+1]*X[2*r*C+1]); sc[1] = X[2*r*C+1]*sc[0]; sc[0] *= -X[2*r*C];
                cblas_zscal(n,sc,compan,n);
                info = LAPACKE_zhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_double *)compan,ldh,(lapack_complex_double *)&Y[2*r*n],(lapack_complex_double *)zz,ldz);
                if (info) { fprintf(stderr,"error in poly2roots_z: lapacke decomposition failed\n"); return 1; }
                //cblas_zcopy(n,roots,1,&Y[2*r*n],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in poly2roots_z: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(compan); free(roots);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
