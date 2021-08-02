//Roots of polynomial using compan matrix, as in Matlab/Octave.
//This does roots for each vector X.
//The roots are from an eig using LAPACKE,
//and are tested to match Octave, except sometimes for sorting order.

//For the complex case with non-contiguous vecs in X and Y,
//the answers are correct except for a mysterious sign change.

#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int poly2roots_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int poly2roots_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int poly2roots_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int poly2roots_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);


int poly2roots_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in poly2roots_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u || Lx<2u) {}
    else if (Lx==2u)
    {
        if (Lx==N) { *Y = -*(X+1); *++Y = 0.0f; }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, ++X, ++Y) { ++X; *Y = -*X; *++Y = 0.0f; }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=B)
                {
                    for (size_t b=0u; b<B; ++b, ++X, ++Y) { *Y = -*(X+K); *++Y = 0.0f; }
                }
            }
        }
    }
    else
    {
        const size_t Ly = Lx - 1u;
        float x0;
        const int job = 'E', compz = 'N';  //eigenvalues only
        const lapack_int ldh = (int)Ly, n = (int)Ly, ldz = 1;
        const lapack_int ilo = 1, ihi = n;  //avoids balancing
        lapack_int info;
        float *compan, *wr, *wi, zz[1];
        if (!(wr=(float *)malloc((size_t)(n)*sizeof(float)))) { fprintf(stderr,"error in poly2roots_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(wi=(float *)malloc((size_t)(n)*sizeof(float)))) { fprintf(stderr,"error in poly2roots_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(compan=(float *)calloc((size_t)(n*n),sizeof(float)))) { fprintf(stderr,"error in poly2roots_s: problem with calloc. "); perror("calloc"); return 1; }

        if (Lx==N)
        {
            ++compan;
            for (size_t l=0u; l<Ly; ++l, compan+=Lx) { *compan = 1.0f; }
            compan -= Lx*Ly + 1u;
            x0 = -*X++;
            for (size_t l=0u; l<Ly; ++l, ++X, compan+=Ly) { *compan = *X / x0; }
            compan -= Ly*Ly;
            info = LAPACKE_shseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
            if (info) { fprintf(stderr,"error in poly2roots_s: lapacke decomposition failed\n"); return 1; }
            for (size_t l=0u; l<Ly; ++l) { *Y++ = *wr++; *Y++ = *wi++; }
            wr -= Ly; wi -= Ly;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v)
                {
                    for (size_t l=0u; l<Ly*Ly; ++l, ++compan) { *compan = 0.0f; }
                    compan -= Lx;
                    for (size_t l=2u; l<Ly; ++l, compan-=Lx) { *compan = 1.0f; }
                    *compan-- = 1.0f;
                    x0 = -*X++;
                    for (size_t l=0u; l<Ly; ++l, ++X, compan+=Ly) { *compan = *X / x0; }
                    compan -= Ly*Ly;
                    info = LAPACKE_shseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
                    if (info) { fprintf(stderr,"error in poly2roots_s: lapacke decomposition failed\n"); return 1; }
                    for (size_t l=0u; l<Ly; ++l) { *Y++ = *wr++; *Y++ = *wi++; }
                    wr -= Ly; wi -= Ly;
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=K*Lx-1u, Y-=2u*K*Ly-2u)
                    {
                        for (size_t l=0u; l<Ly*Ly; ++l, ++compan) { *compan = 0.0f; }
                        compan -= Lx;
                        for (size_t l=2u; l<Ly; ++l, compan-=Lx) { *compan = 1.0f; }
                        *compan-- = 1.0f;
                        x0 = -*X; X += K;
                        for (size_t l=0u; l<Ly; ++l, X+=K, compan+=Ly) { *compan = *X / x0; }
                        compan -= Ly*Ly;
                        info = LAPACKE_shseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
                        if (info) { fprintf(stderr,"error in poly2roots_s: lapacke decomposition failed\n"); return 1; }
                        for (size_t l=0u; l<Ly; ++l, Y+=2u*K-1u) { *Y = *wr++; *++Y = *wi++; }
                        wr -= Ly; wi -= Ly;
                    }
                }
            }
        }

        free(compan); free(wr); free(wi);
    }

    return 0;
}


int poly2roots_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in poly2roots_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u || Lx<2u) {}
    else if (Lx==2u)
    {
        if (Lx==N) { *Y = -*(X+1); *++Y = 0.0; }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, ++X, ++Y) { ++X; *Y = -*X; *++Y = 0.0; }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=B)
                {
                    for (size_t b=0u; b<B; ++b, ++X, ++Y) { *Y = -*(X+K); *++Y = 0.0; }
                }
            }
        }
    }
    else
    {
        const size_t Ly = Lx - 1u;
        double x0;

        const int job = 'E', compz = 'N';  //eigenvalues only
        const lapack_int ldh = (int)Ly, n = (int)Ly, ldz = 1;
        const lapack_int ilo = 1, ihi = n;  //avoids balancing
        lapack_int info;
        double *compan, *wr, *wi, zz[1];
        if (!(wr=(double *)malloc((size_t)(n)*sizeof(double)))) { fprintf(stderr,"error in poly2roots_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(wi=(double *)malloc((size_t)(n)*sizeof(double)))) { fprintf(stderr,"error in poly2roots_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(compan=(double *)calloc((size_t)(n*n),sizeof(double)))) { fprintf(stderr,"error in poly2roots_d: problem with calloc. "); perror("calloc"); return 1; }

        if (Lx==N)
        {
            ++compan;
            for (size_t l=0u; l<Ly; ++l, compan+=Lx) { *compan = 1.0; }
            compan -= Lx*Ly + 1u;
            x0 = -*X++;
            for (size_t l=0u; l<Ly; ++l, ++X, compan+=Ly) { *compan = *X / x0; }
            compan -= Ly*Ly;
            info = LAPACKE_dhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
            if (info) { fprintf(stderr,"error in poly2roots_d: lapacke decomposition failed\n"); return 1; }
            for (size_t l=0u; l<Ly; ++l) { *Y++ = *wr++; *Y++ = *wi++; }
            wr -= Ly; wi -= Ly;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v)
                {
                    for (size_t l=0u; l<Ly*Ly; ++l, ++compan) { *compan = 0.0; }
                    compan -= Lx;
                    for (size_t l=2u; l<Ly; ++l, compan-=Lx) { *compan = 1.0; }
                    *compan-- = 1.0;
                    x0 = -*X++;
                    for (size_t l=0u; l<Ly; ++l, ++X, compan+=Ly) { *compan = *X / x0; }
                    compan -= Ly*Ly;
                    info = LAPACKE_dhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
                    if (info) { fprintf(stderr,"error in poly2roots_d: lapacke decomposition failed\n"); return 1; }
                    for (size_t l=0u; l<Ly; ++l) { *Y++ = *wr++; *Y++ = *wi++; }
                    wr -= Ly; wi -= Ly;
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=K*Lx-1u, Y-=2u*K*Ly-2u)
                    {
                        for (size_t l=0u; l<Ly*Ly; ++l, ++compan) { *compan = 0.0; }
                        compan -= Lx;
                        for (size_t l=2u; l<Ly; ++l, compan-=Lx) { *compan = 1.0; }
                        *compan-- = 1.0;
                        x0 = -*X; X += K;
                        for (size_t l=0u; l<Ly; ++l, X+=K, compan+=Ly) { *compan = *X / x0; }
                        compan -= Ly*Ly;
                        info = LAPACKE_dhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
                        if (info) { fprintf(stderr,"error in poly2roots_d: lapacke decomposition failed\n"); return 1; }
                        for (size_t l=0u; l<Ly; ++l, Y+=2u*K-1u) { *Y = *wr++; *++Y = *wi++; }
                        wr -= Ly; wi -= Ly;
                    }
                }
            }
        }

        free(compan); free(wr); free(wi);
    }

    return 0;
}


int poly2roots_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in poly2roots_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u || Lx<2u) {}
    else
    {
        const size_t Ly = Lx - 1u;
        const int job = 'E', compz = 'N';  //eigenvalues only
        const lapack_int ldh = (int)Ly, n = (int)Ly, ldz = 1;
        const lapack_int ilo = 1, ihi = n;  //avoids balancing
        lapack_int info;
        float *compan, zz[2], scr, sci;
        if (!(compan=(float *)calloc((size_t)(4*n*n),sizeof(float)))) { fprintf(stderr,"error in poly2roots_c: problem with calloc. "); perror("calloc"); return 1; }

        if (Lx==N)
        {
            compan += 2;
            for (size_t l=0u; l<Ly; ++l, compan+=2u*Lx) { *compan = 1.0f; }
            compan -= 2u*Lx*Ly + 2u;
            scr = 1.0f/(*X**X+*(X+1)**(X+1)); sci = *(X+1)*scr; scr *= -*X;
            X += 2;
            for (size_t l=0u; l<Ly; ++l, X+=2, compan+=2u*Ly-1u)
            {
                *compan++ = *X*scr - *(X+1)*sci;
                *compan = *X*sci + *(X+1)*scr;
            }
            compan -= 2u*Ly*Ly;
            info = LAPACKE_chseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_float *)compan,ldh,(lapack_complex_float *)Y,(lapack_complex_float *)zz,ldz);
            if (info) { fprintf(stderr,"error in poly2roots_c: lapacke decomposition failed\n"); return 1; }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, Y+=2u*Ly)
                {
                    if (Ly>1u)
                    {
                        for (size_t l=0u; l<2u*Ly*Ly; ++l, ++compan) { *compan = 0.0f; }
                        compan -= 2u*Lx;
                        for (size_t l=2u; l<Ly; ++l, compan-=2u*Lx) { *compan = 1.0f; }
                        *compan = 1.0f; compan -= 2;
                    }
                    scr = 1.0f/(*X**X+*(X+1)**(X+1)); sci = *(X+1)*scr; scr *= -*X;
                    X += 2;
                    for (size_t l=0u; l<Ly; ++l, X+=2, compan+=2u*Ly-1u)
                    {
                        *compan++ = *X*scr - *(X+1)*sci;
                        *compan = *X*sci + *(X+1)*scr;
                    }
                    compan -= 2u*Ly*Ly;
                    info = LAPACKE_chseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_float *)compan,ldh,(lapack_complex_float *)Y,(lapack_complex_float *)zz,ldz);
                    if (info) { fprintf(stderr,"error in poly2roots_c: lapacke decomposition failed\n"); return 1; }
                }
            }
            else
            {
                float *roots;
                if (!(roots=(float *)malloc((size_t)(2*n)*sizeof(float)))) { fprintf(stderr,"error in poly2roots_c: problem with malloc. "); perror("malloc"); return 1; }
                
                for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=2u*K*Lx-2u, Y-=2u*K*Ly-2u)
                    {
                        if (Ly>1u)
                        {
                            for (size_t l=0u; l<2u*Ly*Ly; ++l, ++compan) { *compan = 0.0f; }
                            compan -= 2u*Lx;
                            for (size_t l=2u; l<Ly; ++l, compan-=2u*Lx) { *compan = 1.0f; }
                            *compan = 1.0f; compan -= 2;
                        }
                        scr = 1.0f/(*X**X+*(X+1)**(X+1)); sci = *(X+1)*scr; scr *= -*X;
                        X += 2u*K;
                        for (size_t l=0u; l<Ly; ++l, X+=2u*K, compan+=2u*Ly-1u)
                        {
                            *compan++ = *X*scr - *(X+1)*sci;
                            *compan = *X*sci + *(X+1)*scr;
                        }
                        compan -= 2u*Ly*Ly;
                        info = LAPACKE_chseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_float *)compan,ldh,(lapack_complex_float *)roots,(lapack_complex_float *)zz,ldz);
                        if (info) { fprintf(stderr,"error in poly2roots_c: lapacke decomposition failed\n"); return 1; }
                        for (size_t l=0u; l<Ly; ++l, ++roots, Y+=2u*K-1u) { *Y++ = *roots++; *Y = *roots; }
                        roots -= 2u*Ly;
                    }
                }
                free(roots);
            }
        }
        free(compan);
    }

    return 0;
}


int poly2roots_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in poly2roots_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u || Lx<2u) {}
    else
    {
        const size_t Ly = Lx - 1u;
        const int job = 'E', compz = 'N';  //eigenvalues only
        const lapack_int ldh = (int)Ly, n = (int)Ly, ldz = 1;
        const lapack_int ilo = 1, ihi = n;  //avoids balancing
        lapack_int info;
        double *compan, zz[2], scr, sci;
        if (!(compan=(double *)calloc((size_t)(4*n*n),sizeof(double)))) { fprintf(stderr,"error in poly2roots_z: problem with calloc. "); perror("calloc"); return 1; }

        if (Lx==N)
        {
            compan += 2;
            for (size_t l=0u; l<Ly; ++l, compan+=2u*Lx) { *compan = 1.0; }
            compan -= 2u*Lx*Ly + 2u;
            scr = 1.0/(*X**X+*(X+1)**(X+1)); sci = *(X+1)*scr; scr *= -*X;
            X += 2;
            for (size_t l=0u; l<Ly; ++l, X+=2, compan+=2u*Ly-1u)
            {
                *compan++ = *X*scr - *(X+1)*sci;
                *compan = *X*sci + *(X+1)*scr;
            }
            compan -= 2u*Ly*Ly;
            info = LAPACKE_zhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_double *)compan,ldh,(lapack_complex_double *)Y,(lapack_complex_double *)zz,ldz);
            if (info) { fprintf(stderr,"error in poly2roots_z: lapacke decomposition failed\n"); return 1; }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, Y+=2u*Ly)
                {
                    if (Ly>1u)
                    {
                        for (size_t l=0u; l<2u*Ly*Ly; ++l, ++compan) { *compan = 0.0; }
                        compan -= 2u*Lx;
                        for (size_t l=2u; l<Ly; ++l, compan-=2u*Lx) { *compan = 1.0; }
                        *compan = 1.0; compan -= 2;
                    }
                    scr = 1.0/(*X**X+*(X+1)**(X+1)); sci = *(X+1)*scr; scr *= -*X;
                    X += 2;
                    for (size_t l=0u; l<Ly; ++l, X+=2, compan+=2u*Ly-1u)
                    {
                        *compan++ = *X*scr - *(X+1)*sci;
                        *compan = *X*sci + *(X+1)*scr;
                    }
                    compan -= 2u*Ly*Ly;
                    info = LAPACKE_zhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_double *)compan,ldh,(lapack_complex_double *)Y,(lapack_complex_double *)zz,ldz);
                    if (info) { fprintf(stderr,"error in poly2roots_z: lapacke decomposition failed\n"); return 1; }
                }
            }
            else
            {
                double *roots;
                if (!(roots=(double *)malloc((size_t)(2*n)*sizeof(double)))) { fprintf(stderr,"error in poly2roots_z: problem with malloc. "); perror("malloc"); return 1; }
                
                for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=2u*K*Lx-2u, Y-=2u*K*Ly-2u)
                    {
                        if (Ly>1u)
                        {
                            for (size_t l=0u; l<2u*Ly*Ly; ++l, ++compan) { *compan = 0.0; }
                            compan -= 2u*Lx;
                            for (size_t l=2u; l<Ly; ++l, compan-=2u*Lx) { *compan = 1.0; }
                            *compan = 1.0; compan -= 2;
                        }
                        scr = 1.0/(*X**X+*(X+1)**(X+1)); sci = *(X+1)*scr; scr *= -*X;
                        X += 2u*K;
                        for (size_t l=0u; l<Ly; ++l, X+=2u*K, compan+=2u*Ly-1u)
                        {
                            *compan++ = *X*scr - *(X+1)*sci;
                            *compan = *X*sci + *(X+1)*scr;
                        }
                        compan -= 2u*Ly*Ly;
                        info = LAPACKE_zhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_double *)compan,ldh,(lapack_complex_double *)roots,(lapack_complex_double *)zz,ldz);
                        if (info) { fprintf(stderr,"error in poly2roots_z: lapacke decomposition failed\n"); return 1; }
                        for (size_t l=0u; l<Ly; ++l, ++roots, Y+=2u*K-1u) { *Y++ = *roots++; *Y = *roots; }
                        roots -= 2u*Ly;
                    }
                }
                free(roots);
            }
        }
        free(compan);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
