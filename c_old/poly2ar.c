//Gets autoregressive (AR) parameters from polynomials along rows or cols of X.
//If the polynomial is a0 a1 a2..., then the AR coeffs are -a1/a0 -a2/a0...
//This is invertible with ar2poly only if a0==1.

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int poly2ar_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int poly2ar_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
//int poly2ar_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
//int poly2ar_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int poly2ar_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in poly2ar_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t Ly = Lx - 1;

    if (N==Lx)
    {
        const float xi = -1.0f / X[0];
        for (size_t l=0; l<Ly; l++) { Y[l] = xi * X[l+1];  }
    }
    else
    {
        const size_t B = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? H*S*C : (dim==1) ? H*S : (dim==2) ? H : 1);
        const size_t G = N/(B*Lx);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? H*S*C : (dim==1) ? H*S : (dim==2) ? H : 1);
        const size_t Jx = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        const size_t Jy = (iscolmajor) ? ((dim==0) ? R-1 : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H-1);
        for (size_t g=0; g<G; g++, X+=B*(Lx-Jx), Y+=B*(Ly-Jy))
        {
            for (size_t b=0; b<B; b++, X+=Jx, Y+=Jy)
            {
                cblas_scopy((int)Ly,&X[K],(int)K,Y,(int)K);
                cblas_sscal((int)Ly,-1.0f/X[0],Y,(int)K);
            }
        }
    }

    return 0;
}


int poly2ar_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in poly2ar_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t Ly = Lx - 1;

    if (N==Lx)
    {
        const double xi = -1.0 / X[0];
        for (size_t l=0; l<Ly; l++) { Y[l] = xi * X[l+1];  }
    }
    else
    {
        const size_t B = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? H*S*C : (dim==1) ? H*S : (dim==2) ? H : 1);
        const size_t G = N/(B*Lx);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? H*S*C : (dim==1) ? H*S : (dim==2) ? H : 1);
        const size_t Jx = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        const size_t Jy = (iscolmajor) ? ((dim==0) ? R-1 : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H-1);
        for (size_t g=0; g<G; g++, X+=B*(Lx-Jx), Y+=B*(Ly-Jy))
        {
            for (size_t b=0; b<B; b++, X+=Jx, Y+=Jy)
            {
                cblas_dcopy((int)Ly,&X[K],(int)K,Y,(int)K);
                cblas_dscal((int)Ly,-1.0/X[0],Y,(int)K);
            }
        }
    }

    return 0;
}


int poly2ar_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in poly2ar_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t Ly = Lx - 1;

    float sc[2] = {0.0f,0.0f};

    if (N==Lx)
    {
        const double xi = -1.0 / X[0];
        for (size_t l=0; l<Ly; l++) { Y[l] = xi * X[l+1];  }
    }
    else
    {
        const size_t B = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? H*S*C : (dim==1) ? H*S : (dim==2) ? H : 1);
        const size_t G = N/(B*Lx);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? H*S*C : (dim==1) ? H*S : (dim==2) ? H : 1);
        const size_t Jx = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        const size_t Jy = (iscolmajor) ? ((dim==0) ? R-1 : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H-1);
        for (size_t g=0; g<G; g++, X+=B*(Lx-Jx), Y+=B*(Ly-Jy))
        {
            for (size_t b=0; b<B; b++, X+=Jx, Y+=Jy)
            {
                cblas_ccopy((int)Ly,&X[K],(int)K,Y,(int)K);
                cblas_cscal((int)Ly,&sc[0],Y,2*(int)K);
            }
        }
    }

    return 0;
}


// int poly2ar_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
// {
//     if (dim>3) { fprintf(stderr,"error in poly2ar_z: dim must be in [0 3]\n"); return 1; }

//     double sc[2] = {0.0,0.0};

//     if (dim==0)
//     {
//         if (iscolmajor)
//         {
//             for (size_t c=0; c<C; c++)
//             {
//                 cblas_zcopy((int)R-1,&X[2*c*R+2],1,&Y[2*c*(R-1)],1);
//                 sc[0] = 1.0 / (X[2*c*R]*X[2*c*R]+X[2*c*R+1]*X[2*c*R+1]);
//                 sc[1] = sc[0] * X[2*c*R+1]; sc[0] *= -X[2*c*R];
//                 cblas_zscal((int)R-1,&sc[0],&Y[2*c*(R-1)],1);
//             }
//         }
//         else
//         {
//             cblas_zcopy((int)((R-1)*C),&X[2*C],1,Y,1);
//             for (size_t c=0; c<C; c++)
//             {
//                 sc[0] = 1.0 / (X[2*c]*X[2*c]+X[2*c+1]*X[2*c+1]);
//                 sc[1] = sc[0] * X[2*c+1]; sc[0] *= -X[2*c];
//                 cblas_zscal((int)R-1,&sc[0],&Y[2*c],2*(int)C);
//             }
//         }
//     }
//     else if (dim==1)
//     {
//         if (iscolmajor)
//         {
//             cblas_zcopy((int)(R*(C-1)),&X[2*R],1,Y,1);
//             for (size_t r=0; r<R; r++)
//             {
//                 sc[0] = 1.0 / (X[2*r]*X[2*r]+X[2*r+1]*X[2*r+1]);
//                 sc[1] = sc[0] * X[2*r+1]; sc[0] *= -X[2*r];
//                 cblas_zscal((int)C-1,&sc[0],&Y[2*r],(int)R);
//             }
//         }
//         else
//         {
//             for (size_t r=0; r<R; r++)
//             {
//                 cblas_zcopy((int)C-1,&X[2*r*C+2],1,&Y[2*r*(C-1)],1);
//                 sc[0] = 1.0 / (X[2*r*C]*X[2*r*C]+X[2*r*C+1]*X[2*r*C+1]);
//                 sc[1] = sc[0] * X[2*r*C+1]; sc[0] *= -X[2*r*C];
//                 cblas_zscal((int)C-1,&sc[0],&Y[2*r*(C-1)],1);
//             }
//         }
//     }
//     else
//     {
//         fprintf(stderr,"error in poly2ar_z: dim must be 0 or 1.\n"); return 1;
//     }

//     //Exit
//     return 0;
// }


#ifdef __cplusplus
}
}
#endif
