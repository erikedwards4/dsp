//Teager-Kaiser Energy Operator (TKEO) according to de Matos [2018],
//which is computed using outputs of smooth_diff and smooth_diffdiff:
//y[n] = x'[n]*x'[n] - x[n]*x''[n]

//For complex input, the output is real, and is equal to the sum
//of the TKEOs of the real and imaginary parts, as proven
//by Hamila et al. [1999], as cited in de Matos [2018].

//See de Matos MC. 2018. Seismic attributes from the complex Teager-Kaiser energy.
//and its excellent reference:
//Holoborodko P. 2008. Smooth noise robust differentiators. www.holoborodko.com.

//For real X, the exact method of de Matos [2018] would be:
//analytic_sig X | tkeo_smooth > Y\n";
//where Y is the "complex TK energy" (but is real-valued).
//(The input to tkeo_smooth will be complex, so sum rule above applies.)

//This paper also gives the complex TKV (variational Teager-Kaiser) method,
//which is a simple extension of this, but would require looking at a few of his references.
//Basically, a further Hilbert transform of the (real-valued) "complex TK energy" gives the TKV.

//I decided to use signal extrapolation by zero-padding,
//really just to use xcorr1 to do the hard work.
//Since using xcorr1, the framing is controlled by cs0, stride, Ly.

#include <stdio.h>
#include <stdlib.h>
#include "xcorr1.c"
#include "codee_dsp.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int tkeo_smooth_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int cs0, const size_t str, const size_t dil, const size_t Ly)
{
    if (dim>3u) { fprintf(stderr,"error in tkeo_smooth_s: dim must be in [0 3]\n"); return 1; }
    if (str<1u) { fprintf(stderr,"error in tkeo_smooth_s: str must be positive\n"); return 1; }
    if (Ly<1u) { fprintf(stderr,"error in tkeo_smooth_s: Ly (length of vecs in Y) must be positive\n"); return 1; }

    int ret = 0;
    const size_t Nx = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ny = Ly*(Nx/Lx);
    if (Lx<11u) { fprintf(stderr,"error in tkeo_smooth_s: Lx (length of vecs in X) must be >=11 (length of smooth_diff)\n"); return 1; }

    //Init 1st-order differentiator (n=2, N=11 of smooth_diff)
    const float D[11] = {-1.0f/512.0f,-8.0f/512.0f,-27.0f/512.0f,-48.0f/512.0f,-42.0f/512.0,0.0f,42.0f/512.0f,48.0f/512.0f,27.0/512.0f,8.0f/512.0f,1.0f/512.0f};

    //Init 2nd-order differentiator (n=5, N=9 of smooth_diffdiff)
    const float DD[9] = {-7.0f/192.0f,12.0f/192.0f,52.0f/192.0f,-12.0f/192.0f,-90.0f/192.0f,-12.0f/192.0f,52.0f/192.0f,12.0f/192.0f,-7.0f/192.0f};

    //Init Xd for intermediate output of diff and diffdiff
    float *Xd = (float *)malloc(Ny*sizeof(float));
    if (!Xd) { fprintf(stderr,"error in tkeo_smooth_s: problem with malloc. "); perror("malloc"); return 1; }

    //Smooth diff of X
    ret = xcorr1_s (Xd,X,D,R,C,S,H,iscolmajor,11u,cs0+5*(int)dil,str,dil,Ly,dim);
    //fprintf(stderr,"Xd = "); for (size_t n=0u; n<Ny; ++n) { fprintf(stderr,"%g ",(double)Xd[n]); } fprintf(stderr,"\n");

    //Get Xd^2 (1st term of smooth TKEO)
    for (size_t n=Ny; n>0u; --n, ++Xd, ++Y) { *Y = *Xd * *Xd; }
    Xd -= Ny; Y -= Ny;

    //Smooth diffdiff of X
    ret = xcorr1_s (Xd,X,DD,R,C,S,H,iscolmajor,9u,cs0+4*(int)dil,str,dil,Ly,dim);
    //fprintf(stderr,"Xdd = "); for (size_t n=0u; n<Ny; ++n) { fprintf(stderr,"%g ",(double)Xd[n]); } fprintf(stderr,"\n");

    //Get smooth TKEO
    X += cs0;
    for (size_t n=Ny; n>0u; --n, X+=str, ++Xd, ++Y) { *Y -= *X * *Xd; }
    Xd -= Ny;

    //Finish
    free(Xd);
    return ret;
}


int tkeo_smooth_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int cs0, const size_t str, const size_t dil, const size_t Ly)
{
    if (dim>3u) { fprintf(stderr,"error in tkeo_smooth_d: dim must be in [0 3]\n"); return 1; }
    if (str<1u) { fprintf(stderr,"error in tkeo_smooth_d: str must be positive\n"); return 1; }
    if (Ly<1u) { fprintf(stderr,"error in tkeo_smooth_d: Ly (length of vecs in Y) must be positive\n"); return 1; }

    int ret = 0;
    const size_t Nx = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ny = Ly*(Nx/Lx);
    if (Lx<11u) { fprintf(stderr,"error in tkeo_smooth_d: L (length of vecs in X) must be >=11 (length of smooth_diff)\n"); return 1; }

    //Init 1st-order differentiator (n=2, N=11 of smooth_diff)
    const double D[11] = {-1.0/512.0,-8.0/512.0,-27.0/512.0,-48.0/512.0,-42.0/512.0,0.0,42.0/512.0,48.0/512.0,27.0/512.0,8.0/512.0,1.0/512.0};

    //Init 2nd-order differentiator (n=5, N=9 of smooth_diffdiff)
    const double DD[9] = {-7.0/192.0,12.0/192.0,52.0/192.0,-12.0/192.0,-90.0/192.0,-12.0/192.0,52.0/192.0,12.0/192.0,-7.0/192.0};

    //Init Xd for intermediate output of diff and diffdiff
    double *Xd = (double *)malloc(Ny*sizeof(double));
    if (!Xd) { fprintf(stderr,"error in tkeo_smooth_d: problem with malloc. "); perror("malloc"); return 1; }

    //Smooth diff of X
    ret = xcorr1_d (Xd,X,D,R,C,S,H,iscolmajor,11u,cs0+5*(int)dil,str,dil,Ly,dim);

    //Get Xd^2 (1st term of smooth TKEO)
    for (size_t n=Ny; n>0u; --n, ++Xd, ++Y) { *Y = *Xd * *Xd; }
    Xd -= Ny; Y -= Ny;

    //Smooth diffdiff of X
    ret = xcorr1_d (Xd,X,DD,R,C,S,H,iscolmajor,9u,cs0+4*(int)dil,str,dil,Ly,dim);

    //Get smooth TKEO
    for (size_t n=Ny; n>0u; --n, ++X, ++Xd, ++Y) { *Y -= *X * *Xd; }
    Xd -= Ny;

    //Finish
    free(Xd);
    return ret;
}


int tkeo_smooth_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int cs0, const size_t str, const size_t dil, const size_t Ly)
{
    if (dim>3u) { fprintf(stderr,"error in tkeo_smooth_c: dim must be in [0 3]\n"); return 1; }
    if (str<1u) { fprintf(stderr,"error in tkeo_smooth_c: str must be positive\n"); return 1; }
    if (Ly<1u) { fprintf(stderr,"error in tkeo_smooth_c: Ly (length of vecs in Y) must be positive\n"); return 1; }

    int ret = 0;
    const size_t Nx = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ny = Ly*(Nx/Lx);
    if (Lx<11u) { fprintf(stderr,"error in tkeo_smooth_c: Lx (length of vecs in X) must be >=11 (length of smooth_diff)\n"); return 1; }

    //These are needed for using real-valued xcorr for complex input
    const size_t R2 = (dim==0u) ? 2u*R : R;
    const size_t C2 = (dim==1u) ? 2u*C : C;
    const size_t S2 = (dim==2u) ? 2u*S : S;
    const size_t H2 = (dim==3u) ? 2u*H : H;

    //Init 1st-order differentiator (n=2, N=11 of smooth_diff)
    const float D[11] = {-1.0f/512.0f,-8.0f/512.0f,-27.0f/512.0f,-48.0f/512.0f,-42.0f/512.0,0.0f,42.0f/512.0f,48.0f/512.0f,27.0/512.0f,8.0f/512.0f,1.0f/512.0f};

    //Init 2nd-order differentiator (n=5, N=9 of smooth_diffdiff)
    const float DD[9] = {-7.0f/192.0f,12.0f/192.0f,52.0f/192.0f,-12.0f/192.0f,-90.0f/192.0f,-12.0f/192.0f,52.0f/192.0f,12.0f/192.0f,-7.0f/192.0f};

    //Init Xd for intermediate output of diff and diffdiff
    float *Xd = (float *)malloc(2u*Ny*sizeof(float));
    if (!Xd) { fprintf(stderr,"error in tkeo_smooth_c: problem with malloc. "); perror("malloc"); return 1; }

    //Smooth diff of real(X)
    ret = xcorr1_s (Xd,X,D,R2,C2,S2,H2,iscolmajor,11u,2*cs0+10*(int)dil,2u*str,2u*dil,Ly,dim);
    //fprintf(stderr,"Xdr = "); for (size_t n=0u; n<Ny; ++n) { fprintf(stderr,"%g ",(double)Xd[n]); } fprintf(stderr,"\n");

    //Get Xd^2 for real part
    for (size_t n=Ny; n>0u; --n, ++Xd, ++Y) { *Y = *Xd * *Xd; }
    Xd -= Ny; Y -= Ny;

    //Smooth diffdiff of real(X)
    ret = xcorr1_s (Xd,X,DD,R2,C2,S2,H2,iscolmajor,9u,2*cs0+8*(int)dil,2u*str,2u*dil,Ly,dim);
    //fprintf(stderr,"Xddr = "); for (size_t n=0u; n<Ny; ++n) { fprintf(stderr,"%g ",(double)Xd[n]); } fprintf(stderr,"\n");

    //Get smooth TKEO for real part
    X += 2*cs0;
    for (size_t n=Ny; n>0u; --n, X+=2u*str, ++Xd, ++Y) { *Y -= *X * *Xd; }
    X -= 2*cs0 + (int)(2u*Ny*str);
    Xd -= Ny; Y -= Ny;

    //Smooth diff of imag(X)
    ret = xcorr1_s (Xd,X,D,R2,C2,S2,H2,iscolmajor,11u,2*cs0+10*(int)dil+1,2u*str,2u*dil,Ly,dim);
    //fprintf(stderr,"Xdi = "); for (size_t n=0u; n<Ny; ++n) { fprintf(stderr,"%g ",(double)Xd[n]); } fprintf(stderr,"\n");

    //Add in Xd^2 for imag part
    for (size_t n=Ny; n>0u; --n, ++Xd, ++Y) { *Y += *Xd * *Xd; }
    Xd -= Ny; Y -= Ny;

    //Smooth diffdiff of imag(X)
    ret = xcorr1_s (Xd,X,DD,R2,C2,S2,H2,iscolmajor,9u,2*cs0+8*(int)dil+1,2u*str,2u*dil,Ly,dim);
    //fprintf(stderr,"Xddi = "); for (size_t n=0u; n<Ny; ++n) { fprintf(stderr,"%g ",(double)Xd[n]); } fprintf(stderr,"\n");

    //Get smooth complex TKEO
    X += 2*cs0 + 1;
    for (size_t n=Ny; n>0u; --n, X+=2u*str, ++Xd, ++Y) { *Y -= *X * *Xd; }
    X -= 2*cs0 + 1 + (int)(2u*Ny*str);
    Xd -= Ny; Y -= Ny;

    //Finish
    free(Xd);
    return ret;
}


int tkeo_smooth_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int cs0, const size_t str, const size_t dil, const size_t Ly)
{
    if (dim>3u) { fprintf(stderr,"error in tkeo_smooth_z: dim must be in [0 3]\n"); return 1; }
    if (str<1u) { fprintf(stderr,"error in tkeo_smooth_z: str must be positive\n"); return 1; }
    if (Ly<1u) { fprintf(stderr,"error in tkeo_smooth_z: Ly (length of vecs in Y) must be positive\n"); return 1; }

    int ret = 0;
    const size_t Nx = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ny = Ly*(Nx/Lx);
    if (Lx<11u) { fprintf(stderr,"error in tkeo_smooth_z: Lx (length of vecs in X) must be >=11 (length of smooth_diff)\n"); return 1; }

    //These are needed for using real-valued xcorr for complex input
    const size_t R2 = (dim==0u) ? 2u*R : R;
    const size_t C2 = (dim==1u) ? 2u*C : C;
    const size_t S2 = (dim==2u) ? 2u*S : S;
    const size_t H2 = (dim==3u) ? 2u*H : H;

    //Init 1st-order differentiator (n=2, N=11 of smooth_diff)
    const double D[11] = {-1.0/512.0,-8.0/512.0,-27.0/512.0,-48.0/512.0,-42.0/512.0,0.0,42.0/512.0,48.0/512.0,27.0/512.0,8.0/512.0,1.0/512.0};

    //Init 2nd-order differentiator (n=5, N=9 of smooth_diffdiff)
    const double DD[9] = {-7.0/192.0,12.0/192.0,52.0/192.0,-12.0/192.0,-90.0/192.0,-12.0/192.0,52.0/192.0,12.0/192.0,-7.0/192.0};

    //Init Xd for intermediate output of diff and diffdiff
    double *Xd = (double *)malloc(2u*Ny*sizeof(double));
    if (!Xd) { fprintf(stderr,"error in tkeo_smooth_z: problem with malloc. "); perror("malloc"); return 1; }

    //Smooth diff of real(X)
    ret = xcorr1_d (Xd,X,D,R2,C2,S2,H2,iscolmajor,11u,2*cs0+10*(int)dil,2u*str,2u*dil,Ly,dim);

    //Get Xd^2 for real part
    for (size_t n=Ny; n>0u; --n, ++Xd, ++Y) { *Y = *Xd * *Xd; }
    Xd -= Ny; Y -= Ny;

    //Smooth diffdiff of real(X)
    ret = xcorr1_d (Xd,X,DD,R2,C2,S2,H2,iscolmajor,9u,2*cs0+8*(int)dil,2u*str,2u*dil,Ly,dim);

    //Get smooth TKEO for real part
    X += 2*cs0;
    for (size_t n=Ny; n>0u; --n, X+=2u*str, ++Xd, ++Y) { *Y -= *X * *Xd; }
    X -= 2*cs0 + (int)(2u*Ny*str);
    Xd -= Ny; Y -= Ny;

    //Smooth diff of imag(X)
    ret = xcorr1_d (Xd,X,D,R2,C2,S2,H2,iscolmajor,11u,2*cs0+10*(int)dil+1,2u*str,2u*dil,Ly,dim);

    //Add in Xd^2 for imag part
    for (size_t n=Ny; n>0u; --n, ++Xd, ++Y) { *Y += *Xd * *Xd; }
    Xd -= Ny; Y -= Ny;

    //Smooth diffdiff of imag(X)
    ret = xcorr1_d (Xd,X,DD,R2,C2,S2,H2,iscolmajor,9u,2*cs0+8*(int)dil+1,2u*str,2u*dil,Ly,dim);

    //Get smooth complex TKEO
    X += 2*cs0 + 1;
    for (size_t n=Ny; n>0u; --n, X+=2u*str, ++Xd, ++Y) { *Y -= *X * *Xd; }
    X -= 2*cs0 + 1 + (int)(2u*Ny*str);
    Xd -= Ny; Y -= Ny;

    //Finish
    free(Xd);
    return ret;
}


#ifdef __cplusplus
}
}
#endif
