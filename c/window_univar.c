//This takes a uniivariate time series and outputs a set of windowed frames.
//The step size (stp) is allowed to be non-integer,
//but the center sample of each frame is always rounded to the nearest integer.

//This istp==stp version does not work frame-by-frame; it is much faster to go length-wise.
//However, the step size (stp) must be integer, so that is tested for here.

//Within the istp==stp version, the cblas_saxpy version is about 2x slower.
//On the other hand, the while loop version is just as fast (see commented code below for dim==1 and iscolmajor).

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int window_univar_s (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const float *W, const int dim, const int c0, const float stp, const char mn0);
int window_univar_d (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const double *W, const int dim, const int c0, const double stp, const char mn0);
int window_univar_c (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const float *W, const int dim, const int c0, const float stp, const char mn0);
int window_univar_z (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const double *W, const int dim, const int c0, const double stp, const char mn0);


int window_univar_s (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const float *W, const int dim, const int c0, const float stp, const char mn0)
{
    const float z = 0.0f, o = 1.0f;
    const int T = (dim==0u) ? C : R;
    const int L = (dim==0u) ? R : C;
    const int Lpre = L/2;      //nsamps before center samp
    const int istp = (int)stp;
    int ss = c0 - Lpre;        //start samp of current frame
    int l, t = 0;
    int Tpre, Tmid, Tpost;
    int ssl = ss, esl = ss + (T-1)*istp;

    //Checks
    if (N<1) { fprintf(stderr,"error in window_univar_s: N (length X) must be positive\n"); return 1; }
    if (R<1) { fprintf(stderr,"error in window_univar_s: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in window_univar_s: C (ncols Y) must be positive\n"); return 1; }
    if (c0<0) { fprintf(stderr,"error in window_univar_s: c0 (center samp of 1st frame) must be nonnegative\n"); return 1; }
    if (c0>N-1) { fprintf(stderr,"error in window_univar_s: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in window_univar_s: L (winlength) must be positive\n"); return 1; }
    if (L>=N) { fprintf(stderr,"error in window_univar_s: L (winlength) must be < N (length X)\n"); return 1; }
    if (stp<=0.0f) { fprintf(stderr,"error in window_univar_s: stp (step size) must be positive\n"); return 1; }

    if (dim==0u)
    {
        if (iscolmajor)
        {
            if (stp==floorf(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_scopy(Tpre,&z,0,&Y[l],L); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_scopy(Tpost,&z,0,&Y[l+L*(Tpre+Tmid)],L); }
                    if (Tmid>0) { cblas_scopy(Tmid,&X[ssl+Tpre*istp],istp,&Y[l+L*Tpre],L); cblas_sscal(Tmid,W[l],&Y[l+L*Tpre],L); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_scopy(-ss,&z,0,&Y[t*L],1);
                    cblas_ssbmv(CblasColMajor,CblasUpper,L+ss,0,1.0f,&W[-ss],1,&X[0],1,0.0f,&Y[t*L-ss],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_ssbmv(CblasColMajor,CblasUpper,L,0,1.0f,&W[0],1,&X[ss],1,0.0f,&Y[t*L],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_ssbmv(CblasColMajor,CblasUpper,N-ss,0,1.0f,&W[0],1,&X[ss],1,0.0f,&Y[t*L],1);
                    cblas_scopy(L-N+ss,&z,0,&Y[t*L+N-ss],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0) { for (t=0; t<T; t++) { cblas_saxpy(L,-cblas_sdot(L,&Y[t*L],1,&o,0)/L,&o,0,&Y[t*L],1); } }
        }
        else
        {
            if (stp==floorf(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_scopy(Tpre,&z,0,&Y[l*T],1); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_scopy(Tpost,&z,0,&Y[l*T+Tpre+Tmid],1); }
                    if (Tmid>0) { cblas_scopy(Tmid,&X[ssl+Tpre*istp],istp,&Y[l*T+Tpre],1); cblas_sscal(Tmid,W[l],&Y[l*T+Tpre],1); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_scopy(-ss,&z,0,&Y[t],T);
                    cblas_ssbmv(CblasRowMajor,CblasUpper,L+ss,0,1.0f,&W[-ss],1,&X[0],1,0.0f,&Y[-T*ss],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_ssbmv(CblasRowMajor,CblasUpper,L,0,1.0f,&W[0],1,&X[ss],1,0.0f,&Y[t],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_ssbmv(CblasRowMajor,CblasUpper,N-ss,0,1.0f,&W[0],1,&X[ss],1,0.0f,&Y[t],T);
                    cblas_scopy(L-N+ss,&z,0,&Y[t+T*(N-ss)],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0) { for (t=0; t<T; t++) { cblas_saxpy(L,-cblas_sdot(L,&Y[t],T,&o,0)/L,&o,0,&Y[t],T); } }
        }
    }
    else if (dim==1u)
    {
        if (iscolmajor)
        {
            if (stp==floorf(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T; //Tpre = 1 - ssl/istp; Tpost = 1 + (esl-N)/istp;
                    //cblas_scopy(T,&z,0,&Y[l*T],1); cblas_saxpy(T-Tpre,W[l],&X[ssl+Tpre*istp],istp,&Y[l*T+Tpre],1);
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_scopy(Tpre,&z,0,&Y[l*T],1); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_scopy(Tpost,&z,0,&Y[l*T+Tpre+Tmid],1); }
                    if (Tmid>0) { cblas_scopy(Tmid,&X[ssl+Tpre*istp],istp,&Y[l*T+Tpre],1); cblas_sscal(Tmid,W[l],&Y[l*T+Tpre],1); }
                    //t = 0; n = l*T; cs = ssl;
                    //while (t<Tpre) { Y[n] = 0.0f; t++; n++; cs+=stp; }
                    //while (t<T-Tpost) { Y[n] = W[l]*X[cs]; t++; n++; cs+=stp; }
                    //while (t<T) { Y[n] = 0.0f; t++; n++; }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_scopy(-ss,&z,0,&Y[t],T);
                    cblas_ssbmv(CblasColMajor,CblasUpper,L+ss,0,1.0f,&W[-ss],1,&X[0],1,0.0f,&Y[-T*ss],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_ssbmv(CblasColMajor,CblasUpper,L,0,1.0f,&W[0],1,&X[ss],1,0.0f,&Y[t],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_ssbmv(CblasColMajor,CblasUpper,N-ss,0,1.0f,&W[0],1,&X[ss],1,0.0f,&Y[t],T);
                    cblas_scopy(L-N+ss,&z,0,&Y[t+T*(N-ss)],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0) { for (t=0; t<T; t++) { cblas_saxpy(L,-cblas_sdot(L,&Y[t],T,&o,0)/L,&o,0,&Y[t],T); } }
        }
        else
        {
            if (stp==floorf(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_scopy(Tpre,&z,0,&Y[l],L); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_scopy(Tpost,&z,0,&Y[l+L*(Tpre+Tmid)],L); }
                    if (Tmid>0) { cblas_scopy(Tmid,&X[ssl+Tpre*istp],istp,&Y[l+L*Tpre],L); cblas_sscal(Tmid,W[l],&Y[l+L*Tpre],L); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_scopy(-ss,&z,0,&Y[t*L],1);
                    cblas_ssbmv(CblasRowMajor,CblasUpper,L+ss,0,1.0f,&W[-ss],1,&X[0],1,0.0f,&Y[t*L-ss],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_ssbmv(CblasRowMajor,CblasUpper,L,0,1.0f,&W[0],1,&X[ss],1,0.0f,&Y[t*L],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_ssbmv(CblasRowMajor,CblasUpper,N-ss,0,1.0f,&W[0],1,&X[ss],1,0.0f,&Y[t*L],1);
                    cblas_scopy(L-N+ss,&z,0,&Y[t*L+N-ss],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0) { for (t=0; t<T; t++) { cblas_saxpy(L,-cblas_sdot(L,&Y[t*L],1,&o,0)/L,&o,0,&Y[t*L],1); } }
        }
    }
    else
    {
        fprintf(stderr,"error in window_univar_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int window_univar_d (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const double *W, const int dim, const int c0, const double stp, const char mn0)
{
    const double z = 0.0, o = 1.0;
    const int T = (dim==0u) ? C : R;
    const int L = (dim==0u) ? R : C;
    const int Lpre = L/2;      //nsamps before center samp
    const int istp = (int)stp;
    int ss = c0 - Lpre;        //start samp of current frame
    int l, t = 0;
    int Tpre, Tmid, Tpost;
    int ssl = ss, esl = ss + (T-1)*istp;

    //Checks
    if (N<1) { fprintf(stderr,"error in window_univar_d: N (length X) must be positive\n"); return 1; }
    if (R<1) { fprintf(stderr,"error in window_univar_d: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in window_univar_d: C (ncols Y) must be positive\n"); return 1; }
    if (c0<0) { fprintf(stderr,"error in window_univar_d: c0 (center samp of 1st frame) must be nonnegative\n"); return 1; }
    if (c0>N-1) { fprintf(stderr,"error in window_univar_d: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in window_univar_d: L (winlength) must be positive\n"); return 1; }
    if (L>=N) { fprintf(stderr,"error in window_univar_d: L (winlength) must be < N (length X)\n"); return 1; }
    if (stp<=0.0) { fprintf(stderr,"error in window_univar_d: stp (step size) must be positive\n"); return 1; }

    if (dim==0u)
    {
        if (iscolmajor)
        {
            if (stp==floor(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_dcopy(Tpre,&z,0,&Y[l],L); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_dcopy(Tpost,&z,0,&Y[l+L*(Tpre+Tmid)],L); }
                    if (Tmid>0) { cblas_dcopy(Tmid,&X[ssl+Tpre*istp],istp,&Y[l+L*Tpre],L); cblas_dscal(Tmid,W[l],&Y[l+L*Tpre],L); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_dcopy(-ss,&z,0,&Y[t*L],1);
                    cblas_dsbmv(CblasColMajor,CblasUpper,L+ss,0,1.0,&W[-ss],1,&X[0],1,0.0,&Y[t*L-ss],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_dsbmv(CblasColMajor,CblasUpper,L,0,1.0,&W[0],1,&X[ss],1,0.0,&Y[t*L],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_dsbmv(CblasColMajor,CblasUpper,N-ss,0,1.0,&W[0],1,&X[ss],1,0.0,&Y[t*L],1);
                    cblas_dcopy(L-N+ss,&z,0,&Y[t*L+N-ss],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0) { for (t=0; t<T; t++) { cblas_daxpy(L,-cblas_ddot(L,&Y[t*L],1,&o,0)/L,&o,0,&Y[t*L],1); } }
        }
        else
        {
            if (stp==floor(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_dcopy(Tpre,&z,0,&Y[l*T],1); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_dcopy(Tpost,&z,0,&Y[l*T+Tpre+Tmid],1); }
                    if (Tmid>0) { cblas_dcopy(Tmid,&X[ssl+Tpre*istp],istp,&Y[l*T+Tpre],1); cblas_dscal(Tmid,W[l],&Y[l*T+Tpre],1); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_dcopy(-ss,&z,0,&Y[t],T);
                    cblas_dsbmv(CblasRowMajor,CblasUpper,L+ss,0,1.0,&W[-ss],1,&X[0],1,0.0,&Y[-T*ss],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_dsbmv(CblasRowMajor,CblasUpper,L,0,1.0,&W[0],1,&X[ss],1,0.0,&Y[t],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_dsbmv(CblasRowMajor,CblasUpper,N-ss,0,1.0,&W[0],1,&X[ss],1,0.0,&Y[t],T);
                    cblas_dcopy(L-N+ss,&z,0,&Y[t+T*(N-ss)],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0) { for (t=0; t<T; t++) { cblas_daxpy(L,-cblas_ddot(L,&Y[t],T,&o,0)/L,&o,0,&Y[t],T); } }
        }
    }
    else if (dim==1u)
    {
        if (iscolmajor)
        {
            if (stp==floor(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_dcopy(Tpre,&z,0,&Y[l*T],1); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_dcopy(Tpost,&z,0,&Y[l*T+Tpre+Tmid],1); }
                    if (Tmid>0) { cblas_dcopy(Tmid,&X[ssl+Tpre*istp],istp,&Y[l*T+Tpre],1); cblas_dscal(Tmid,W[l],&Y[l*T+Tpre],1); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_dcopy(-ss,&z,0,&Y[t],T);
                    cblas_dsbmv(CblasColMajor,CblasUpper,L+ss,0,1.0,&W[-ss],1,&X[0],1,0.0,&Y[-T*ss],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_dsbmv(CblasColMajor,CblasUpper,L,0,1.0,&W[0],1,&X[ss],1,0.0,&Y[t],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_dsbmv(CblasColMajor,CblasUpper,N-ss,0,1.0,&W[0],1,&X[ss],1,0.0,&Y[t],T);
                    cblas_dcopy(L-N+ss,&z,0,&Y[t+T*(N-ss)],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0) { for (t=0; t<T; t++) { cblas_daxpy(L,-cblas_ddot(L,&Y[t],T,&o,0)/L,&o,0,&Y[t],T); } }
        }
        else
        {
            if (stp==floor(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_dcopy(Tpre,&z,0,&Y[l],L); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_dcopy(Tpost,&z,0,&Y[l+L*(Tpre+Tmid)],L); }
                    if (Tmid>0) { cblas_dcopy(Tmid,&X[ssl+Tpre*istp],istp,&Y[l+L*Tpre],L); cblas_dscal(Tmid,W[l],&Y[l+L*Tpre],L); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_dcopy(-ss,&z,0,&Y[t*L],1);
                    cblas_dsbmv(CblasRowMajor,CblasUpper,L+ss,0,1.0,&W[-ss],1,&X[0],1,0.0,&Y[t*L-ss],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_dsbmv(CblasRowMajor,CblasUpper,L,0,1.0,&W[0],1,&X[ss],1,0.0,&Y[t*L],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_dsbmv(CblasRowMajor,CblasUpper,N-ss,0,1.0,&W[0],1,&X[ss],1,0.0,&Y[t*L],1);
                    cblas_dcopy(L-N+ss,&z,0,&Y[t*L+N-ss],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0) { for (t=0; t<T; t++) { cblas_daxpy(L,-cblas_ddot(L,&Y[t*L],1,&o,0)/L,&o,0,&Y[t*L],1); } }
        }
    }
    else
    {
        fprintf(stderr,"error in window_univar_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int window_univar_c (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const float *W, const int dim, const int c0, const float stp, const char mn0)
{
    const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
    const int T = (dim==0u) ? C : R;
    const int L = (dim==0u) ? R : C;
    const int Lpre = L/2;      //nsamps before center samp
    const int istp = (int)stp;
    int ss = c0 - Lpre;        //start samp of current frame
    int l, t = 0;
    int Tpre, Tmid, Tpost;
    int ssl = ss, esl = ss + (T-1)*istp;

    //Checks
    if (N<1) { fprintf(stderr,"error in window_univar_c: N (length X) must be positive\n"); return 1; }
    if (R<1) { fprintf(stderr,"error in window_univar_c: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in window_univar_c: C (ncols Y) must be positive\n"); return 1; }
    if (c0<0) { fprintf(stderr,"error in window_univar_c: c0 (center samp of 1st frame) must be nonnegative\n"); return 1; }
    if (c0>N-1) { fprintf(stderr,"error in window_univar_c: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in window_univar_c: L (winlength) must be positive\n"); return 1; }
    if (L>=N) { fprintf(stderr,"error in window_univar_c: L (winlength) must be < N (length X)\n"); return 1; }
    if (stp<=0.0f) { fprintf(stderr,"error in window_univar_c: stp (step size) must be positive\n"); return 1; }

    if (dim==0u)
    {
        if (iscolmajor)
        {
            if (stp==floorf(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_ccopy(Tpre,&z[0],0,&Y[2*l],L); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_ccopy(Tpost,&z[0],0,&Y[2*(l+L*(Tpre+Tmid))],L); }
                    if (Tmid>0) { cblas_ccopy(Tmid,&X[2*(ssl+Tpre*istp)],istp,&Y[2*(l+L*Tpre)],L); cblas_cscal(Tmid,&W[2*l],&Y[2*(l+L*Tpre)],L); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_ccopy(-ss,&z[0],0,&Y[2*t*L],1);
                    cblas_chbmv(CblasColMajor,CblasUpper,L+ss,0,&o[0],&W[-2*ss],1,&X[0],1,&z[0],&Y[2*(t*L-ss)],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_chbmv(CblasColMajor,CblasUpper,L,0,&o[0],&W[0],1,&X[2*ss],1,&z[0],&Y[2*t*L],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_chbmv(CblasColMajor,CblasUpper,N-ss,0,&o[0],&W[0],1,&X[2*ss],1,&z[0],&Y[2*t*L],1);
                    cblas_ccopy(L-N+ss,&z[0],0,&Y[2*(t*L+N-ss)],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0)
            {
                for (t=0; t<T; t++)
                {
                    cblas_saxpy(L,-cblas_sdot(L,&Y[2*t*L],2,&o[0],0)/L,&o[0],0,&Y[2*t*L],2);
                    cblas_saxpy(L,-cblas_sdot(L,&Y[2*t*L+1],2,&o[0],0)/L,&o[0],0,&Y[2*t*L+1],2);
                }
            }
        }
        else
        {
            if (stp==floorf(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_ccopy(Tpre,&z[0],0,&Y[2*l*T],1); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_ccopy(Tpost,&z[0],0,&Y[2*(l*T+Tpre+Tmid)],1); }
                    if (Tmid>0) { cblas_ccopy(Tmid,&X[2*(ssl+Tpre*istp)],istp,&Y[2*(l*T+Tpre)],1); cblas_cscal(Tmid,&W[2*l],&Y[2*(l*T+Tpre)],1); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_ccopy(-ss,&z[0],0,&Y[2*t],T);
                    cblas_chbmv(CblasRowMajor,CblasUpper,L+ss,0,&o[0],&W[-2*ss],1,&X[0],1,&z[0],&Y[-2*T*ss],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_chbmv(CblasRowMajor,CblasUpper,L,0,&o[0],&W[0],1,&X[2*ss],1,&z[0],&Y[2*t],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_chbmv(CblasRowMajor,CblasUpper,N-ss,0,&o[0],&W[0],1,&X[2*ss],1,&z[0],&Y[2*t],T);
                    cblas_ccopy(L-N+ss,&z[0],0,&Y[2*(t+T*(N-ss))],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0)
            {
                for (t=0; t<T; t++)
                {
                    cblas_saxpy(L,-cblas_sdot(L,&Y[2*t],2,&o[0],0)/L,&o[0],0,&Y[2*t],2);
                    cblas_saxpy(L,-cblas_sdot(L,&Y[2*t+1],2,&o[0],0)/L,&o[0],0,&Y[2*t+1],2);
                }
            }
        }
    }
    else if (dim==1u)
    {
        if (iscolmajor)
        {
            if (stp==floorf(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_ccopy(Tpre,&z[0],0,&Y[2*l*T],1); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_ccopy(Tpost,&z[0],0,&Y[2*(l*T+Tpre+Tmid)],1); }
                    if (Tmid>0) { cblas_ccopy(Tmid,&X[2*(ssl+Tpre*istp)],istp,&Y[2*(l*T+Tpre)],1); cblas_cscal(Tmid,&W[2*l],&Y[2*(l*T+Tpre)],1); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_ccopy(-ss,&z[0],0,&Y[2*t],T);
                    cblas_chbmv(CblasColMajor,CblasUpper,L+ss,0,&o[0],&W[-2*ss],1,&X[0],1,&z[0],&Y[-2*T*ss],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_chbmv(CblasColMajor,CblasUpper,L,0,&o[0],&W[0],1,&X[2*ss],1,&z[0],&Y[2*t],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_chbmv(CblasColMajor,CblasUpper,N-ss,0,&o[0],&W[0],1,&X[2*ss],1,&z[0],&Y[2*t],T);
                    cblas_ccopy(L-N+ss,&z[0],0,&Y[2*(t+T*(N-ss))],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0)
            {
                for (t=0; t<T; t++)
                {
                    cblas_saxpy(L,-cblas_sdot(L,&Y[2*t],2,&o[0],0)/L,&o[0],0,&Y[2*t],2);
                    cblas_saxpy(L,-cblas_sdot(L,&Y[2*t+1],2,&o[0],0)/L,&o[0],0,&Y[2*t+1],2);
                }
            }
        }
        else
        {
            if (stp==floorf(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_ccopy(Tpre,&z[0],0,&Y[2*l],L); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_ccopy(Tpost,&z[0],0,&Y[2*(l+L*(Tpre+Tmid))],L); }
                    if (Tmid>0) { cblas_ccopy(Tmid,&X[2*(ssl+Tpre*istp)],istp,&Y[2*(l+L*Tpre)],L); cblas_cscal(Tmid,&W[2*l],&Y[2*(l+L*Tpre)],L); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_ccopy(-ss,&z[0],0,&Y[2*t*L],1);
                    cblas_chbmv(CblasRowMajor,CblasUpper,L+ss,0,&o[0],&W[-2*ss],1,&X[0],1,&z[0],&Y[2*(t*L-ss)],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_chbmv(CblasRowMajor,CblasUpper,L,0,&o[0],&W[0],1,&X[2*ss],1,&z[0],&Y[2*t*L],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_chbmv(CblasRowMajor,CblasUpper,N-ss,0,&o[0],&W[0],1,&X[2*ss],1,&z[0],&Y[2*t*L],1);
                    cblas_ccopy(L-N+ss,&z[0],0,&Y[2*(t*L+N-ss)],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0)
            {
                for (t=0; t<T; t++)
                {
                    cblas_saxpy(L,-cblas_sdot(L,&Y[2*t*L],2,&o[0],0)/L,&o[0],0,&Y[2*t*L],2);
                    cblas_saxpy(L,-cblas_sdot(L,&Y[2*t*L+1],2,&o[0],0)/L,&o[0],0,&Y[2*t*L+1],2);
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in window_univar_c: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int window_univar_z (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const double *W, const int dim, const int c0, const double stp, const char mn0)
{
    const double z[2] = {0.0,0.0}, o[2] = {1.0,0.0};
    const int T = (dim==0u) ? C : R;
    const int L = (dim==0u) ? R : C;
    const int Lpre = L/2;      //nsamps before center samp
    const int istp = (int)stp;
    int ss = c0 - Lpre;        //start samp of current frame
    int l, t = 0;
    int Tpre, Tmid, Tpost;
    int ssl = ss, esl = ss + (T-1)*istp;

    //Checks
    if (N<1) { fprintf(stderr,"error in window_univar_z: N (length X) must be positive\n"); return 1; }
    if (R<1) { fprintf(stderr,"error in window_univar_z: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in window_univar_z: C (ncols Y) must be positive\n"); return 1; }
    if (c0<0) { fprintf(stderr,"error in window_univar_z: c0 (center samp of 1st frame) must be nonnegative\n"); return 1; }
    if (c0>N-1) { fprintf(stderr,"error in window_univar_z: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in window_univar_z: L (winlength) must be positive\n"); return 1; }
    if (L>=N) { fprintf(stderr,"error in window_univar_z: L (winlength) must be < N (length X)\n"); return 1; }
    if (stp<=0.0) { fprintf(stderr,"error in window_univar_z: stp (step size) must be positive\n"); return 1; }

    if (dim==0u)
    {
        if (iscolmajor)
        {
            if (stp==floor(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_zcopy(Tpre,&z[0],0,&Y[2*l],L); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_zcopy(Tpost,&z[0],0,&Y[2*(l+L*(Tpre+Tmid))],L); }
                    if (Tmid>0) { cblas_zcopy(Tmid,&X[2*(ssl+Tpre*istp)],istp,&Y[2*(l+L*Tpre)],L); cblas_zscal(Tmid,&W[2*l],&Y[2*(l+L*Tpre)],L); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_zcopy(-ss,&z[0],0,&Y[2*t*L],1);
                    cblas_zhbmv(CblasColMajor,CblasUpper,L+ss,0,&o[0],&W[-2*ss],1,&X[0],1,&z[0],&Y[2*(t*L-ss)],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_zhbmv(CblasColMajor,CblasUpper,L,0,&o[0],&W[0],1,&X[2*ss],1,&z[0],&Y[2*t*L],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_zhbmv(CblasColMajor,CblasUpper,N-ss,0,&o[0],&W[0],1,&X[2*ss],1,&z[0],&Y[2*t*L],1);
                    cblas_zcopy(L-N+ss,&z[0],0,&Y[2*(t*L+N-ss)],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0)
            {
                for (t=0; t<T; t++)
                {
                    cblas_daxpy(L,-cblas_ddot(L,&Y[2*t*L],2,&o[0],0)/L,&o[0],0,&Y[2*t*L],2);
                    cblas_daxpy(L,-cblas_ddot(L,&Y[2*t*L+1],2,&o[0],0)/L,&o[0],0,&Y[2*t*L+1],2);
                }
            }
        }
        else
        {
            if (stp==floor(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_zcopy(Tpre,&z[0],0,&Y[2*l*T],1); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_zcopy(Tpost,&z[0],0,&Y[2*(l*T+Tpre+Tmid)],1); }
                    if (Tmid>0) { cblas_zcopy(Tmid,&X[2*(ssl+Tpre*istp)],istp,&Y[2*(l*T+Tpre)],1); cblas_zscal(Tmid,&W[2*l],&Y[2*(l*T+Tpre)],1); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_zcopy(-ss,&z[0],0,&Y[2*t],T);
                    cblas_zhbmv(CblasRowMajor,CblasUpper,L+ss,0,&o[0],&W[-2*ss],1,&X[0],1,&z[0],&Y[-2*T*ss],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_zhbmv(CblasRowMajor,CblasUpper,L,0,&o[0],&W[0],1,&X[2*ss],1,&z[0],&Y[2*t],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_zhbmv(CblasRowMajor,CblasUpper,N-ss,0,&o[0],&W[0],1,&X[2*ss],1,&z[0],&Y[2*t],T);
                    cblas_zcopy(L-N+ss,&z[0],0,&Y[2*(t+T*(N-ss))],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0)
            {
                for (t=0; t<T; t++)
                {
                    cblas_daxpy(L,-cblas_ddot(L,&Y[2*t],2,&o[0],0)/L,&o[0],0,&Y[2*t],2);
                    cblas_daxpy(L,-cblas_ddot(L,&Y[2*t+1],2,&o[0],0)/L,&o[0],0,&Y[2*t+1],2);
                }
            }
        }
    }
    else if (dim==1u)
    {
        if (iscolmajor)
        {
            if (stp==floor(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_zcopy(Tpre,&z[0],0,&Y[2*l*T],1); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_zcopy(Tpost,&z[0],0,&Y[2*(l*T+Tpre+Tmid)],1); }
                    if (Tmid>0) { cblas_zcopy(Tmid,&X[2*(ssl+Tpre*istp)],istp,&Y[2*(l*T+Tpre)],1); cblas_zscal(Tmid,&W[2*l],&Y[2*(l*T+Tpre)],1); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_zcopy(-ss,&z[0],0,&Y[2*t],T);
                    cblas_zhbmv(CblasColMajor,CblasUpper,L+ss,0,&o[0],&W[-2*ss],1,&X[0],1,&z[0],&Y[-2*T*ss],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_zhbmv(CblasColMajor,CblasUpper,L,0,&o[0],&W[0],1,&X[2*ss],1,&z[0],&Y[2*t],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_zhbmv(CblasColMajor,CblasUpper,N-ss,0,&o[0],&W[0],1,&X[2*ss],1,&z[0],&Y[2*t],T);
                    cblas_zcopy(L-N+ss,&z[0],0,&Y[2*(t+T*(N-ss))],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0)
            {
                for (t=0; t<T; t++)
                {
                    cblas_daxpy(L,-cblas_ddot(L,&Y[2*t],2,&o[0],0)/L,&o[0],0,&Y[2*t],2);
                    cblas_daxpy(L,-cblas_ddot(L,&Y[2*t+1],2,&o[0],0)/L,&o[0],0,&Y[2*t+1],2);
                }
            }
        }
        else
        {
            if (stp==floor(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_zcopy(Tpre,&z[0],0,&Y[2*l],L); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_zcopy(Tpost,&z[0],0,&Y[2*(l+L*(Tpre+Tmid))],L); }
                    if (Tmid>0) { cblas_zcopy(Tmid,&X[2*(ssl+Tpre*istp)],istp,&Y[2*(l+L*Tpre)],L); cblas_zscal(Tmid,&W[2*l],&Y[2*(l+L*Tpre)],L); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_zcopy(-ss,&z[0],0,&Y[2*t*L],1);
                    cblas_zhbmv(CblasRowMajor,CblasUpper,L+ss,0,&o[0],&W[-2*ss],1,&X[0],1,&z[0],&Y[2*(t*L-ss)],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_zhbmv(CblasRowMajor,CblasUpper,L,0,&o[0],&W[0],1,&X[2*ss],1,&z[0],&Y[2*t*L],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_zhbmv(CblasRowMajor,CblasUpper,N-ss,0,&o[0],&W[0],1,&X[2*ss],1,&z[0],&Y[2*t*L],1);
                    cblas_zcopy(L-N+ss,&z[0],0,&Y[2*(t*L+N-ss)],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0)
            {
                for (t=0; t<T; t++)
                {
                    cblas_daxpy(L,-cblas_ddot(L,&Y[2*t*L],2,&o[0],0)/L,&o[0],0,&Y[2*t*L],2);
                    cblas_daxpy(L,-cblas_ddot(L,&Y[2*t*L+1],2,&o[0],0)/L,&o[0],0,&Y[2*t*L+1],2);
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in window_univar_z: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif

