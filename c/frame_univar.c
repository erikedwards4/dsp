//This takes a continuous uniivariate time series and outputs a set of frames.
//The frame rate (fr) and sample rate (fs) are allowed to be non-integer,
//but the center sample of each frame is always rounded to the nearest integer.

//Input c0 is the center-sample of the first frame (usually c0==0).

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int frame_univar_s (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const int dim, const int c0, const float stp);
int frame_univar_d (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const int dim, const int c0, const double stp);
int frame_univar_c (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const int dim, const int c0, const float stp);
int frame_univar_z (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const int dim, const int c0, const double stp);


int frame_univar_s (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const int dim, const int c0, const float stp)
{
    const float z = 0.0f;
    const int T = (dim==0) ? C : R;
    const int L = (dim==0) ? R : C;
    const int Lpre = L/2;      //nsamps before center samp
    const int istp = (int)stp;
    int ss = c0 - Lpre;        //start samp of current frame
    int l, t = 0;
    int Tpre, Tmid, Tpost;
    int ssl = ss, esl = ss + (T-1)*istp;

    //Checks
    if (N<1) { fprintf(stderr,"error in frame_univar_s: N (length X) must be positive\n"); return 1; }
    if (R<1) { fprintf(stderr,"error in frame_univar_s: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in frame_univar_s: C (ncols Y) must be positive\n"); return 1; }
    if (c0<0) { fprintf(stderr,"error in frame_univar_s: c0 (center samp of 1st frame) must be nonnegative\n"); return 1; }
    if (c0>N-1) { fprintf(stderr,"error in frame_univar_s: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in frame_univar_s: L must be positive\n"); return 1; }
    if (L>=N) { fprintf(stderr,"error in frame_univar_s: L must be < N (length X)\n"); return 1; }
    if (stp<=0.0f) { fprintf(stderr,"error in frame_univar_s: stp (step size) must be positive\n"); return 1; }

    if (dim==0)
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
                    if (Tmid>0) { cblas_scopy(Tmid,&X[ssl+Tpre*istp],istp,&Y[l+L*Tpre],L); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_scopy(-ss,&z,0,&Y[t*L],1);
                    cblas_scopy(L+ss,&X[0],1,&Y[t*L-ss],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_scopy(L,&X[ss],1,&Y[t*L],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_scopy(N-ss,&X[ss],1,&Y[t*L],1);
                    cblas_scopy(L-N+ss,&z,0,&Y[t*L+N-ss],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
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
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_scopy(Tpre,&z,0,&Y[l*T],1); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_scopy(Tpost,&z,0,&Y[l*T+Tpre+Tmid],1); }
                    if (Tmid>0) { cblas_scopy(Tmid,&X[ssl+Tpre*istp],istp,&Y[l*T+Tpre],1); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_scopy(-ss,&z,0,&Y[t],T);
                    cblas_scopy(L+ss,&X[0],1,&Y[-T*ss],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_scopy(L,&X[ss],1,&Y[t],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_scopy(N-ss,&X[ss],1,&Y[t],T);
                    cblas_scopy(L-N+ss,&z,0,&Y[t+T*(N-ss)],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            if (stp==floorf(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_scopy(Tpre,&z,0,&Y[l*T],1); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_scopy(Tpost,&z,0,&Y[l*T+Tpre+Tmid],1); }
                    if (Tmid>0) { cblas_scopy(Tmid,&X[ssl+Tpre*istp],istp,&Y[l*T+Tpre],1); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_scopy(-ss,&z,0,&Y[t],T);
                    cblas_scopy(L+ss,&X[0],1,&Y[-T*ss],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_scopy(L,&X[ss],1,&Y[t],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_scopy(N-ss,&X[ss],1,&Y[t],T);
                    cblas_scopy(L-N+ss,&z,0,&Y[t+T*(N-ss)],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
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
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_scopy(Tpre,&z,0,&Y[l],L); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_scopy(Tpost,&z,0,&Y[l+L*(Tpre+Tmid)],L); }
                    if (Tmid>0) { cblas_scopy(Tmid,&X[ssl+Tpre*istp],istp,&Y[l+L*Tpre],L); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_scopy(-ss,&z,0,&Y[t*L],1);
                    cblas_scopy(L+ss,&X[0],1,&Y[t*L-ss],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_scopy(L,&X[ss],1,&Y[t*L],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_scopy(N-ss,&X[ss],1,&Y[t*L],1);
                    cblas_scopy(L-N+ss,&z,0,&Y[t*L+N-ss],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in frame_univar_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int frame_univar_d (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const int dim, const int c0, const double stp)
{
    const double z = 0.0;
    const int T = (dim==0) ? C : R;
    const int L = (dim==0) ? R : C;
    const int Lpre = L/2;      //nsamps before center samp
    const int istp = (int)stp;
    int ss = c0 - Lpre;        //start samp of current frame
    int l, t = 0;
    int Tpre, Tmid, Tpost;
    int ssl = ss, esl = ss + (T-1)*istp;

    //Checks
    if (N<1) { fprintf(stderr,"error in frame_univar_d: N (length X) must be positive\n"); return 1; }
    if (R<1) { fprintf(stderr,"error in frame_univar_d: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in frame_univar_d: C (ncols Y) must be positive\n"); return 1; }
    if (c0<0) { fprintf(stderr,"error in frame_univar_d: c0 (center samp of 1st frame) must be nonnegative\n"); return 1; }
    if (c0>N-1) { fprintf(stderr,"error in frame_univar_d: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in frame_univar_d: L must be positive\n"); return 1; }
    if (L>=N) { fprintf(stderr,"error in frame_univar_d: L must be < N (length X)\n"); return 1; }
    if (stp<=0.0) { fprintf(stderr,"error in frame_univar_d: stp (step size) must be positive\n"); return 1; }

    if (dim==0)
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
                    if (Tmid>0) { cblas_dcopy(Tmid,&X[ssl+Tpre*istp],istp,&Y[l+L*Tpre],L); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_dcopy(-ss,&z,0,&Y[t*L],1);
                    cblas_dcopy(L+ss,&X[0],1,&Y[t*L-ss],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_dcopy(L,&X[ss],1,&Y[t*L],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_dcopy(N-ss,&X[ss],1,&Y[t*L],1);
                    cblas_dcopy(L-N+ss,&z,0,&Y[t*L+N-ss],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
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
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_dcopy(Tpre,&z,0,&Y[l*T],1); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_dcopy(Tpost,&z,0,&Y[l*T+Tpre+Tmid],1); }
                    if (Tmid>0) { cblas_dcopy(Tmid,&X[ssl+Tpre*istp],istp,&Y[l*T+Tpre],1); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_dcopy(-ss,&z,0,&Y[t],T);
                    cblas_dcopy(L+ss,&X[0],1,&Y[-T*ss],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_dcopy(L,&X[ss],1,&Y[t],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_dcopy(N-ss,&X[ss],1,&Y[t],T);
                    cblas_dcopy(L-N+ss,&z,0,&Y[t+T*(N-ss)],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
            }
        }
    }
    else if (dim==1)
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
                    if (Tmid>0) { cblas_dcopy(Tmid,&X[ssl+Tpre*istp],istp,&Y[l*T+Tpre],1); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_dcopy(-ss,&z,0,&Y[t],T);
                    cblas_dcopy(L+ss,&X[0],1,&Y[-T*ss],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_dcopy(L,&X[ss],1,&Y[t],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_dcopy(N-ss,&X[ss],1,&Y[t],T);
                    cblas_dcopy(L-N+ss,&z,0,&Y[t+T*(N-ss)],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
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
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_dcopy(Tpre,&z,0,&Y[l],L); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_dcopy(Tpost,&z,0,&Y[l+L*(Tpre+Tmid)],L); }
                    if (Tmid>0) { cblas_dcopy(Tmid,&X[ssl+Tpre*istp],istp,&Y[l+L*Tpre],L); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_dcopy(-ss,&z,0,&Y[t*L],1);
                    cblas_dcopy(L+ss,&X[0],1,&Y[t*L-ss],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_dcopy(L,&X[ss],1,&Y[t*L],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_dcopy(N-ss,&X[ss],1,&Y[t*L],1);
                    cblas_dcopy(L-N+ss,&z,0,&Y[t*L+N-ss],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in frame_univar_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int frame_univar_c (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const int dim, const int c0, const float stp)
{
    const float z[2] = {0.0f,0.0f};
    const int T = (dim==0) ? C : R;
    const int L = (dim==0) ? R : C;
    const int Lpre = L/2;      //nsamps before center samp
    const int istp = (int)stp;
    int ss = c0 - Lpre;        //start samp of current frame
    int l, t = 0;
    int Tpre, Tmid, Tpost;
    int ssl = ss, esl = ss + (T-1)*istp;

    //Checks
    if (N<1) { fprintf(stderr,"error in frame_univar_c: N (length X) must be positive\n"); return 1; }
    if (R<1) { fprintf(stderr,"error in frame_univar_c: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in frame_univar_c: C (ncols Y) must be positive\n"); return 1; }
    if (c0<0) { fprintf(stderr,"error in frame_univar_c: c0 (center samp of 1st frame) must be nonnegative\n"); return 1; }
    if (c0>N-1) { fprintf(stderr,"error in frame_univar_c: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in frame_univar_c: L must be positive\n"); return 1; }
    if (L>=N) { fprintf(stderr,"error in frame_univar_c: L must be < N (length X)\n"); return 1; }
    if (stp<=0.0f) { fprintf(stderr,"error in frame_univar_c: stp (step size) must be positive\n"); return 1; }

    if (dim==0)
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
                    if (Tmid>0) { cblas_ccopy(Tmid,&X[2*(ssl+Tpre*istp)],istp,&Y[2*(l+L*Tpre)],L); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_ccopy(-ss,&z[0],0,&Y[2*t*L],1);
                    cblas_ccopy(L+ss,&X[0],1,&Y[2*(t*L-ss)],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_ccopy(L,&X[2*ss],1,&Y[2*t*L],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_ccopy(N-ss,&X[2*ss],1,&Y[2*t*L],1);
                    cblas_ccopy(L-N+ss,&z[0],0,&Y[2*(t*L+N-ss)],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
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
                    if (Tmid>0) { cblas_ccopy(Tmid,&X[2*(ssl+Tpre*istp)],istp,&Y[2*(l*T+Tpre)],1); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_ccopy(-ss,&z[0],0,&Y[2*t],T);
                    cblas_ccopy(L+ss,&X[0],1,&Y[-2*T*ss],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_ccopy(L,&X[2*ss],1,&Y[2*t],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_ccopy(N-ss,&X[2*ss],1,&Y[2*t],T);
                    cblas_ccopy(L-N+ss,&z[0],0,&Y[2*(t+T*(N-ss))],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
            }
        }
    }
    else if (dim==1)
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
                    if (Tmid>0) { cblas_ccopy(Tmid,&X[2*(ssl+Tpre*istp)],istp,&Y[2*(l*T+Tpre)],1); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_ccopy(-ss,&z[0],0,&Y[2*t],T);
                    cblas_ccopy(L+ss,&X[0],1,&Y[-2*T*ss],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_ccopy(L,&X[2*ss],1,&Y[2*t],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_ccopy(N-ss,&X[2*ss],1,&Y[2*t],T);
                    cblas_ccopy(L-N+ss,&z[0],0,&Y[2*(t+T*(N-ss))],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
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
                    if (Tmid>0) { cblas_ccopy(Tmid,&X[2*(ssl+Tpre*istp)],istp,&Y[2*(l+L*Tpre)],L); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_ccopy(-ss,&z[0],0,&Y[2*t*L],1);
                    cblas_ccopy(L+ss,&X[0],1,&Y[2*(t*L-ss)],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_ccopy(L,&X[2*ss],1,&Y[2*t*L],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_ccopy(N-ss,&X[2*ss],1,&Y[2*t*L],1);
                    cblas_ccopy(L-N+ss,&z[0],0,&Y[2*(t*L+N-ss)],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in frame_univar_c: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int frame_univar_z (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const int dim, const int c0, const double stp)
{
    const double z[2] = {0.0,0.0};
    const int T = (dim==0) ? C : R;
    const int L = (dim==0) ? R : C;
    const int Lpre = L/2;      //nsamps before center samp
    const int istp = (int)stp;
    int ss = c0 - Lpre;        //start samp of current frame
    int l, t = 0;
    int Tpre, Tmid, Tpost;
    int ssl = ss, esl = ss + (T-1)*istp;

    //Checks
    if (N<1) { fprintf(stderr,"error in frame_univar_z: N (length X) must be positive\n"); return 1; }
    if (R<1) { fprintf(stderr,"error in frame_univar_z: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in frame_univar_z: C (ncols Y) must be positive\n"); return 1; }
    if (c0<0) { fprintf(stderr,"error in frame_univar_z: c0 (center samp of 1st frame) must be nonnegative\n"); return 1; }
    if (c0>N-1) { fprintf(stderr,"error in frame_univar_z: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in frame_univar_z: L must be positive\n"); return 1; }
    if (L>=N) { fprintf(stderr,"error in frame_univar_z: L must be < N (length X)\n"); return 1; }
    if (stp<=0.0) { fprintf(stderr,"error in frame_univar_z: stp (step size) must be positive\n"); return 1; }

    if (dim==0)
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
                    if (Tmid>0) { cblas_zcopy(Tmid,&X[2*(ssl+Tpre*istp)],istp,&Y[2*(l+L*Tpre)],L); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_zcopy(-ss,&z[0],0,&Y[2*t*L],1);
                    cblas_zcopy(L+ss,&X[0],1,&Y[2*(t*L-ss)],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_zcopy(L,&X[2*ss],1,&Y[2*t*L],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_zcopy(N-ss,&X[2*ss],1,&Y[2*t*L],1);
                    cblas_zcopy(L-N+ss,&z[0],0,&Y[2*(t*L+N-ss)],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
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
                    if (Tmid>0) { cblas_zcopy(Tmid,&X[2*(ssl+Tpre*istp)],istp,&Y[2*(l*T+Tpre)],1); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_zcopy(-ss,&z[0],0,&Y[2*t],T);
                    cblas_zcopy(L+ss,&X[0],1,&Y[-2*T*ss],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_zcopy(L,&X[2*ss],1,&Y[2*t],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_zcopy(N-ss,&X[2*ss],1,&Y[2*t],T);
                    cblas_zcopy(L-N+ss,&z[0],0,&Y[2*(t+T*(N-ss))],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
            }
        }
    }
    else if (dim==1)
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
                    if (Tmid>0) { cblas_zcopy(Tmid,&X[2*(ssl+Tpre*istp)],istp,&Y[2*(l*T+Tpre)],1); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_zcopy(-ss,&z[0],0,&Y[2*t],T);
                    cblas_zcopy(L+ss,&X[0],1,&Y[-2*T*ss],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_zcopy(L,&X[2*ss],1,&Y[2*t],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_zcopy(N-ss,&X[2*ss],1,&Y[2*t],T);
                    cblas_zcopy(L-N+ss,&z[0],0,&Y[2*(t+T*(N-ss))],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
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
                    if (Tmid>0) { cblas_zcopy(Tmid,&X[2*(ssl+Tpre*istp)],istp,&Y[2*(l+L*Tpre)],L); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_zcopy(-ss,&z[0],0,&Y[2*t*L],1);
                    cblas_zcopy(L+ss,&X[0],1,&Y[2*(t*L-ss)],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_zcopy(L,&X[2*ss],1,&Y[2*t*L],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_zcopy(N-ss,&X[2*ss],1,&Y[2*t*L],1);
                    cblas_zcopy(L-N+ss,&z[0],0,&Y[2*(t*L+N-ss)],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in frame_univar_z: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif

