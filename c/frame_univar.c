//This takes a continuous uniivariate time series and outputs a set of frames.
//The frame rate (fr) and sample rate (fs) are allowed to be non-integer,
//but the center sample of each frame is always rounded to the nearest integer.

//Input c0 is the center-sample of the first frame (usually c0==0).

#include <stdio.h>
#include <float.h>
#include <math.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int frame_univar_s (float *Y, const float *X, const size_t N, const size_t L, const float c0, const float stp);
int frame_univar_int_s (float *Y, const float *X, const size_t N, const size_t L, const size_t c0, const size_t stp);


int frame_univar_s (float *Y, const float *X, const size_t N, const size_t L, const float c0, const float stp)
{
    if (L>=N) { fprintf(stderr,"error in frame_univar_s: L must be < N (length X)\n"); return 1; }
    if (L<1u) { fprintf(stderr,"error in frame_univar_s: L must be positive\n"); return 1; }
    if (N<1u) { fprintf(stderr,"error in frame_univar_s: N (length X) must be positive\n"); return 1; }
    if (c0>N-1u) { fprintf(stderr,"error in frame_univar_s: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (stp<FLT_EPSILON) { fprintf(stderr,"error in frame_univar_s: stp (step size) must be positive\n"); return 1; }

    const size_t Lpre = L/2u;               //nsamps before center samp
    const size_t Lpost = L-L/2u-1u;         //nsamps after center samp
    float st, ct, et;                       //start, center, end times of current frame (floats, but still in units of samps)
    size_t ss, cs, es;                      //start, center, end samps of current frame (after rounding)
    const float Tpre = (float)Lpre, Tpost = (float)Lpost;

    //const size_t istp = (size_t)stp;
    //size_t ss = c0 - Lpre;         //start samp of current frame
    //size_t t = 0;
    //size_t Tpre, Tmid, Tpost;
    //size_t ssl = ss, esl = ss + (T-1u)*istp;

    //float ct = (float)c0;

    et = c0 + Tpost;
    while (et<0.5f)
    {
        for (size_t l=0u; l<L; ++l, ++Y) { *Y = 0.0f; }
        et += stp;
    }
    ct = et - Tpost;
    st = ct - Tpre;
    while (st<0.5f)
    {
        for (size_t l=0u; l<)
    }
    cs = roundf(et-Tpost);
    while (et<(float)N-0.5f)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = *X; }
        et += stp;
        cs = (size_t)(roundf(et-Tpost));
        es = cs + Lpost;
        
        prevcs = cs; cs = (size_t)(roundf(et-Tpost));
        prevss = ss; ss = cs - Lpre;
        preves = es; es = cs + Lpost;
        X -= (size_t)(roundf(et-Tpost)-roundf(prevct));
        cs = roundf(et-Tpost);
    }

    cs = c0;
    while (cs<N)
    {
        if (cs<Lpre)
        {
            for (size_t l=0u; l<Lpre-cs; ++l, ++Y) { *Y = 0.0f; }
            for (size_t l=Lpre-cs; l<L; ++l, ++X, ++Y) { *Y = *X; }
        }
        else if (N-cs<Lpost) {}
        else
        {
            for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = *X; }
        }
        ct += stp; cs = (size_t)stp;
    }

    if (c0<Lpre)
    {
        for (size_t l=0u; l<Lpre-c0; ++l, ++Y) { *Y = 0.0f; }
        for (size_t l=Lpre-c0; l<L; ++l, ++X, ++Y) { *Y = *X; }
    }

    n = c0;

    if (stp==floorf(stp))
    {
        for (size_t l=0u; l<L; ++l, ++ssl, ++esl)
        {
            Tmid = T;
            if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_scopy(Tpre,&z,0,&Y[l],L); } else { Tpre = 0; }
            if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_scopy(Tpost,&z,0,&Y[l+L*(Tpre+Tmid)],L); }
            if (Tmid>0) { cblas_scopy(Tmid,&X[ssl+Tpre*istp],istp,&Y[l+L*Tpre],L); }
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

    if (L==1u)
    {
        for (size_t n=c0; )
    }
}



#ifdef __cplusplus
}
}
#endif
