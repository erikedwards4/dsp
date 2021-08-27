//This takes a continuous univariate time series (vector X),
//and outputs a set of frames (matrix Y).

//There is an established convention from HTK through Kaldi to Librosa,
//involving a setting "snip-edges", which usually has default "true".
//For compatibility, ease of use, and to avoid introducing a new set of conventions,
//the number of frames (W) is controlled here by "snip-edges",
//using the identical formula as used in Kaldi (NumFrames in feature-window.cc).
//Here: W=num_frames, N=num_samples, L=frame_length, and stp=frame_shift

//Also for compatibility to Kaldi and Librosa,
//frames that overlap the edge of X are filled in by flipping the edge of X,
//e.g., Y <- X[3] X[2] X[1] X[0] X[1] X[2] X[3] X[4] X[5] ...
//Use frame_univar_flt to use zeros outside of [0 N-1].

//Also for compatibility to Kaldi, for snip_edges=false,
//the start-samp of the first frame is at stp/2 - L/2.

//The following framing convention is forced here:
//Samples from one frame are always contiguous in memory, regardless of row- vs. col-major.
//So, if Y is row-major, then it has size W x L;
//but if Y is col-major, then it has size L x W.
//This avoids confusion, maximizes speed, minimizes code length,
//and remains compatible with the only viable way to stream data online.

//For more flexibility, see also frame_univar_flt.c.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int frame_univar_s (float *Y, const float *X, const size_t N, const size_t L, const size_t stp, const int snip_edges);
int frame_univar_d (double *Y, const double *X, const size_t N, const size_t L, const size_t stp, const int snip_edges);
int frame_univar_c (float *Y, const float *X, const size_t N, const size_t L, const size_t stp, const int snip_edges);
int frame_univar_z (double *Y, const double *X, const size_t N, const size_t L, const size_t stp, const int snip_edges);


int frame_univar_s (float *Y, const float *X, const size_t N, const size_t L, const size_t stp, const int snip_edges)
{
    if (L<1u) { fprintf(stderr,"error in frame_univar_s: L must be positive\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in frame_univar_s: stp must be positive\n"); return 1; }
    if (snip_edges && L>N) { fprintf(stderr,"error in frame_univar_s: L must be < N if snip_edges\n"); return 1; }

    //Set number of frames (W)
    const size_t W = (snip_edges) ? 1u+(N-L)/stp : (N+stp/2u)/stp;

    if (W==0u) {}
    else if (snip_edges)
    {
        const int xd = (int)L - (int)stp;           //X inc after each frame
        
        for (size_t w=W; w>0u; --w, X-=xd)
        {
            for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
        }
    }
    else
    {
        const int xd = (int)L - (int)stp;           //X inc after each frame
        const size_t Lpre = L/2u;                   //nsamps before center samp
        int ss = (int)(stp/2u) - (int)Lpre;         //start-samp of current frame
        int n, prev_n = 0;                          //current/prev samps in X

        for (size_t w=W; w>0u; --w, ss+=stp)
        {
            if (ss<0 || ss>(int)N-(int)L)
            {
                for (int s=ss; s<ss+(int)L; ++s, ++Y)
                {
                    //This ensures extrapolation by signal reversal to any length
                    n = s;
                    while (n<0 || n>=(int)N)
                    {
                        n = (n<0) ? -n-1 : (n<(int)N) ? n : 2*(int)N-1-n;
                    }
                    X += n - prev_n;
                    *Y = *X;
                    prev_n = n;
                }
            }
            else
            {
                X += ss - prev_n;
                for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
                X -= xd;
                prev_n = ss + (int)stp;
            }
        }
    }

    return 0;
}


int frame_univar_d (double *Y, const double *X, const size_t N, const size_t L, const size_t stp, const int snip_edges)
{
    if (L<1u) { fprintf(stderr,"error in frame_univar_d: L must be positive\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in frame_univar_d: stp must be positive\n"); return 1; }
    if (snip_edges && L>N) { fprintf(stderr,"error in frame_univar_d: L must be < N if snip_edges\n"); return 1; }

    //Set number of frames (W)
    const size_t W = (snip_edges) ? 1u+(N-L)/stp : (N+stp/2u)/stp;

    if (W==0u) {}
    else if (snip_edges)
    {
        const int xd = (int)L - (int)stp;           //X inc after each frame

        for (size_t w=W; w>0u; --w, X-=xd)
        {
            for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
        }
    }
    else
    {
        const int xd = (int)L - (int)stp;           //X inc after each frame
        const size_t Lpre = L/2u;                   //nsamps before center samp
        int ss = (int)(stp/2u) - (int)Lpre;         //start-samp of current frame
        int n, prev_n = 0;                          //current/prev samps in X

        for (size_t w=W; w>0u; --w, ss+=stp)
        {
            if (ss<0 || ss>(int)N-(int)L)
            {
                for (int s=ss; s<ss+(int)L; ++s, ++Y)
                {
                    //This ensures extrapolation by signal reversal to any length
                    n = s;
                    while (n<0 || n>=(int)N)
                    {
                        n = (n<0) ? -n-1 : (n<(int)N) ? n : 2*(int)N-1-n;
                    }
                    X += n - prev_n;
                    *Y = *X;
                    prev_n = n;
                }
            }
            else
            {
                X += ss - prev_n;
                for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
                X -= xd;
                prev_n = ss + (int)stp;
            }
        }
    }

    return 0;
}


int frame_univar_c (float *Y, const float *X, const size_t N, const size_t L, const size_t stp, const int snip_edges)
{
    if (L<1u) { fprintf(stderr,"error in frame_univar_c: L must be positive\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in frame_univar_c: stp must be positive\n"); return 1; }
    if (snip_edges && L>N) { fprintf(stderr,"error in frame_univar_c: L must be < N if snip_edges\n"); return 1; }

    //Set number of frames (W)
    const size_t W = (snip_edges) ? 1u+(N-L)/stp : (N+stp/2u)/stp;

    if (W==0u) {}
    else if (snip_edges)
    {
        const int xd = 2*((int)L-(int)stp);         //X inc after each frame

        for (size_t w=W; w>0u; --w, X-=xd)
        {
            for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; *++Y = *++X; }
        }
    }
    else
    {
        const int xd = 2*((int)L-(int)stp);         //X inc after each frame
        const size_t Lpre = L/2u;                   //nsamps before center samp
        int ss = (int)(stp/2u) - (int)Lpre;         //start-samp of current frame
        int n, prev_n = 0;                          //current/prev samps in X

        for (size_t w=W; w>0u; --w, ss+=stp)
        {
            if (ss<0 || ss>(int)N-(int)L)
            {
                for (int s=ss; s<ss+(int)L; ++s, ++Y)
                {
                    //This ensures extrapolation by signal reversal to any length
                    n = s;
                    while (n<0 || n>=(int)N)
                    {
                        n = (n<0) ? -n-1 : (n<(int)N) ? n : 2*(int)N-1-n;
                    }
                    X += 2*(n-prev_n);
                    *Y = *X; *++Y = *(X+1);
                    prev_n = n;
                }
            }
            else
            {
                X += 2*(ss-prev_n);
                for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; *++Y = *++X; }
                X -= xd;
                prev_n = ss + (int)stp;
            }
        }
    }

    return 0;
}


int frame_univar_z (double *Y, const double *X, const size_t N, const size_t L, const size_t stp, const int snip_edges)
{
    if (L<1u) { fprintf(stderr,"error in frame_univar_z: L must be positive\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in frame_univar_z: stp must be positive\n"); return 1; }
    if (snip_edges && L>N) { fprintf(stderr,"error in frame_univar_z: L must be < N if snip_edges\n"); return 1; }

    //Set number of frames (W)
    const size_t W = (snip_edges) ? 1u+(N-L)/stp : (N+stp/2u)/stp;

    if (W==0u) {}
    else if (snip_edges)
    {
        const int xd = 2*((int)L-(int)stp);         //X inc after each frame

        for (size_t w=W; w>0u; --w, X-=xd)
        {
            for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; *++Y = *++X; }
        }
    }
    else
    {
        const int xd = 2*((int)L-(int)stp);         //X inc after each frame
        const size_t Lpre = L/2u;                   //nsamps before center samp
        int ss = (int)(stp/2u) - (int)Lpre;         //start-samp of current frame
        int n, prev_n = 0;                          //current/prev samps in X

        for (size_t w=W; w>0u; --w, ss+=stp)
        {
            if (ss<0 || ss>(int)N-(int)L)
            {
                for (int s=ss; s<ss+(int)L; ++s, ++Y)
                {
                    //This ensures extrapolation by signal reversal to any length
                    n = s;
                    while (n<0 || n>=(int)N)
                    {
                        n = (n<0) ? -n-1 : (n<(int)N) ? n : 2*(int)N-1-n;
                    }
                    X += 2*(n-prev_n);
                    *Y = *X; *++Y = *(X+1);
                    prev_n = n;
                }
            }
            else
            {
                X += 2*(ss-prev_n);
                for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; *++Y = *++X; }
                X -= xd;
                prev_n = ss + (int)stp;
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
