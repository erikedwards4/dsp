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
//e.g., Y <- X[3] X[2] X[1] X[0] X[1] X[2] X[3] X[4] X[5] ... X[N-1]
//Use frame_univar_float to use zeros outside of [0 N-1].

//The following framing convention is forced here:
//Samples from one frame are always contiguous in memory, regardless of row- vs. col-major.
//So, if Y is row-major, then it has size W x L;
//but if Y is col-major, then it has size L x W.
//This avoids confusion, maximizes speed, minimizes code length,
//and remains compatible with the only viable way to stream data online.

//For more flexibility, see also frame_univar_float.c.

#include <stdio.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int frame_univar_s (float *Y, const float *X, const size_t N, const size_t W, const size_t L, const size_t stp);
int frame_univar_d (double *Y, const double *X, const size_t N, const size_t W, const size_t L, const size_t stp);
int frame_univar_c (float *Y, const float *X, const size_t N, const size_t W, const size_t L, const size_t stp);
int frame_univar_z (double *Y, const double *X, const size_t N, const size_t W, const size_t L, const size_t stp);


int frame_univar_s (float *Y, const float *X, const size_t N, const size_t W, const size_t L, const size_t stp)
{
    if (L<1u) { fprintf(stderr,"error in frame_univar_s: L must be positive\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in frame_univar_s: stp must be positive\n"); return 1; }

    //Establish snip_edges from W, N, L, stp
    //This also serves as an input check
    int snip_edges;
    if (W==(N+stp/2u)/stp) { snip_edges = 0; }
    else if (L>N) { fprintf(stderr,"error in frame_univar_s: L must be < N if snip_edges\n"); return 1; }
    else if (W==1u+(N-L)/stp) { snip_edges = 1; }
    else { fprintf(stderr,"error in frame_univar_s: W (num_frames) must be 1u+(N-L)/stp or (N+stp/2u)/stp\n"); return 1; }

    if (W==0u) {}
    else
    {
        const size_t xd = (L<stp) ? stp-L : L-stp;      //amount to decrement X by after each frame

        if (snip_edges)
        {
            for (size_t w=0u; w<W; ++w, X-=xd)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = *X; }
            }
        }
        else
        {
            const size_t Lpre = L/2u;                   //nsamps before center samp
            const size_t Lpost = L-L/2u-1u;             //nsamps after center samp
            size_t cs = 0u, w = 0u;                     //current center-samp and frame
            
            while (cs<Lpre && w<W)
            {
                X += Lpre - cs;
                for (size_t l=0u; l<Lpre-cs; ++l, --X, ++Y) { *Y = *X; }
                //for (size_t l=0u; l<Lpre-cs; ++l, ++Y) { *Y = 0.0f; }
                for (size_t l=Lpre-cs; l<L; ++l, ++X, ++Y) { *Y = *X; }
                X -= L-Lpre+cs; X += stp; cs += stp; ++w;
            }
            while(cs+Lpost<N && w<W)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = *X; }
                X -= xd; cs += stp; ++w;
            }
            while(cs<N+Lpre && w<W)
            {
                for (size_t l=0u; l<L-(N+Lpre-cs); ++l, ++X, ++Y) { *Y = *X; }
                //for (size_t l=N+Lpre-cs; l<L; ++l, ++Y) { *Y = 0.0f; }
                for (size_t l=L-(N+Lpre-cs); l<L; ++l, ++Y) { *Y = *--X; }
                X -= L-(N+Lpre-cs); X += N+Lpre-cs; X += stp; cs += stp; ++w;
            }
            while (w<W) //this should not occur in general
            {
                for (size_t l=0u; l<L; ++l, ++Y) { *Y = 0.0f; }
                ++w;
            }
        }
    }

    return 0;
}


int frame_univar_d (double *Y, const double *X, const size_t N, const size_t W, const size_t L, const size_t stp)
{
    if (L<1u) { fprintf(stderr,"error in frame_univar_d: L must be positive\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in frame_univar_d: stp must be positive\n"); return 1; }

    //Establish snip_edges from W, N, L, stp
    //This also serves as an input check
    int snip_edges;
    if (W==(N+stp/2u)/stp) { snip_edges = 0; }
    else if (L>N) { fprintf(stderr,"error in frame_univar_d: L must be < N if snip_edges\n"); return 1; }
    else if (W==1u+(N-L)/stp) { snip_edges = 1; }
    else { fprintf(stderr,"error in frame_univar_d: W (num_frames) must be 1u+(N-L)/stp or (N+stp/2u)/stp\n"); return 1; }

    if (W==0u) {}
    else
    {
        const size_t xd = (L<stp) ? stp-L : L-stp;      //amount to decrement X by after each frame

        if (snip_edges)
        {
            for (size_t w=0u; w<W; ++w, X-=xd)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = *X; }
            }
        }
        else
        {
            const size_t Lpre = L/2u;                   //nsamps before center samp
            const size_t Lpost = L-L/2u-1u;             //nsamps after center samp
            size_t cs = 0u, w = 0u;                     //current center-samp and frame
            
            while (cs<Lpre && w<W)
            {
                for (size_t l=0u; l<Lpre-cs; ++l, --X, ++Y) { *Y = *X; }
                //for (size_t l=0u; l<Lpre-cs; ++l, ++Y) { *Y = 0.0; }
                for (size_t l=Lpre-cs; l<L; ++l, ++X, ++Y) { *Y = *X; }
                X -= L-Lpre+cs; X += stp; cs += stp; ++w;
            }
            while(cs+Lpost<N && w<W)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = *X; }
                X -= xd; cs += stp; ++w;
            }
            while(cs<N+Lpre && w<W)
            {
                for (size_t l=0u; l<N+Lpre-cs; ++l, ++X, ++Y) { *Y = *X; }
                for (size_t l=N+Lpre-cs; l<L; ++l, ++Y) { *Y = 0.0; }
                X -= N+Lpre-cs; X += stp; cs += stp; ++w;
            }
            while (w<W)
            {
                for (size_t l=0u; l<L; ++l, ++Y) { *Y = 0.0; }
                ++w;
            }
        }
    }

    return 0;
}


int frame_univar_c (float *Y, const float *X, const size_t N, const size_t W, const size_t L, const size_t stp)
{
    if (L<1u) { fprintf(stderr,"error in frame_univar_c: L must be positive\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in frame_univar_c: stp must be positive\n"); return 1; }

    //Establish snip_edges from W, N, L, stp
    //This also serves as an input check
    int snip_edges;
    if (W==(N+stp/2u)/stp) { snip_edges = 0; }
    else if (L>N) { fprintf(stderr,"error in frame_univar_c: L must be < N if snip_edges\n"); return 1; }
    else if (W==1u+(N-L)/stp) { snip_edges = 1; }
    else { fprintf(stderr,"error in frame_univar_c: W (num_frames) must be 1u+(N-L)/stp or (N+stp/2u)/stp\n"); return 1; }

    if (W==0u) {}
    else
    {
        const size_t xd = (L<stp) ? 2u*(stp-L) : 2u*(L-stp);    //amount to decrement X by after each frame

        if (snip_edges)
        {
            for (size_t w=0u; w<W; ++w, X-=xd)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            }
        }
        else
        {
            const size_t Lpre = L/2u;                   //nsamps before center samp
            const size_t Lpost = L-L/2u-1u;             //nsamps after center samp
            size_t cs = 0u, w = 0u;                     //current center-samp and frame
            
            while (cs<Lpre && w<W)
            {
                for (size_t l=0u; l<Lpre-cs; ++l, ++Y) { *Y = 0.0f; *++Y = 0.0f; }
                for (size_t l=Lpre-cs; l<L; ++l, ++X, ++Y) { *Y = *X; *++Y = *++X; }
                X -= 2u*(L-Lpre+cs); X += 2u*stp; cs += stp; ++w;
            }
            while(cs+Lpost<N && w<W)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = *X; *++Y = *++X; }
                X -= xd; cs += stp; ++w;
            }
            while(cs<N+Lpre && w<W)
            {
                for (size_t l=0u; l<N+Lpre-cs; ++l, ++X, ++Y) { *Y = *X; *++Y = *++X; }
                for (size_t l=N+Lpre-cs; l<L; ++l, ++Y) { *Y = 0.0f; *++Y = 0.0f; }
                X -= 2u*(N+Lpre-cs); X += 2u*stp; cs += stp; ++w;
            }
            while (w<W)
            {
                for (size_t l=0u; l<L; ++l, ++Y) { *Y = 0.0f; *++Y = 0.0f; }
                ++w;
            }
        }
    }

    return 0;
}


int frame_univar_z (double *Y, const double *X, const size_t N, const size_t W, const size_t L, const size_t stp)
{
    if (L<1u) { fprintf(stderr,"error in frame_univar_z: L must be positive\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in frame_univar_z: stp must be positive\n"); return 1; }

    //Establish snip_edges from W, N, L, stp
    //This also serves as an input check
    int snip_edges;
    if (W==(N+stp/2u)/stp) { snip_edges = 0; }
    else if (L>N) { fprintf(stderr,"error in frame_univar_z: L must be < N if snip_edges\n"); return 1; }
    else if (W==1u+(N-L)/stp) { snip_edges = 1; }
    else { fprintf(stderr,"error in frame_univar_z: W (num_frames) must be 1u+(N-L)/stp or (N+stp/2u)/stp\n"); return 1; }

    if (W==0u) {}
    else
    {
        const size_t xd = (L<stp) ? 2u*(stp-L) : 2u*(L-stp);    //amount to decrement X by after each frame

        if (snip_edges)
        {
            for (size_t w=0u; w<W; ++w, X-=xd)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            }
        }
        else
        {
            const size_t Lpre = L/2u;                   //nsamps before center samp
            const size_t Lpost = L-L/2u-1u;             //nsamps after center samp
            size_t cs = 0u, w = 0u;                     //current center-samp and frame
            
            while (cs<Lpre && w<W)
            {
                for (size_t l=0u; l<Lpre-cs; ++l, ++Y) { *Y = 0.0; *++Y = 0.0; }
                for (size_t l=Lpre-cs; l<L; ++l, ++X, ++Y) { *Y = *X; *++Y = *++X; }
                X -= 2u*(L-Lpre+cs); X += 2u*stp; cs += stp; ++w;
            }
            while(cs+Lpost<N && w<W)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = *X; *++Y = *++X; }
                X -= xd; cs += stp; ++w;
            }
            while(cs<N+Lpre && w<W)
            {
                for (size_t l=0u; l<N+Lpre-cs; ++l, ++X, ++Y) { *Y = *X; *++Y = *++X; }
                for (size_t l=N+Lpre-cs; l<L; ++l, ++Y) { *Y = 0.0; *++Y = 0.0; }
                X -= 2u*(N+Lpre-cs); X += 2u*stp; cs += stp; ++w;
            }
            while (w<W)
            {
                for (size_t l=0u; l<L; ++l, ++Y) { *Y = 0.0; *++Y = 0.0; }
                ++w;
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
