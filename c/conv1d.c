//1D convolution of each vector in X1 by X2.
//Each vector in X1 has length L1. X2 has length L2 (kernel_size).

//This is identical to xcorr1d, except X2 is flipped.
//Note that some "convolution" functions out there actually do cross-correlation.
//Note that this does NOT attempt to emmulate Octave's conv (but is similar in principle).

//This emulates PyTorch Conv1d, which actually does cross-correlation.
//X2 is the convolving kernel with length kernel_size.
//Each vector in X1 has length L_in, and each vector in Y has length L_out, set by:
//L_out = floor[1 + (L_in + 2*pad - dil*(L2-1) - 1)/stride].
//L_out = floor[1 + (N_in - N_out)/stride].
//where:
//N_in is the full length of input vectors including padding.
//N_out is the full length of X2 including dilation.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int conv1d_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int pad, const size_t str, const size_t dil, const size_t dim);
int conv1d_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int pad, const size_t str, const size_t dil, const size_t dim);
int conv1d_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int pad, const size_t str, const size_t dil, const size_t dim);
int conv1d_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int pad, const size_t str, const size_t dil, const size_t dim);


int conv1d_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int pad, const size_t str, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv1d_s: dim must be in [0 3]\n"); return 1; }
    if (str<1u) { fprintf(stderr,"error in conv1d_s: str (stride) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in conv1d_s: dil (dilation) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in conv1d_s: L2 (length of X2) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in conv1d_s: N (total length of X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in conv1d_s: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if ((int)L1+2*pad<=(int)(dil*(L2-1u))) { fprintf(stderr,"error in conv1d_s: L1+2*pad must be > dil*(L2-1)\n"); return 1; }

    const int N1 = (int)L1 + 2*pad;             //full length of vecs in X1 including padding
    const int N2 = (int)(dil*(L2-1u)) + 1;      //full length of vec in X2 including dilation
    const int inc = (int)str - (int)(dil*L2);   //fixed increment for X1 below
    size_t w=0u, W;                             //current frame and total frames
    int ss, es;                                 //current start-samp, end-samp
    float sm;                                   //intermediate sum

    //Set W (L_out)
    W = 1u + (size_t)(N1-N2)/str;

    if (W==0u) {}
    else if (L1==N)
    {
        //Init ss, es
        ss = -pad; es = ss + N2 - 1;

        //X2 before first samp of X1
        while (es<0 && w<W) { *Y++ = 0.0f; es+=str; ++w; }
        ss = es - N2 + 1;

        //X2 overlaps first samp of X1
        while (ss<0 && es<(int)L1 && w<W)
        {
            sm = 0.0f; X1 += es;
            for (int n=es; n>=0; n-=dil, X1-=dil, ++X2) { sm += *X1 * *X2; }
            *Y++ = sm;
            X1 += dil - (size_t)es%dil;
            X2 -= 1 + es/(int)dil;
            ss+=str; es+=str; ++w;
        }
        X1 += ss; X2 += L2 - 1u;

        if (N2>(int)L1) //X1 fully within X2
        {
            while (ss<0 && w<W)
            {
                sm = 0.0f;
                for (int n=(int)dil*(-ss/(int)dil); n<(int)L1-ss; n+=dil)
                {
                    sm += X1[n+ss] * X2[-n/(int)dil];
                }
                *Y++ = sm;
                ss+=str; ++w;
            }
            X1 += ss;
            es = ss + N2 - 1;
        }
        else            //X2 fully within X1
        {        
            while (es<(int)L1 && w<W)
            {
                sm = 0.0f;
                for (size_t l=L2; l>0u; --l, X1+=dil, --X2) { sm += *X1 * *X2; }
                *Y++ = sm;
                X1 += inc; X2 += L2;
                es+=str; ++w;
            }
            ss = es - N2 + 1;
        }

        //X2 overlaps end samp of X1
        while (ss<(int)L1 && w<W)
        {
            size_t l = 0u; sm = 0.0f;
            for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=dil, --X2) { sm += *X1 * *X2; }
            *Y++ = sm;
            X1 += (int)str - (int)(dil*l); X2 += l;
            ss+=str; es+=str; ++w;
        }

        //X2 past end samp of X1
        while (w<W) { *Y++ = 0.0f; ++w; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L1, G = V/B;

        for (size_t g=G; g>0u; --g, X1+=B*(L1-1u), Y+=B*(W-1u))
        {
            for (size_t b=B; b>0u; --b, ++X1, Y-=K*W-1u)
            {
                //Init ss, es, w
                ss = -pad; es = ss + N2 - 1; w = 0u;

                //X2 before first samp of X1
                while (es<0 && w<W) { *Y = 0.0f; Y+=K; es+=str; ++w; }
                ss = es - N2 + 1;

                //X2 overlaps first samp of X1
                while (ss<0 && es<(int)L1 && w<W)
                {
                    sm = 0.0f; X1 += es*(int)K;
                    for (int n=es; n>=0; n-=dil, X1-=dil*K, ++X2) { sm += *X1 * *X2; }
                    *Y = sm; Y += K;
                    X1 += K*(dil-(size_t)es%dil);
                    X2 -= 1 + es/(int)dil;
                    ss+=str; es+=str; ++w;
                }
                X1 += ss*(int)K; X2 += L2 - 1u;

                if (N2>(int)L1) //X1 fully within X2
                {
                    while (ss<0 && w<W)
                    {
                        sm = 0.0f;
                        for (int n=(int)dil*(-ss/(int)dil); n<(int)L1-ss; n+=dil)
                        {
                            sm += X1[(int)K*(n+ss)] * X2[-n/(int)dil];
                        }
                        *Y = sm; Y += K;
                        ss+=str; ++w;
                    }
                    X1 += ss*(int)K;
                    es = ss + N2 - 1;
                }
                else            //X2 fully within X1
                {        
                    while (es<(int)L1 && w<W)
                    {
                        sm = 0.0f;
                        for (size_t l=L2; l>0u; --l, X1+=dil*K, --X2) { sm += *X1 * *X2; }
                        *Y = sm; Y += K;
                        X1 += inc*(int)K; X2 += L2;
                        es+=str; ++w;
                    }
                    ss = es - N2 + 1;
                }

                //X2 overlaps end samp of X1
                while (ss<(int)L1 && w<W)
                {
                    size_t l = 0u; sm = 0.0f;
                    for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=dil*K, --X2) { sm += *X1 * *X2; }
                    *Y = sm; Y += K;
                    X1 += (int)K*((int)str-(int)(dil*l));
                    X2 += l;
                    ss+=str; es+=str; ++w;
                }

                //X2 past end samp of X1
                while (w<W) { *Y = 0.0f; Y += K; ++w; }

                //Reset X1, X2
                X1 -= (int)K*ss; X2 -= L2 - 1u;
            }
        }
    }

    return 0;
}


int conv1d_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int pad, const size_t str, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv1d_d: dim must be in [0 3]\n"); return 1; }
    if (str<1u) { fprintf(stderr,"error in conv1d_d: str (stride) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in conv1d_d: dil (dilation) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in conv1d_d: L2 (length of X2) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in conv1d_d: N (total length of X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in conv1d_d: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if ((int)L1+2*pad<=(int)(dil*(L2-1u))) { fprintf(stderr,"error in conv1d_d: L1+2*pad must be > dil*(L2-1)\n"); return 1; }

    const int N1 = (int)L1 + 2*pad;             //full length of vecs in X1 including padding
    const int N2 = (int)(dil*(L2-1u)) + 1;      //full length of vec in X2 including dilation
    const int inc = (int)str - (int)(dil*L2);   //fixed increment for X1 below
    size_t w=0u, W;                             //current frame and total frames
    int ss, es;                                 //current start-samp, end-samp
    double sm;                                   //intermediate sum

    //Set W (L_out)
    W = 1u + (size_t)(N1-N2)/str;

    if (W==0u) {}
    else if (L1==N)
    {
        //Init ss, es
        ss = -pad; es = ss + N2 - 1;

        //X2 before first samp of X1
        while (es<0 && w<W) { *Y++ = 0.0; es+=str; ++w; }
        ss = es - N2 + 1;

        //X2 overlaps first samp of X1
        while (ss<0 && es<(int)L1 && w<W)
        {
            sm = 0.0; X1 += es;
            for (int n=es; n>=0; n-=dil, X1-=dil, ++X2) { sm += *X1 * *X2; }
            *Y++ = sm;
            X1 += dil - (size_t)es%dil;
            X2 -= 1 + es/(int)dil;
            ss+=str; es+=str; ++w;
        }
        X1 += ss; X2 += L2 - 1u;

        if (N2>(int)L1) //X1 fully within X2
        {
            while (ss<0 && w<W)
            {
                sm = 0.0;
                for (int n=(int)dil*(-ss/(int)dil); n<(int)L1-ss; n+=dil)
                {
                    sm += X1[n+ss] * X2[-n/(int)dil];
                }
                *Y++ = sm;
                ss+=str; ++w;
            }
            X1 += ss;
            es = ss + N2 - 1;
        }
        else            //X2 fully within X1
        {        
            while (es<(int)L1 && w<W)
            {
                sm = 0.0;
                for (size_t l=L2; l>0u; --l, X1+=dil, --X2) { sm += *X1 * *X2; }
                *Y++ = sm;
                X1 += inc; X2 += L2;
                es+=str; ++w;
            }
            ss = es - N2 + 1;
        }

        //X2 overlaps end samp of X1
        while (ss<(int)L1 && w<W)
        {
            size_t l = 0u; sm = 0.0;
            for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=dil, --X2) { sm += *X1 * *X2; }
            *Y++ = sm;
            X1 += (int)str - (int)(dil*l); X2 += l;
            ss+=str; es+=str; ++w;
        }

        //X2 past end samp of X1
        while (w<W) { *Y++ = 0.0; ++w; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L1, G = V/B;

        for (size_t g=G; g>0u; --g, X1+=B*(L1-1u), Y+=B*(W-1u))
        {
            for (size_t b=B; b>0u; --b, ++X1, Y-=K*W-1u)
            {
                //Init ss, es, w
                ss = -pad; es = ss + N2 - 1; w = 0u;

                //X2 before first samp of X1
                while (es<0 && w<W) { *Y = 0.0; Y+=K; es+=str; ++w; }
                ss = es - N2 + 1;

                //X2 overlaps first samp of X1
                while (ss<0 && es<(int)L1 && w<W)
                {
                    sm = 0.0; X1 += es*(int)K;
                    for (int n=es; n>=0; n-=dil, X1-=dil*K, ++X2) { sm += *X1 * *X2; }
                    *Y = sm; Y += K;
                    X1 += K*(dil-(size_t)es%dil);
                    X2 -= 1 + es/(int)dil;
                    ss+=str; es+=str; ++w;
                }
                X1 += ss*(int)K; X2 += L2 - 1u;

                if (N2>(int)L1) //X1 fully within X2
                {
                    while (ss<0 && w<W)
                    {
                        sm = 0.0;
                        for (int n=(int)dil*(-ss/(int)dil); n<(int)L1-ss; n+=dil)
                        {
                            sm += X1[(int)K*(n+ss)] * X2[-n/(int)dil];
                        }
                        *Y = sm; Y += K;
                        ss+=str; ++w;
                    }
                    X1 += ss*(int)K;
                    es = ss + N2 - 1;
                }
                else            //X2 fully within X1
                {        
                    while (es<(int)L1 && w<W)
                    {
                        sm = 0.0;
                        for (size_t l=L2; l>0u; --l, X1+=dil*K, --X2) { sm += *X1 * *X2; }
                        *Y = sm; Y += K;
                        X1 += inc*(int)K; X2 += L2;
                        es+=str; ++w;
                    }
                    ss = es - N2 + 1;
                }

                //X2 overlaps end samp of X1
                while (ss<(int)L1 && w<W)
                {
                    size_t l = 0u; sm = 0.0;
                    for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=dil*K, --X2) { sm += *X1 * *X2; }
                    *Y = sm; Y += K;
                    X1 += (int)K*((int)str-(int)(dil*l));
                    X2 += l;
                    ss+=str; es+=str; ++w;
                }

                //X2 past end samp of X1
                while (w<W) { *Y = 0.0; Y += K; ++w; }

                //Reset X1, X2
                X1 -= (int)K*ss; X2 -= L2 - 1u;
            }
        }
    }

    return 0;
}


int conv1d_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int pad, const size_t str, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv1d_c: dim must be in [0 3]\n"); return 1; }
    if (str<1u) { fprintf(stderr,"error in conv1d_c: str (stride) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in conv1d_c: dil (dilation) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in conv1d_c: L2 (length of X2) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in conv1d_c: N (total length of X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in conv1d_c: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if ((int)L1+2*pad<=(int)(dil*(L2-1u))) { fprintf(stderr,"error in conv1d_c: L1+2*pad must be > dil*(L2-1)\n"); return 1; }

    const int N1 = (int)L1 + 2*pad;             //full length of vecs in X1 including padding
    const int N2 = (int)(dil*(L2-1u)) + 1;      //full length of vec in X2 including dilation
    const int inc = (int)str - (int)(dil*L2);   //fixed increment for X1 below
    size_t w=0u, W;                             //current frame and total frames
    int ss, es;                                 //current start-samp, end-samp
    float smr, smi;                             //intermediate sums

    //Set W (L_out)
    W = 1u + (size_t)(N1-N2)/str;

    if (W==0u) {}
    else if (L1==N)
    {
        //Init ss, es
        ss = -pad; es = ss + N2 - 1;

        //X2 before first samp of X1
        while (es<0 && w<W) { *Y++ = 0.0f; *Y++ = 0.0f; es+=str; ++w; }
        ss = es - N2 + 1;

        //X2 overlaps first samp of X1
        while (ss<0 && es<(int)L1 && w<W)
        {
            smr = smi = 0.0f; X1 += 2*es;
            for (int n=es; n>=0; n-=dil, X1-=2u*dil, X2+=2)
            {
                smr += *X1**X2 - *(X1+1)**(X2+1);
                smi += *X1**(X2+1) + *(X1+1)**X2;
            }
            *Y++ = smr; *Y++ = smi;
            X1 += 2u*(dil - (size_t)es%dil);
            X2 -= 2*(1 + es/(int)dil);
            ss+=str; es+=str; ++w;
        }
        X1 += 2*ss; X2 += 2u*(L2-1u);

        if (N2>(int)L1) //X1 fully within X2
        {
            while (ss<0 && w<W)
            {
                smr = smi = 0.0f;
                for (int n=(int)dil*(-ss/(int)dil); n<(int)L1-ss; n+=dil)
                {
                    smr += X1[2*(n+ss)]*X2[-2*(n/(int)dil)] - X1[2*(n+ss)+1]*X2[-2*(n/(int)dil)+1];
                    smi += X1[2*(n+ss)]*X2[-2*(n/(int)dil)+1] + X1[2*(n+ss)+1]*X2[-2*(n/(int)dil)];
                }
                *Y++ = smr; *Y++ = smi;
                ss+=str; ++w;
            }
            X1 += 2*ss;
            es = ss + N2 - 1;
        }
        else            //X2 fully within X1
        {        
            while (es<(int)L1 && w<W)
            {
                smr = smi = 0.0f;
                for (size_t l=L2; l>0u; --l, X1+=2u*dil, X2-=2)
                {
                    smr += *X1**X2 - *(X1+1)**(X2+1);
                    smi += *X1**(X2+1) + *(X1+1)**X2;
                }
                *Y++ = smr; *Y++ = smi;
                X1 += 2*inc; X2 += 2u*L2;
                es+=str; ++w;
            }
            ss = es - N2 + 1;
        }

        //X2 overlaps end samp of X1
        while (ss<(int)L1 && w<W)
        {
            size_t l = 0u; smr = smi = 0.0f;
            for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=2u*dil, X2-=2)
            {
                smr += *X1**X2 - *(X1+1)**(X2+1);
                smi += *X1**(X2+1) + *(X1+1)**X2;
            }
            *Y++ = smr; *Y++ = smi;
            X1 += 2*((int)str-(int)(dil*l));
            X2 += 2u*l;
            ss+=str; es+=str; ++w;
        }

        //X2 past end samp of X1
        while (w<W) { *Y++ = 0.0f; *Y++ = 0.0f; ++w; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L1, G = V/B;

        for (size_t g=G; g>0u; --g, X1+=2u*B*(L1-1u), Y+=2u*B*(W-1u))
        {
            for (size_t b=B; b>0u; --b, X1+=2, Y-=2u*K*W-2u)
            {
                //Init ss, es, w
                ss = -pad; es = ss + N2 - 1; w = 0u;

                //X2 before first samp of X1
                while (es<0 && w<W) { *Y = *(Y+1) = 0.0f; Y+=2u*K; es+=str; ++w; }
                ss = es - N2 + 1;

                //X2 overlaps first samp of X1
                while (ss<0 && es<(int)L1 && w<W)
                {
                    smr = smi = 0.0f; X1 += 2*es*(int)K;
                    for (int n=es; n>=0; n-=dil, X1-=2u*dil*K, X2+=2)
                    {
                        smr += *X1**X2 - *(X1+1)**(X2+1);
                        smi += *X1**(X2+1) + *(X1+1)**X2;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    X1 += 2u*K*(dil-(size_t)es%dil);
                    X2 -= 2*(1+es/(int)dil);
                    ss+=str; es+=str; ++w;
                }
                X1 += 2*ss*(int)K; X2 += 2u*(L2-1u);

                if (N2>(int)L1) //X1 fully within X2
                {
                    while (ss<0 && w<W)
                    {
                        smr = smi = 0.0f;
                        for (int n=(int)dil*(-ss/(int)dil); n<(int)L1-ss; n+=dil)
                        {
                            smr += X1[2*(int)K*(n+ss)]*X2[-2*(n/(int)dil)] - X1[2*(int)K*(n+ss)+1]*X2[-2*(n/(int)dil)+1];
                            smi += X1[2*(int)K*(n+ss)]*X2[-2*(n/(int)dil)+1] + X1[2*(int)K*(n+ss)+1]*X2[-2*(n/(int)dil)];
                        }
                        *Y = smr; *(Y+1) = smi; Y += 2u*K;
                        ss+=str; ++w;
                    }
                    X1 += 2*ss*(int)K;
                    es = ss + N2 - 1;
                }
                else            //X2 fully within X1
                {        
                    while (es<(int)L1 && w<W)
                    {
                        smr = smi = 0.0f;
                        for (size_t l=L2; l>0u; --l, X1+=2u*dil*K, X2-=2)
                        {
                            smr += *X1**X2 - *(X1+1)**(X2+1);
                            smi += *X1**(X2+1) + *(X1+1)**X2;
                        }
                        *Y = smr; *(Y+1) = smi; Y += 2u*K;
                        X1 += 2*inc*(int)K; X2 += 2u*L2;
                        es+=str; ++w;
                    }
                    ss = es - N2 + 1;
                }

                //X2 overlaps end samp of X1
                while (ss<(int)L1 && w<W)
                {
                    size_t l = 0u; smr = smi = 0.0f;
                    for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=2u*dil*K, X2-=2)
                    {
                        smr += *X1**X2 - *(X1+1)**(X2+1);
                        smi += *X1**(X2+1) + *(X1+1)**X2;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    X1 += 2*(int)K*((int)str-(int)(dil*l));
                    X2 += 2u*l;
                    ss+=str; es+=str; ++w;
                }

                //X2 past end samp of X1
                while (w<W) { *Y = *(Y+1) = 0.0f; Y += 2u*K; ++w; }

                //Reset X1, X2
                X1 -= 2*(int)K*ss; X2 -= 2u*L2 - 2u;
            }
        }
    }

    return 0;
}


int conv1d_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int pad, const size_t str, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv1d_z: dim must be in [0 3]\n"); return 1; }
    if (str<1u) { fprintf(stderr,"error in conv1d_z: str (stride) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in conv1d_z: dil (dilation) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in conv1d_z: L2 (length of X2) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in conv1d_z: N (total length of X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in conv1d_z: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if ((int)L1+2*pad<=(int)(dil*(L2-1u))) { fprintf(stderr,"error in conv1d_z: L1+2*pad must be > dil*(L2-1)\n"); return 1; }

    const int N1 = (int)L1 + 2*pad;             //full length of vecs in X1 including padding
    const int N2 = (int)(dil*(L2-1u)) + 1;      //full length of vec in X2 including dilation
    const int inc = (int)str - (int)(dil*L2);   //fixed increment for X1 below
    size_t w=0u, W;                             //current frame and total frames
    int ss, es;                                 //current start-samp, end-samp
    double smr, smi;                            //intermediate sums

    //Set W (L_out)
    W = 1u + (size_t)(N1-N2)/str;

    if (W==0u) {}
    else if (L1==N)
    {
        //Init ss, es
        ss = -pad; es = ss + N2 - 1;

        //X2 before first samp of X1
        while (es<0 && w<W) { *Y++ = 0.0; *Y++ = 0.0; es+=str; ++w; }
        ss = es - N2 + 1;

        //X2 overlaps first samp of X1
        while (ss<0 && es<(int)L1 && w<W)
        {
            smr = smi = 0.0; X1 += 2*es;
            for (int n=es; n>=0; n-=dil, X1-=2u*dil, X2+=2)
            {
                smr += *X1**X2 - *(X1+1)**(X2+1);
                smi += *X1**(X2+1) + *(X1+1)**X2;
            }
            *Y++ = smr; *Y++ = smi;
            X1 += 2u*(dil - (size_t)es%dil);
            X2 -= 2*(1 + es/(int)dil);
            ss+=str; es+=str; ++w;
        }
        X1 += 2*ss; X2 += 2u*(L2-1u);

        if (N2>(int)L1) //X1 fully within X2
        {
            while (ss<0 && w<W)
            {
                smr = smi = 0.0;
                for (int n=(int)dil*(-ss/(int)dil); n<(int)L1-ss; n+=dil)
                {
                    smr += X1[2*(n+ss)]*X2[-2*(n/(int)dil)] - X1[2*(n+ss)+1]*X2[-2*(n/(int)dil)+1];
                    smi += X1[2*(n+ss)]*X2[-2*(n/(int)dil)+1] + X1[2*(n+ss)+1]*X2[-2*(n/(int)dil)];
                }
                *Y++ = smr; *Y++ = smi;
                ss+=str; ++w;
            }
            X1 += 2*ss;
            es = ss + N2 - 1;
        }
        else            //X2 fully within X1
        {        
            while (es<(int)L1 && w<W)
            {
                smr = smi = 0.0;
                for (size_t l=L2; l>0u; --l, X1+=2u*dil, X2-=2)
                {
                    smr += *X1**X2 - *(X1+1)**(X2+1);
                    smi += *X1**(X2+1) + *(X1+1)**X2;
                }
                *Y++ = smr; *Y++ = smi;
                X1 += 2*inc; X2 += 2u*L2;
                es+=str; ++w;
            }
            ss = es - N2 + 1;
        }

        //X2 overlaps end samp of X1
        while (ss<(int)L1 && w<W)
        {
            size_t l = 0u; smr = smi = 0.0;
            for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=2u*dil, X2-=2)
            {
                smr += *X1**X2 - *(X1+1)**(X2+1);
                smi += *X1**(X2+1) + *(X1+1)**X2;
            }
            *Y++ = smr; *Y++ = smi;
            X1 += 2*((int)str-(int)(dil*l));
            X2 += 2u*l;
            ss+=str; es+=str; ++w;
        }

        //X2 past end samp of X1
        while (w<W) { *Y++ = 0.0; *Y++ = 0.0; ++w; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L1, G = V/B;

        for (size_t g=G; g>0u; --g, X1+=2u*B*(L1-1u), Y+=2u*B*(W-1u))
        {
            for (size_t b=B; b>0u; --b, X1+=2, Y-=2u*K*W-2u)
            {
                //Init ss, es, w
                ss = -pad; es = ss + N2 - 1; w = 0u;

                //X2 before first samp of X1
                while (es<0 && w<W) { *Y = *(Y+1) = 0.0; Y+=2u*K; es+=str; ++w; }
                ss = es - N2 + 1;

                //X2 overlaps first samp of X1
                while (ss<0 && es<(int)L1 && w<W)
                {
                    smr = smi = 0.0; X1 += 2*es*(int)K;
                    for (int n=es; n>=0; n-=dil, X1-=2u*dil*K, X2+=2)
                    {
                        smr += *X1**X2 - *(X1+1)**(X2+1);
                        smi += *X1**(X2+1) + *(X1+1)**X2;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    X1 += 2u*K*(dil-(size_t)es%dil);
                    X2 -= 2*(1+es/(int)dil);
                    ss+=str; es+=str; ++w;
                }
                X1 += 2*ss*(int)K; X2 += 2u*(L2-1u);

                if (N2>(int)L1) //X1 fully within X2
                {
                    while (ss<0 && w<W)
                    {
                        smr = smi = 0.0;
                        for (int n=(int)dil*(-ss/(int)dil); n<(int)L1-ss; n+=dil)
                        {
                            smr += X1[2*(int)K*(n+ss)]*X2[-2*(n/(int)dil)] - X1[2*(int)K*(n+ss)+1]*X2[-2*(n/(int)dil)+1];
                            smi += X1[2*(int)K*(n+ss)]*X2[-2*(n/(int)dil)+1] + X1[2*(int)K*(n+ss)+1]*X2[-2*(n/(int)dil)];
                        }
                        *Y = smr; *(Y+1) = smi; Y += 2u*K;
                        ss+=str; ++w;
                    }
                    X1 += 2*ss*(int)K;
                    es = ss + N2 - 1;
                }
                else            //X2 fully within X1
                {        
                    while (es<(int)L1 && w<W)
                    {
                        smr = smi = 0.0;
                        for (size_t l=L2; l>0u; --l, X1+=2u*dil*K, X2-=2)
                        {
                            smr += *X1**X2 - *(X1+1)**(X2+1);
                            smi += *X1**(X2+1) + *(X1+1)**X2;
                        }
                        *Y = smr; *(Y+1) = smi; Y += 2u*K;
                        X1 += 2*inc*(int)K; X2 += 2u*L2;
                        es+=str; ++w;
                    }
                    ss = es - N2 + 1;
                }

                //X2 overlaps end samp of X1
                while (ss<(int)L1 && w<W)
                {
                    size_t l = 0u; smr = smi = 0.0;
                    for (int n=ss; n<(int)L1; n+=dil, ++l, X1+=2u*dil*K, X2-=2)
                    {
                        smr += *X1**X2 - *(X1+1)**(X2+1);
                        smi += *X1**(X2+1) + *(X1+1)**X2;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    X1 += 2*(int)K*((int)str-(int)(dil*l));
                    X2 += 2u*l;
                    ss+=str; es+=str; ++w;
                }

                //X2 past end samp of X1
                while (w<W) { *Y = *(Y+1) = 0.0; Y += 2u*K; ++w; }

                //Reset X1, X2
                X1 -= 2*(int)K*ss; X2 -= 2u*L2 - 2u;
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
