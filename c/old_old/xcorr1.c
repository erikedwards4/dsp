//1D cross-correlation of each vector in X1 by X2.
//Each vector in X1 has length L1. X2 has length L2 (kernel_size).

//This is identical to conv1, except X2 is NOT flipped.

//The framing is controlled by inputs str, es0, Ly (stride, 1st start sample, length of vecs in Y).

//This is like xcorr1d in allowing stride (str), but doesn't allow padding and dilation.
//This is more convenient in many use cases which are BLAS-like (need stride, start, and num only).
//Unlike BLAS, there is no stride for the output Y, but that could easily be added.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int xcorr1_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int es0, const size_t str, const size_t Ly, const size_t dim);
int xcorr1_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int es0, const size_t str, const size_t Ly, const size_t dim);
int xcorr1_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int es0, const size_t str, const size_t Ly, const size_t dim);
int xcorr1_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int es0, const size_t str, const size_t Ly, const size_t dim);


int xcorr1_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int es0, const size_t str, const size_t Ly, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in xcorr1_s: dim must be in [0 3]\n"); return 1; }
    if (str<1u) { fprintf(stderr,"error in xcorr1_s: str (stride) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in xcorr1_s: L2 (length of X2) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in xcorr1_s: N (total length of X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in xcorr1_s: L1 (length of vecs in X1) must be positive\n"); return 1; }

    const int inc = (int)str - (int)L2;   //fixed increment for X1 below
    size_t w = 0u;                        //current frame (in [0 Ly-1])
    int ss, es = es0;                     //current start-samp, end-samp
    float sm;                             //intermediate sum

    //Don't flip X2
    X2 += L2 - 1u;

    if (Ly==0u) {}
    else if (L1==N)
    {
        //X2 before first samp of X1
        while (es<0 && w<Ly) { *Y++ = 0.0f; es+=str; ++w; }
        ss = es - L2 + 1;

        //X2 overlaps first samp of X1
        while (ss<0 && es<(int)L1 && w<Ly)
        {
            sm = 0.0f; X1 += es;
            for (int n=es; n>0; --n, --X1, --X2) { sm += *X1 * *X2; }
            *Y++ = sm + *X1**X2;
            X2 += es;
            ss+=str; es+=str; ++w;
        }
        X1 += ss; X2 -= L2 - 1u;

        if (L2>L1)  //X1 fully within X2
        {
            while (ss<0 && w<Ly)
            {
                sm = 0.0f;
                for (size_t l=L1; l>0u; --l, ++X1, ++X2) { sm += *X1 * *X2; }
                *Y++ = sm;
                X1 -= L1; X2 -= L1+1u;
                ss+=str; ++w;
            }
            es = ss + L2 - 1;
        }
        else        //X2 fully within X1
        {
            while (es<(int)L1 && w<Ly)
            {
                sm = 0.0f;
                for (size_t l=L2; l>0u; --l, ++X1, ++X2) { sm += *X1 * *X2; }
                *Y++ = sm;
                X1 += inc; X2 -= L2;
                es+=str; ++w;
            }
            ss = es - L2 + 1;
        }

        //X2 overlaps end samp of X1
        while (ss<(int)L1 && w<Ly)
        {
            sm = 0.0f;
            for (int n=ss; n<(int)L1; ++n, ++X1, ++X2) { sm += *X1 * *X2; }
            *Y++ = sm;
            X1 += (int)str - (int)L1 + ss; X2 -= (int)L1 - ss;
            ss+=str; ++w;
        }

        //X2 past end samp of X1
        while (w<Ly) { *Y++ = 0.0f; ++w; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L1, G = V/B;

        for (size_t g=G; g>0u; --g, X1+=B*(L1-1u), Y+=B*(Ly-1u))
        {
            for (size_t b=B; b>0u; --b, ++X1, Y-=K*Ly-1u)
            {
                //Reset framing
                es = es0; w = 0u;

                //X2 before first samp of X1
                while (es<0 && w<Ly) { *Y = 0.0f; Y+=K; es+=str; ++w; }
                ss = es - L2 + 1;

                //X2 overlaps first samp of X1
                while (ss<0 && es<(int)L1 && w<Ly)
                {
                    sm = 0.0f; X1 += es*(int)K;
                    for (int n=es; n>0; --n, X1-=K, --X2) { sm += *X1 * *X2; }
                    *Y = sm + *X1**X2; Y += K;
                    X2 += es;
                    ss+=str; es+=str; ++w;
                }
                X1 += ss*(int)K; X2 -= L2 - 1u;

                if (L2>L1) //X1 fully within X2
                {
                    while (ss<0 && w<Ly)
                    {
                        sm = 0.0f;
                        for (size_t l=L1; l>0u; --l, X1+=K, ++X2) { sm += *X1 * *X2; }
                        *Y = sm; Y += K;
                        X1 -= K*L1; X2 -= L1+1u;
                        ss+=str; ++w;
                    }
                    es = ss + L2 - 1;
                }
                else            //X2 fully within X1
                {        
                    while (es<(int)L1 && w<Ly)
                    {
                        sm = 0.0f;
                        for (size_t l=L2; l>0u; --l, X1+=K, ++X2) { sm += *X1 * *X2; }
                        *Y = sm; Y += K;
                        X1 += inc*(int)K; X2 -= L2;
                        es+=str; ++w;
                    }
                    ss = es - L2 + 1;
                }

                //X2 overlaps end samp of X1
                while (ss<(int)L1 && w<Ly)
                {
                    sm = 0.0f;
                    for (int n=ss; n<(int)L1; ++n, X1+=K, ++X2) { sm += *X1 * *X2; }
                    *Y = sm; Y += K;
                    X1 += (int)K*((int)str-(int)L1+ss);
                    X2 -= (int)L1 - ss;
                    ss+=str; ++w;
                }

                //X2 past end samp of X1
                while (w<Ly) { *Y = 0.0f; Y += K; ++w; }

                //Reset X1, X2
                X1 -= (int)K*ss; X2 += L2 - 1u;
            }
        }
    }

    return 0;
}


int xcorr1_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int es0, const size_t str, const size_t Ly, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in xcorr1_d: dim must be in [0 3]\n"); return 1; }
    if (str<1u) { fprintf(stderr,"error in xcorr1_d: str (stride) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in xcorr1_d: L2 (length of X2) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in xcorr1_d: N (total length of X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in xcorr1_d: L1 (length of vecs in X1) must be positive\n"); return 1; }

    const int inc = (int)str - (int)L2;   //fixed increment for X1 below
    size_t w = 0u;                        //current frame (in [0 Ly-1])
    int ss, es = es0;                     //current start-samp, end-samp
    double sm;                            //intermediate sum

    //Don't flip X2
    X2 += L2 - 1u;

    if (Ly==0u) {}
    else if (L1==N)
    {
        //X2 before first samp of X1
        while (es<0 && w<Ly) { *Y++ = 0.0; es+=str; ++w; }
        ss = es - L2 + 1;

        //X2 overlaps first samp of X1
        while (ss<0 && es<(int)L1 && w<Ly)
        {
            sm = 0.0; X1 += es;
            for (int n=es; n>0; --n, --X1, --X2) { sm += *X1 * *X2; }
            *Y++ = sm + *X1**X2;
            X2 += es;
            ss+=str; es+=str; ++w;
        }
        X1 += ss; X2 -= L2 - 1u;

        if (L2>L1)  //X1 fully within X2
        {
            while (ss<0 && w<Ly)
            {
                sm = 0.0;
                for (size_t l=L1; l>0u; --l, ++X1, ++X2) { sm += *X1 * *X2; }
                *Y++ = sm;
                X1 -= L1; X2 -= L1+1u;
                ss+=str; ++w;
            }
            es = ss + L2 - 1;
        }
        else        //X2 fully within X1
        {        
            while (es<(int)L1 && w<Ly)
            {
                sm = 0.0;
                for (size_t l=L2; l>0u; --l, ++X1, ++X2) { sm += *X1 * *X2; }
                *Y++ = sm;
                X1 += inc; X2 -= L2;
                es+=str; ++w;
            }
            ss = es - L2 + 1;
        }

        //X2 overlaps end samp of X1
        while (ss<(int)L1 && w<Ly)
        {
            sm = 0.0;
            for (int n=ss; n<(int)L1; ++n, ++X1, ++X2) { sm += *X1 * *X2; }
            *Y++ = sm;
            X1 += (int)str - (int)L1 + ss; X2 -= (int)L1 - ss;
            ss+=str; ++w;
        }

        //X2 past end samp of X1
        while (w<Ly) { *Y++ = 0.0; ++w; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L1, G = V/B;

        for (size_t g=G; g>0u; --g, X1+=B*(L1-1u), Y+=B*(Ly-1u))
        {
            for (size_t b=B; b>0u; --b, ++X1, Y-=K*Ly-1u)
            {
                //Reset framing
                es = es0; w = 0u;

                //X2 before first samp of X1
                while (es<0 && w<Ly) { *Y = 0.0; Y+=K; es+=str; ++w; }
                ss = es - L2 + 1;

                //X2 overlaps first samp of X1
                while (ss<0 && es<(int)L1 && w<Ly)
                {
                    sm = 0.0; X1 += es*(int)K;
                    for (int n=es; n>0; --n, X1-=K, --X2) { sm += *X1 * *X2; }
                    *Y = sm + *X1**X2; Y += K;
                    X2 += es;
                    ss+=str; es+=str; ++w;
                }
                X1 += ss*(int)K; X2 -= L2 - 1u;

                if (L2>L1) //X1 fully within X2
                {
                    while (ss<0 && w<Ly)
                    {
                        sm = 0.0;
                        for (size_t l=L1; l>0u; --l, X1+=K, ++X2) { sm += *X1 * *X2; }
                        *Y = sm; Y += K;
                        X1 -= K*L1; X2 -= L1+1u;
                        ss+=str; ++w;
                    }
                    es = ss + L2 - 1;
                }
                else            //X2 fully within X1
                {        
                    while (es<(int)L1 && w<Ly)
                    {
                        sm = 0.0;
                        for (size_t l=L2; l>0u; --l, X1+=K, ++X2) { sm += *X1 * *X2; }
                        *Y = sm; Y += K;
                        X1 += inc*(int)K; X2 -= L2;
                        es+=str; ++w;
                    }
                    ss = es - L2 + 1;
                }

                //X2 overlaps end samp of X1
                while (ss<(int)L1 && w<Ly)
                {
                    sm = 0.0;
                    for (int n=ss; n<(int)L1; ++n, X1+=K, ++X2) { sm += *X1 * *X2; }
                    *Y = sm; Y += K;
                    X1 += (int)K*((int)str-(int)L1+ss);
                    X2 -= (int)L1 - ss;
                    ss+=str; ++w;
                }

                //X2 past end samp of X1
                while (w<Ly) { *Y = 0.0; Y += K; ++w; }

                //Reset X1, X2
                X1 -= (int)K*ss; X2 += L2 - 1u;
            }
        }
    }

    return 0;
}


int xcorr1_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int es0, const size_t str, const size_t Ly, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in xcorr1_c: dim must be in [0 3]\n"); return 1; }
    if (str<1u) { fprintf(stderr,"error in xcorr1_c: str (stride) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in xcorr1_c: L2 (length of X2) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in xcorr1_c: N (total length of X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in xcorr1_c: L1 (length of vecs in X1) must be positive\n"); return 1; }

    const int inc = 2*((int)str-(int)L2);   //fixed increment for X1 below
    size_t w = 0u;                          //current frame (in [0 Ly-1])
    int ss, es =es0;                        //current start-samp, end-samp
    float smr, smi;                         //intermediate sums

    //Don't flip X2
    X2 += 2u*(L2-1u);

    if (Ly==0u) {}
    else if (L1==N)
    {
        //X2 before first samp of X1
        while (es<0 && w<Ly) { *Y++ = 0.0f; *Y++ = 0.0f; es+=str; ++w; }
        ss = es - L2 + 1;

        //X2 overlaps first samp of X1
        while (ss<0 && es<(int)L1 && w<Ly)
        {
            smr = smi = 0.0f; X1 += 2*es;
            for (int n=es; n>0; --n, X1-=2, X2-=2)
            {
                smr += *X1**X2 - *(X1+1)**(X2+1);
                smi += *X1**(X2+1) + *(X1+1)**X2;
            }
            *Y++ = smr + *X1**X2 - *(X1+1)**(X2+1);
            *Y++ = smi + *X1**(X2+1) + *(X1+1)**X2;
            X2 += 2*es;
            ss+=str; es+=str; ++w;
        }
        X1 += 2*ss; X2 -= 2u*(L2-1u);

        if (L2>L1)  //X1 fully within X2
        {
            while (ss<0 && w<Ly)
            {
                smr = smi = 0.0f;
                for (size_t l=L1; l>0u; --l, X1+=2, X2+=2)
                {
                    smr += *X1**X2 - *(X1+1)**(X2+1);
                    smi += *X1**(X2+1) + *(X1+1)**X2;
                }
                *Y++ = smr; *Y++ = smi;
                X1 -= 2u*L1; X2 -= 2u*(L1+1u);
                ss+=str; ++w;
            }
            es = ss + L2 - 1;
        }
        else        //X2 fully within X1
        {
            while (es<(int)L1 && w<Ly)
            {
                smr = smi = 0.0f;
                for (size_t l=L2; l>0u; --l, X1+=2, X2+=2)
                {
                    smr += *X1**X2 - *(X1+1)**(X2+1);
                    smi += *X1**(X2+1) + *(X1+1)**X2;
                }
                *Y++ = smr; *Y++ = smi;
                X1 += inc; X2 -= 2u*L2;
                es+=str; ++w;
            }
            ss = es - L2 + 1;
        }

        //X2 overlaps end samp of X1
        while (ss<(int)L1 && w<Ly)
        {
            smr = smi = 0.0f;
            for (int n=ss; n<(int)L1; ++n, X1+=2, X2+=2)
            {
                smr += *X1**X2 - *(X1+1)**(X2+1);
                smi += *X1**(X2+1) + *(X1+1)**X2;
            }
            *Y++ = smr; *Y++ = smi;
            X1 += 2*((int)str-(int)L1+ss); X2 -= 2*((int)L1-ss);
            ss+=str; ++w;
        }

        //X2 past end samp of X1
        while (w<Ly) { *Y++ = 0.0f; *Y++ = 0.0f; ++w; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L1, G = V/B;

        for (size_t g=G; g>0u; --g, X1+=2u*B*(L1-1u), Y+=2u*B*(Ly-1u))
        {
            for (size_t b=B; b>0u; --b, X1+=2, Y-=2u*K*Ly-2u)
            {
                //Reset framing
                es = es0; w = 0u;

                //X2 before first samp of X1
                while (es<0 && w<Ly) { *Y = 0.0f; *(Y+1) = 0.0f; Y+=2u*K; es+=str; ++w; }
                ss = es - L2 + 1;

                //X2 overlaps first samp of X1
                while (ss<0 && es<(int)L1 && w<Ly)
                {
                    smr = smi = 0.0f; X1 += 2*es*(int)K;
                    for (int n=es; n>0; --n, X1-=2u*K, X2-=2)
                    {
                        smr += *X1**X2 - *(X1+1)**(X2+1);
                        smi += *X1**(X2+1) + *(X1+1)**X2;
                    }
                    *Y = smr + *X1**X2 - *(X1+1)**(X2+1);
                    *(Y+1) = smi + *X1**(X2+1) + *(X1+1)**X2;
                    Y += 2u*K;
                    X2 += 2*es;
                    ss+=str; es+=str; ++w;
                }
                X1 += 2*ss*(int)K; X2 -= 2u*(L2-1u);

                if (L2>L1) //X1 fully within X2
                {
                    while (ss<0 && w<Ly)
                    {
                        smr = smi = 0.0f;
                        for (size_t l=L1; l>0u; --l, X1+=2u*K, X2+=2)
                        {
                            smr += *X1**X2 - *(X1+1)**(X2+1);
                            smi += *X1**(X2+1) + *(X1+1)**X2;
                        }
                        *Y = smr; *(Y+1) = smi; Y += 2u*K;
                        X1 -= 2u*K*L1; X2 -= 2u*(L1+1u);
                        ss+=str; ++w;
                    }
                    es = ss + (int)L2 - 1;
                }
                else            //X2 fully within X1
                {
                    while (es<(int)L1 && w<Ly)
                    {
                        smr = smi = 0.0f;
                        for (size_t l=L2; l>0u; --l, X1+=2u*K, X2+=2)
                        {
                            smr += *X1**X2 - *(X1+1)**(X2+1);
                            smi += *X1**(X2+1) + *(X1+1)**X2;
                        }
                        *Y = smr; *(Y+1) = smi; Y += 2u*K;
                        X1 += inc*(int)K; X2 -= 2u*L2;
                        es+=str; ++w;
                    }
                    ss = es - (int)L2 + 1;
                }

                //X2 overlaps end samp of X1
                while (ss<(int)L1 && w<Ly)
                {
                    smr = smi = 0.0f;
                    for (int n=ss; n<(int)L1; ++n, X1+=2u*K, X2+=2)
                    {
                        smr += *X1**X2 - *(X1+1)**(X2+1);
                        smi += *X1**(X2+1) + *(X1+1)**X2;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    X1 += 2*(int)K*((int)str-(int)L1+ss);
                    X2 -= 2*((int)L1-ss);
                    ss+=str; ++w;
                }

                //X2 past end samp of X1
                while (w<Ly) { *Y = 0.0f; *(Y+1) = 0.0f; Y += 2u*K; ++w; }

                //Reset X1, X2
                X1 -= 2*(int)K*ss; X2 += 2u*(L2-1u);
            }
        }
    }

    return 0;
}


int xcorr1_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int es0, const size_t str, const size_t Ly, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in xcorr1_z: dim must be in [0 3]\n"); return 1; }
    if (str<1u) { fprintf(stderr,"error in xcorr1_z: str (stride) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in xcorr1_z: L2 (length of X2) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in xcorr1_z: N (total length of X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in xcorr1_z: L1 (length of vecs in X1) must be positive\n"); return 1; }

    const int inc = 2*((int)str-(int)L2);   //fixed increment for X1 below
    size_t w = 0u;                          //current frame (in [0 Ly-1])
    int ss, es = es0;                       //current start-samp, end-samp
    double smr, smi;                        //intermediate sums

    //Don't flip X2
    X2 += 2u*(L2-1u);

    if (Ly==0u) {}
    else if (L1==N)
    {
        //X2 before first samp of X1
        while (es<0 && w<Ly) { *Y++ = 0.0; *Y++ = 0.0; es+=str; ++w; }
        ss = es - L2 + 1;

        //X2 overlaps first samp of X1
        while (ss<0 && es<(int)L1 && w<Ly)
        {
            smr = smi = 0.0; X1 += 2*es;
            for (int n=es; n>0; --n, X1-=2, X2-=2)
            {
                smr += *X1**X2 - *(X1+1)**(X2+1);
                smi += *X1**(X2+1) + *(X1+1)**X2;
            }
            *Y++ = smr + *X1**X2 - *(X1+1)**(X2+1);
            *Y++ = smi + *X1**(X2+1) + *(X1+1)**X2;
            X2 += 2*es;
            ss+=str; es+=str; ++w;
        }
        X1 += 2*ss; X2 -= 2u*(L2-1u);

        if (L2>L1)  //X1 fully within X2
        {
            while (ss<0 && w<Ly)
            {
                smr = smi = 0.0;
                for (size_t l=L1; l>0u; --l, X1+=2, X2+=2)
                {
                    smr += *X1**X2 - *(X1+1)**(X2+1);
                    smi += *X1**(X2+1) + *(X1+1)**X2;
                }
                *Y++ = smr; *Y++ = smi;
                X1 -= 2u*L1; X2 -= 2u*(L1+1u);
                ss+=str; ++w;
            }
            es = ss + L2 - 1;
        }
        else        //X2 fully within X1
        {
            while (es<(int)L1 && w<Ly)
            {
                smr = smi = 0.0;
                for (size_t l=L2; l>0u; --l, X1+=2, X2+=2)
                {
                    smr += *X1**X2 - *(X1+1)**(X2+1);
                    smi += *X1**(X2+1) + *(X1+1)**X2;
                }
                *Y++ = smr; *Y++ = smi;
                X1 += inc; X2 -= 2u*L2;
                es+=str; ++w;
            }
            ss = es - L2 + 1;
        }

        //X2 overlaps end samp of X1
        while (ss<(int)L1 && w<Ly)
        {
            smr = smi = 0.0;
            for (int n=ss; n<(int)L1; ++n, X1+=2, X2+=2)
            {
                smr += *X1**X2 - *(X1+1)**(X2+1);
                smi += *X1**(X2+1) + *(X1+1)**X2;
            }
            *Y++ = smr; *Y++ = smi;
            X1 += 2*((int)str-(int)L1+ss); X2 -= 2*((int)L1-ss);
            ss+=str; ++w;
        }

        //X2 past end samp of X1
        while (w<Ly) { *Y++ = 0.0; *Y++ = 0.0; ++w; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L1, G = V/B;

        for (size_t g=G; g>0u; --g, X1+=2u*B*(L1-1u), Y+=2u*B*(Ly-1u))
        {
            for (size_t b=B; b>0u; --b, X1+=2, Y-=2u*K*Ly-2u)
            {
                //Reset framing
                es = es0; w = 0u;

                //X2 before first samp of X1
                while (es<0 && w<Ly) { *Y = 0.0; *(Y+1) = 0.0; Y+=2u*K; es+=str; ++w; }
                ss = es - L2 + 1;

                //X2 overlaps first samp of X1
                while (ss<0 && es<(int)L1 && w<Ly)
                {
                    smr = smi = 0.0; X1 += 2*es*(int)K;
                    for (int n=es; n>0; --n, X1-=2u*K, X2-=2)
                    {
                        smr += *X1**X2 - *(X1+1)**(X2+1);
                        smi += *X1**(X2+1) + *(X1+1)**X2;
                    }
                    *Y = smr + *X1**X2 - *(X1+1)**(X2+1);
                    *(Y+1) = smi + *X1**(X2+1) + *(X1+1)**X2;
                    Y += 2u*K;
                    X2 += 2*es;
                    ss+=str; es+=str; ++w;
                }
                X1 += 2*ss*(int)K; X2 -= 2u*(L2-1u);

                if (L2>L1) //X1 fully within X2
                {
                    while (ss<0 && w<Ly)
                    {
                        smr = smi = 0.0;
                        for (size_t l=L1; l>0u; --l, X1+=2u*K, X2+=2)
                        {
                            smr += *X1**X2 - *(X1+1)**(X2+1);
                            smi += *X1**(X2+1) + *(X1+1)**X2;
                        }
                        *Y = smr; *(Y+1) = smi; Y += 2u*K;
                        X1 -= 2u*K*L1; X2 -= 2u*(L1+1u);
                        ss+=str; ++w;
                    }
                    es = ss + (int)L2 - 1;
                }
                else            //X2 fully within X1
                {
                    while (es<(int)L1 && w<Ly)
                    {
                        smr = smi = 0.0;
                        for (size_t l=L2; l>0u; --l, X1+=2u*K, X2+=2)
                        {
                            smr += *X1**X2 - *(X1+1)**(X2+1);
                            smi += *X1**(X2+1) + *(X1+1)**X2;
                        }
                        *Y = smr; *(Y+1) = smi; Y += 2u*K;
                        X1 += inc*(int)K; X2 -= 2u*L2;
                        es+=str; ++w;
                    }
                    ss = es - (int)L2 + 1;
                }

                //X2 overlaps end samp of X1
                while (ss<(int)L1 && w<Ly)
                {
                    smr = smi = 0.0;
                    for (int n=ss; n<(int)L1; ++n, X1+=2u*K, X2+=2)
                    {
                        smr += *X1**X2 - *(X1+1)**(X2+1);
                        smi += *X1**(X2+1) + *(X1+1)**X2;
                    }
                    *Y = smr; *(Y+1) = smi; Y += 2u*K;
                    X1 += 2*(int)K*((int)str-(int)L1+ss);
                    X2 -= 2*((int)L1-ss);
                    ss+=str; ++w;
                }

                //X2 past end samp of X1
                while (w<Ly) { *Y = 0.0; *(Y+1) = 0.0; Y += 2u*K; ++w; }

                //Reset X1, X2
                X1 -= 2*(int)K*ss; X2 += 2u*(L2-1u);
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
