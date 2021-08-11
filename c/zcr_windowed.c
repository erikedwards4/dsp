//This gets ZCs as usual, and then applies window W to get rate.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int zcr_windowed_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const int dim, const int c0, const float stp, const int going);
int zcr_windowed_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const int dim, const int c0, const double stp, const int going);

int zcr_windowed_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const int dim, const int c0, const float stp, const int going)
{
    if (dim>3u) { fprintf(stderr,"error in zcr_windowed_s: dim must be in [0 3]\n"); return 1; }
    if (Lw<1u) { fprintf(stderr,"error in zcr_windowed_s: Lw (winlength) must be positive\n"); return 1; }
    if (stp<FLT_EPSILON) { fprintf(stderr,"error in zcr_windowed_s: stp (step size) must be positive\n"); return 1; }
    if (going!=0 && going!=1 && going!=-1) { fprintf(stderr,"error in zcr_windowed_s: going must be in {-1,0,1}\n"); return 1; }

    const size_t N = R*C*S*H;
    if (N==0u) {fprintf(stderr,"error in zcr_windowed_s: X1 is empty (size of at least one dim is 0)\n"); return 1; }

    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lx<1u) { fprintf(stderr,"error in zcr_windowed_s: Lx (length of vecs in X1) must be positive\n"); return 1; }
    if (Lw>Lx) { fprintf(stderr,"error in zcr_windowed_s: Lw (winlength) must be <= Lx (length of vecs in X1)\n"); return 1; }
    if (W>Lx) { fprintf(stderr,"error in zcr_windowed_s: W must be <= Lx (length of vecs in X1)\n"); return 1; }
    if (c0>(float)(Lx-1u)) { fprintf(stderr,"error in zcr_windowed_s: c0 (center samp of 1st frame) must be < Lx (length of vecs in X1)\n"); return 1; }
    
    const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
    const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
    const size_t V = N/L, G = V/B;
    if (K!=1u || (G!=1u && B!=1u)) { fprintf(stderr,"error in zcr_windowed_s: vecs in X1 must be contiguous in memory\n"); return 1; }

    if (W==0u || V==0u) {}
    else
    {
        const size_t Lpre = Lw/2u;                  //nsamps before center samp
        const size_t Lpost = Lw-Lpre-1u;            //nsamps after center samp
        size_t w = 0u;                              //current frame
        float cc = c0;                              //current exact center-samp
        int cs = (int)roundf(c0);                   //current rounded center-samp
        int ss, es = cs + (int)Lpost;               //current rounded start-samp, end-samp
        int prev_cs;                                //previous rounded center-samp
        int *Z, s, sp, sm;                          //ints for ZCR

        //Allocate
        if (!(Z=(int *)malloc(Lx*sizeof(int)))) { fprintf(stderr,"error in zcr_windowed_s: problem with malloc. "); perror("malloc"); return 1; }
        
        //For each vec in X1
        for (size_t v=0u; v<V; ++v)
        {
            //Windows before first samp
            while (es<0 && w<W)
            {
                *Y++ = 0.0f;
                ++w; cc += stp; cs = (int)roundf(cc);
                es = cs + (int)Lpost;
            }
            ss = cs - (int)Lpre;
            prev_cs = cs;

            //Windows overlapping first samp
            while (ss<0 && w<W)
            {
                sm = 0.0f;
                X2 -= ss;
                for (size_t l=(size_t)(-ss); l<Lw; ++l, ++X1, ++X2) { sm += *X1 * *X2; }
                *Y++ = sm;
                ++w; cc += stp; cs = (int)roundf(cc);
                X1 -= (int)Lw + ss; X2 -= Lw;
                ss = cs - (int)Lpre; prev_cs = cs;
            }
            X1 += ss;
            es = cs + (int)Lpost;

            //Windows fully within sig
            while (es<(int)N && w<W)
            {
                if (going==0)
                {
                    sp = (*X1++<0.0f); *Z++ = sm = 0;
                    for (size_t l=1u; l<L; ++l, ++X, ++Z) { s = (*X1<0.0f); sm += (s!=sp); *Z = sm; sp = s; }
                    for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                }
                else if (going==1)
                {
                    sp = (*X1++>=0.0f); *Z++ = sm = 0;
                    for (size_t l=1u; l<L; ++l, ++X1, ++Z) { s = (*X1>=0.0f); sm += s*(s!=sp); *Z = sm; sp = s; }
                    for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                }
                else if (going==-1)
                {
                    sp = (*X1++<0.0f); *Z++ = sm = 0;
                    for (size_t l=1u; l<L; ++l, ++X1, ++Z) { s = (*X1<0.0f); sm += s*(s!=sp); *Z = sm; sp = s; }
                    for (size_t l=0u; l<L-Lw+1u; ++l) { --Z; *Z -= *(Z-Lw); }
                }
                sm = 0.0f;
                for (size_t l=0u; l<Lw; ++l, ++X1, ++X2) { sm += *X1 * *X2; }
                *Y++ = sm;
                ++w; cc += stp; cs = (int)roundf(cc);
                X1 += cs - prev_cs - (int)Lw; X2 -= Lw;
                es = cs + (int)Lpost; prev_cs = cs;
            }
            ss = cs - (int)Lpre;

            //Windows overlapping last samp
            while (ss<(int)N && w<W)
            {
                sm = 0.0f;
                for (size_t l=0u; l<N-(size_t)ss; ++l, ++X1, ++X2) { sm += *X1 * *X2; }
                *Y++ = sm;
                ++w; cc += stp; cs = (int)roundf(cc);
                X1 += cs - prev_cs - (int)N + ss; X2 += ss - (int)N;
                ss = cs - (int)Lpre; prev_cs = cs;
            }

            //Windows after last samp
            while (w<W)
            {
                *Y++ = 0.0f; ++w;
            }
        }

        free(Z);
    }

    return 0;
}

// int zcr_windowed_d(double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const int dim, const int c0, const double stp, const int going)
// {
//     const double z = 0.0;
//     const int N = R*C;
//     const int T = (dim==0) ? 1 + (int)(floor(((int)R-1-c0)/stp)) : 1 + (int)(floor(((int)C-1-c0)/stp));
//     const int Lpre = L/2;          //nsamps before center samp
//     const int Lpost = L - L/2 - 1; //nsamps after center samp
//     int cs = c0, t = 0;
//     int c, r, l, n = -1;
//     char *Z;
//     double *W2;

//     //Checks
//     if (R<1) { fprintf(stderr,"error in zcr_windowed_d: R (nrows X) must be positive\n"); return 1; }
//     if (C<1) { fprintf(stderr,"error in zcr_windowed_d: C (ncols X) must be positive\n"); return 1; }
//     if (L<1) { fprintf(stderr,"error in zcr_windowed_d: L (winlength) must be positive\n"); return 1; }

//     //Initialize Z
//     if (!(Z=(char *)malloc((size_t)(N)*sizeof(char)))) { fprintf(stderr,"error in zcr_windowed_d: problem with malloc. "); perror("malloc"); return 1; }

//     //Initialize W2 (same as W if L even, else 0.5*W added to itself with lag 1)
//     if (!(W2=(double *)malloc((size_t)(L)*sizeof(double)))) { fprintf(stderr,"error in zcr_windowed_d: problem with malloc. "); perror("malloc"); return 1; }
//     if (L%2) { W2[0] = 0.5*W[0]; for (l=1; l<L; l++) { W2[l] = 0.5*(W[l]+W[l-1]); } }
//     else { for (l=0; l<L; l++) { W2[l] = W[l]; } }

//     //Raw material
//     if (going>0) { while (++n<N) { Z[n] = (X[n]>=0.0); } }
//     else { while (++n<N) { Z[n] = (X[n]<0.0); } }

//     if (dim==0)
//     {
//         cblas_dcopy(T*C,&z,0,&Y[0],1);
//         if (iscolmajor)
//         {
//             if (going) { while (--n>0) { Z[n] *= (Z[n]!=Z[n-1]); } }
//             else { while (--n>0) { Z[n] = (Z[n]!=Z[n-1]); } }
//             while (n<N) { Z[n] = 0; n += R; }
//             while (t<T && cs<=Lpre)
//             {
//                 for (c=0; c<C; c++) { for (l=Lpre-cs; l<L; l++) { if (Z[c*R+cs+l-Lpre]) { Y[t+c*T] += W2[l]; } } }
//                 t++; cs = (int)(round(t*stp)) + c0;
//             }
//             while (t<T && cs<R-Lpost)
//             {
//                 for (c=0; c<C; c++) { for (l=0; l<L; l++) { if (Z[c*R+cs+l-Lpre]) { Y[t+c*T] += W2[l]; } } }
//                 t++; cs = (int)(round(t*stp)) + c0;
//             }
//             while (t<T)
//             {
//                 for (c=0; c<C; c++) { for (l=0; l<Lpre+R-cs; l++) { if (Z[c*R+cs+l-Lpre]) { Y[t+c*T] += W2[l]; } } }
//                 t++; cs = (int)(round(t*stp)) + c0;
//             }
//         }
//         else
//         {
//             if (going) { while (--n>=C) { Z[n] *= (Z[n]!=Z[n-C]); } }
//             else { while (--n>=C) { Z[n] = (Z[n]!=Z[n-C]); } }
//             while (n>=0) { Z[n--] = 0; }
//             while (t<T && cs<=Lpre)
//             {
//                 for (c=0; c<C; c++) { for (l=Lpre-cs; l<L; l++) { if (Z[c+(cs+l-Lpre)*C]) { Y[c+t*C] += W2[l]; } } }
//                 t++; cs = (int)(round(t*stp)) + c0;
//             }
//             while (t<T && cs<R-Lpost)
//             {
//                 for (c=0; c<C; c++) { for (l=0; l<L; l++) { if (Z[c+(cs+l-Lpre)*C]) { Y[c+t*C] += W2[l]; } } }
//                 t++; cs = (int)(round(t*stp)) + c0;
//             }
//             while (t<T)
//             {
//                 for (c=0; c<C; c++) { for (l=0; l<Lpre+R-cs; l++) { if (Z[c+(cs+l-Lpre)*C]) { Y[c+t*C] += W2[l]; } } }
//                 t++; cs = (int)(round(t*stp)) + c0;
//             }
//         }
//     }
//     else if (dim==1)
//     {
//         cblas_dcopy(R*T,&z,0,&Y[0],1);
//         if (iscolmajor)
//         {
//             if (going) { while (--n>=R) { Z[n] *= (Z[n]!=Z[n-R]); } }
//             else { while (--n>=R) { Z[n] = (Z[n]!=Z[n-R]); } }
//             while (n>=0) { Z[n--] = 0; }
//             while (t<T && cs<=Lpre)
//             {
//                 for (r=0; r<R; r++) { for (l=Lpre-cs; l<L; l++) { if (Z[r+(cs+l-Lpre)*R]) { Y[r+t*R] += W2[l]; } } }
//                 t++; cs = (int)(round(t*stp)) + c0;
//             }
//             while (t<T && cs<R-Lpost)
//             {
//                 for (r=0; r<R; r++) { for (l=0; l<L; l++) { if (Z[r+(cs+l-Lpre)*R]) { Y[r+t*R] += W2[l]; } } }
//                 t++; cs = (int)(round(t*stp)) + c0;
//             }
//             while (t<T)
//             {
//                 for (r=0; r<R; r++) { for (l=0; l<Lpre+C-cs; l++) { if (Z[r+(cs+l-Lpre)*R]) { Y[r+t*R] += W2[l]; } } }
//                 t++; cs = (int)(round(t*stp)) + c0;
//             }
//         }
//         else
//         {
//             if (going) { while (--n>0) { Z[n] *= (Z[n]!=Z[n-1]); } }
//             else { while (--n>0) { Z[n] = (Z[n]!=Z[n-1]); } }
//             while (n<N) { Z[n] = 0; n += C; }
//             while (t<T && cs<=Lpre)
//             {
//                 for (r=0; r<R; r++) { for (l=Lpre-cs; l<L; l++) { if (Z[r*C+cs+l-Lpre]) { Y[t+r*T] += W2[l]; } } }
//                 t++; cs = (int)(round(t*stp)) + c0;
//             }
//             while (t<T && cs<R-Lpost)
//             {
//                 for (r=0; r<R; r++) { for (l=0; l<L; l++) { if (Z[r*C+cs+l-Lpre]) { Y[t+r*T] += W2[l]; } } }
//                 t++; cs = (int)(round(t*stp)) + c0;
//             }
//             while (t<T)
//             {
//                 for (r=0; r<R; r++) { for (l=0; l<Lpre+C-cs; l++) { if (Z[r*C+cs+l-Lpre]) { Y[t+r*T] += W2[l]; } } }
//                 t++; cs = (int)(round(t*stp)) + c0;
//             }
//         }
//     }
//     else
//     {
//         fprintf(stderr,"error in zcr_windowed_d: dim must be 0 or 1.\n"); return 1;
//     }

//     free(Z); free(W2);
//     return 0;
// }


#ifdef __cplusplus
}
}
#endif
