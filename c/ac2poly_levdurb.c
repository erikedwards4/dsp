//Gets polynomial params from the autocorrelation (AC) function for each vector in X.
//Uses a Levinson-Durbin recursion from the AC values.

//The complex-valued case here uses multiply and divide with no conjugation.
//This matches the real-valued case (i.e., make a complex-valued input with 0 imaginary parts),
//but does not match Octave. However, I can't even match the very first number with Octave,
//even after trying every conj or no-conj combination, so something else is going on.
//Thus, if really requiring the complex-valued case, reconsider first.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ac2poly_levdurb_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2poly_levdurb_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2poly_levdurb_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2poly_levdurb_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);


int ac2poly_levdurb_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ac2poly_levdurb_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (Lx==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 1.0f; }
    }
    else if (Lx==2u)
    {
        if (Lx==N)
        {
            *Y++ = 1.0f; *Y++ = -*(X+1) / *X;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, X+=2)
                {
                    *Y++ = 1.0f; *Y++ = -*(X+1) / *X;
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=B, Y+=B)
                {
                    for (size_t b=0u; b<B; ++b, ++X, ++Y)
                    {
                        *Y = 1.0f; *(Y+K) = -*(X+K) / *X;
                    }
                }
            }
        }
    }
    else
    {
        const size_t P = Lx - 1u;
        float *A, *Atmp, a, e;
        if (!(Atmp=(float *)malloc((P-1u)*sizeof(float)))) { fprintf(stderr,"error in ac2poly_levdurb_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A=(float *)malloc((P)*sizeof(float)))) { fprintf(stderr,"error in ac2poly_levdurb_s: problem with malloc. "); perror("malloc"); return 1; }
        
        struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

        if (Lx==N)
        {
            *Y++ = 1.0f;
            a = -*(X+1) / *X;
            *Y = a; //*rcs++ = a;
            e = *X++; e += a * *X++;
            for (size_t p=1u; p<P; ++p, X+=p)
            {
                a = *X;
                for (size_t q=0u; q<p; ++q, ++Y) { --X; a += *X * *Y; }
                a /= -e;
                *Y = a; //*rcs++ = a;
                for (size_t q=0u; q<p; ++q, ++A) { --Y; *A = *Y; }
                Y += p;
                for (size_t q=0u; q<p; ++q) { --A; --Y; *Y += a * *A; }
                e *= 1.0f - a*a;
            }

            // *Y++ = 1.0f;
            // a = -*(X+1) / *X;
            // *Y++ = *A++ = a;
            // e = *X++; e += a * *X++;
            // for (size_t p=1u; p<P; ++p, ++X, ++Y, A+=p)
            // {
            //     a = *X; X -= p;
            //     for (size_t q=0u; q<p; ++q, ++X) { --Y; a += *X * *Y; }
            //     a /= -e;
            //     Y += p; *Y = a;
            //     for (size_t q=0u; q<p; ++q) { *--A = *--Y; }
            //     A += p;
            //     for (size_t q=0u; q<p; ++q, ++Y) { --A; *Y += a * *A; }
            //     e *= 1.0f - a*a;
            // }
            // A -= P;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v)
                {
                    a = -*(X+1) / *X;
                    *A++ = 1.0f; *A++ = *Atmp++ = a;
                    *Y++ = -a;
                    e = *X++; e += a * *X++;
                    for (size_t p=2u; p<P; ++p, ++X, ++Y)
                    {
                        a = *X;
                        X -= p - 1u;
                        for (size_t q=1u; q<p; ++q, ++X) { --A; a += *X * *A; }
                        a /= -e; *Y = -a; *(A+p-1u) = a;
                        for (size_t q=1u; q<p; ++q, ++A) { --Atmp; *A += a * *Atmp; }
                        A -= p - 1u;
                        for (size_t q=0u; q<p; ++q, ++A, ++Atmp) { *Atmp = *A; }
                        e *= 1.0f - a*a;
                    }
                    a = *X;
                    X -= P - 1u;
                    for (size_t q=1u; q<P; ++q, ++X) { --A; a += *X * *A; }
                    Atmp -= P - 1u; --A;
                    *Y++ = a/e; ++X;
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=B*(P-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=K*Lx-K-1u, Y-=K*P-K-1u)
                    {
                        a = -*(X+K) / *X;
                        *A++ = 1.0f; *A++ = *Atmp++ = a;
                        *Y = -a; Y += K;
                        e = *X; X += K; e += a * *X; X += K;
                        for (size_t p=2u; p<P; ++p, X+=K, Y+=K)
                        {
                            a = *X;
                            X -= K*(p-1u);
                            for (size_t q=1u; q<p; ++q, X+=K) { --A; a += *X * *A; }
                            a /= -e; *Y = -a; *(A+p-1u) = a;
                            for (size_t q=1u; q<p; ++q, ++A) { --Atmp; *A += a * *Atmp; }
                            A -= p - 1u;
                            for (size_t q=0u; q<p; ++q, ++A, ++Atmp) { *Atmp = *A; }
                            e *= 1.0f - a*a;
                        }
                        a = *X;
                        X -= K*(P-1u);
                        for (size_t q=1u; q<P; ++q, X+=K) { --A; a += *X * *A; }
                        Atmp -= P - 1u; --A;
                        *Y = a/e;
                    }
                }
            }
        }
        clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
        free(A); free(Atmp);
    }
	
	return 0;
}

int ac2poly_levdurb_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
	if (dim>3u) { fprintf(stderr,"error in ac2poly_levdurb_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u || Lx<2u) {}
    else
    {
        const size_t P = Lx - 1u;
        double *A, *Atmp, a, e;
        if (!(A=(double *)malloc(P*sizeof(double)))) { fprintf(stderr,"error in ac2poly_levdurb_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Atmp=(double *)malloc((P-1u)*sizeof(double)))) { fprintf(stderr,"error in ac2poly_levdurb_d: problem with malloc. "); perror("malloc"); return 1; }
        
        if (Lx==N)
        {
            a = -*(X+1) / *X;
            *A++ = 1.0; *A++ = *Atmp++ = a;
            *Y++ = -a;
            e = *X++; e += a * *X++;
            for (size_t p=2u; p<P; ++p, ++X, ++Y)
            {
                a = *X;
                X -= p - 1u;
                for (size_t q=1u; q<p; ++q, ++X) { --A; a += *X * *A; }
                a /= -e; *Y = -a; *(A+p-1u) = a;
                for (size_t q=1u; q<p; ++q, ++A) { --Atmp; *A += a * *Atmp; }
                A -= p - 1u;
                for (size_t q=0u; q<p; ++q, ++A, ++Atmp) { *Atmp = *A; }
                e *= 1.0 - a*a;
            }
            a = *X;
            X -= P - 1u;
            for (size_t q=1u; q<P; ++q, ++X) { --A; a += *X * *A; }
            Atmp -= P - 1u; --A;
            *Y = a/e;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v)
                {
                    a = -*(X+1) / *X;
                    *A++ = 1.0; *A++ = *Atmp++ = a;
                    *Y++ = -a;
                    e = *X++; e += a * *X++;
                    for (size_t p=2u; p<P; ++p, ++X, ++Y)
                    {
                        a = *X;
                        X -= p - 1u;
                        for (size_t q=1u; q<p; ++q, ++X) { --A; a += *X * *A; }
                        a /= -e; *Y = -a; *(A+p-1u) = a;
                        for (size_t q=1u; q<p; ++q, ++A) { --Atmp; *A += a * *Atmp; }
                        A -= p - 1u;
                        for (size_t q=0u; q<p; ++q, ++A, ++Atmp) { *Atmp = *A; }
                        e *= 1.0 - a*a;
                    }
                    a = *X;
                    X -= P - 1u;
                    for (size_t q=1u; q<P; ++q, ++X) { --A; a += *X * *A; }
                    Atmp -= P - 1u; --A;
                    *Y++ = a/e; ++X;
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=B*(P-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=K*Lx-K-1u, Y-=K*P-K-1u)
                    {
                        a = -*(X+K) / *X;
                        *A++ = 1.0; *A++ = *Atmp++ = a;
                        *Y = -a; Y += K;
                        e = *X; X += K; e += a * *X; X += K;
                        for (size_t p=2u; p<P; ++p, X+=K, Y+=K)
                        {
                            a = *X;
                            X -= K*(p-1u);
                            for (size_t q=1u; q<p; ++q, X+=K) { --A; a += *X * *A; }
                            a /= -e; *Y = -a; *(A+p-1u) = a;
                            for (size_t q=1u; q<p; ++q, ++A) { --Atmp; *A += a * *Atmp; }
                            A -= p - 1u;
                            for (size_t q=0u; q<p; ++q, ++A, ++Atmp) { *Atmp = *A; }
                            e *= 1.0 - a*a;
                        }
                        a = *X;
                        X -= K*(P-1u);
                        for (size_t q=1u; q<P; ++q, X+=K) { --A; a += *X * *A; }
                        Atmp -= P - 1u; --A;
                        *Y = a/e;
                    }
                }
            }
        }
        free(A); free(Atmp);
    }
	
	return 0;
}


int ac2poly_levdurb_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ac2poly_levdurb_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u || Lx<2u) {}
    else
    {
        const size_t P = Lx - 1u;
        float *A, *Atmp, ar, ai, aa, aar, aai, er, ei, ea;
        if (!(A=(float *)malloc(2u*P*sizeof(float)))) { fprintf(stderr,"error in ac2poly_levdurb_c: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Atmp=(float *)malloc(2u*(P-1u)*sizeof(float)))) { fprintf(stderr,"error in ac2poly_levdurb_c: problem with malloc. "); perror("malloc"); return 1; }

        if (Lx==N)
        {
            aa = *X**X + *(X+1)**(X+1);
            ar = -(*(X+2)**X+*(X+3)**(X+1)) / aa;
            ai = -(*(X+3)**X-*(X+2)**(X+1)) / aa;
            *A++ = 1.0f; *A++ = 0.0f;
            *A++ = *Atmp++ = ar;
            *A++ = *Atmp++ = ai;
            *Y++ = -ar; *Y++ = -ai;
            er = *X++; ei = *X++;
            er += ar**X - ai**(X+1);
            ei += ar**(X+1) + ai**X;
            X += 2;
            for (size_t p=2u; p<P; ++p, X+=2)
            {
                ar = *X++; ai = *X--;
                X -= 2u*(p-1u);
                for (size_t q=1u; q<p; ++q, X+=2)
                {
                    A -= 2;
                    ar += *A**X - *(A+1)**(X+1);
                    ai += *A**(X+1) + *(A+1)**X;
                }
                ea = er*er + ei*ei;
                aa = (ar*er+ai*ei) / ea;
                ai = (ar*ei-ai*er) / ea;
                ar = -aa;
                *Y++ = -ar; *Y++ = -ai;
                *(A+2u*p-2u) = ar; *(A+2u*p-1u) = ai;
                for (size_t q=1u; q<p; ++q)
                {
                    Atmp -= 2;
                    *A++ += ar**Atmp - ai**(Atmp+1);
                    *A++ = ar**(Atmp+1) + ai**Atmp;
                }
                A -= 2u*(p-1u);
                for (size_t q=0u; q<2u*p; ++q, ++A, ++Atmp) { *Atmp = *A; }
                aar = 1.0f - ar*ar + ai*ai;
                aai = - ar*ai - ar*ai;
                ea = er*aar - ei*aai;
                ei = er*aai + ei*aar;
                er = ea;
            }
            ar = *X; ai = *(X+1);
            X -= 2u*(P-1u);
            for (size_t q=1u; q<P; ++q, X+=2)
            {
                A -= 2;
                ar += *A**X - *(A+1)**(X+1);
                ai += *A**(X+1) + *(A+1)**X;
            }
            Atmp -= 2u*(P-1u); A -= 2;
            ea = er*er + ei*ei;
            *Y++ = (ar*er+ai*ei) / ea;
            *Y++ = (ai*er-ar*ei) / ea;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, X+=2)
                {
                    aa = *X**X + *(X+1)**(X+1);
                    ar = -(*(X+2)**X+*(X+3)**(X+1)) / aa;
                    ai = -(*(X+3)**X-*(X+2)**(X+1)) / aa;
                    *A++ = 1.0f; *A++ = 0.0f;
                    *A++ = *Atmp++ = ar;
                    *A++ = *Atmp++ = ai;
                    *Y++ = -ar; *Y++ = -ai;
                    er = *X++; ei = *X++;
                    er += ar**X - ai**(X+1);
                    ei += ar**(X+1) + ai**X;
                    X += 2;
                    for (size_t p=2u; p<P; ++p, X+=2)
                    {
                        ar = *X++; ai = *X--;
                        X -= 2u*(p-1u);
                        for (size_t q=1u; q<p; ++q, X+=2)
                        {
                            A -= 2;
                            ar += *A**X - *(A+1)**(X+1);
                            ai += *A**(X+1) + *(A+1)**X;
                        }
                        ea = er*er + ei*ei;
                        aa = (ar*er+ai*ei) / ea;
                        ai = (ar*ei-ai*er) / ea;
                        ar = -aa;
                        *Y++ = -ar; *Y++ = -ai;
                        *(A+2u*p-2u) = ar; *(A+2u*p-1u) = ai;
                        for (size_t q=1u; q<p; ++q)
                        {
                            Atmp -= 2;
                            *A++ += ar**Atmp - ai**(Atmp+1);
                            *A++ = ar**(Atmp+1) + ai**Atmp;
                        }
                        A -= 2u*(p-1u);
                        for (size_t q=0u; q<2u*p; ++q, ++A, ++Atmp) { *Atmp = *A; }
                        aar = 1.0f - ar*ar + ai*ai;
                        aai = - ar*ai - ar*ai;
                        ea = er*aar - ei*aai;
                        ei = er*aai + ei*aar;
                        er = ea;
                    }
                    ar = *X; ai = *(X+1);
                    X -= 2u*(P-1u);
                    for (size_t q=1u; q<P; ++q, X+=2)
                    {
                        A -= 2;
                        ar += *A**X - *(A+1)**(X+1);
                        ai += *A**(X+1) + *(A+1)**X;
                    }
                    Atmp -= 2u*(P-1u); A -= 2;
                    ea = er*er + ei*ei;
                    *Y++ = (ar*er+ai*ei) / ea;
                    *Y++ = (ai*er-ar*ei) / ea;
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*(P-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=2u*(K*Lx-K-1u), Y-=2u*(K*P-K-1u))
                    {
                        aa = *X**X + *(X+1)**(X+1);
                        ar = -(*(X+2u*K)**X+*(X+2u*K+1u)**(X+1)) / aa;
                        ai = -(*(X+2u*K+1u)**X-*(X+2u*K)**(X+1)) / aa;
                        *A++ = 1.0f; *A++ = 0.0f;
                        *A++ = *Atmp++ = ar;
                        *A++ = *Atmp++ = ai;
                        *Y = -ar; *(Y+1) = -ai; Y += 2u*K;
                        er = *X; ei = *(X+1); X += 2u*K;
                        er += ar**X - ai**(X+1);
                        ei += ar**(X+1) + ai**X;
                        X += 2u*K;
                        for (size_t p=2u; p<P; ++p, X+=2u*K)
                        {
                            ar = *X++; ai = *X--;
                            X -= 2u*(p-1u)*K;
                            for (size_t q=1u; q<p; ++q, X+=2u*K)
                            {
                                A -= 2;
                                ar += *A**X - *(A+1)**(X+1);
                                ai += *A**(X+1) + *(A+1)**X;
                            }
                            ea = er*er + ei*ei;
                            aa = (ar*er+ai*ei) / ea;
                            ai = (ar*ei-ai*er) / ea;
                            ar = -aa;
                            *Y = -ar; *(Y+1) = -ai; Y += 2u*K;
                            *(A+2u*p-2u) = ar; *(A+2u*p-1u) = ai;
                            for (size_t q=1u; q<p; ++q)
                            {
                                Atmp -= 2;
                                *A++ += ar**Atmp - ai**(Atmp+1);
                                *A++ = ar**(Atmp+1) + ai**Atmp;
                            }
                            A -= 2u*(p-1u);
                            for (size_t q=0u; q<2u*p; ++q, ++A, ++Atmp) { *Atmp = *A; }
                            aar = 1.0f - ar*ar + ai*ai;
                            aai = - ar*ai - ar*ai;
                            ea = er*aar - ei*aai;
                            ei = er*aai + ei*aar;
                            er = ea;
                        }
                        ar = *X; ai = *(X+1);
                        X -= 2u*(P-1u)*K;
                        for (size_t q=1u; q<P; ++q, X+=2u*K)
                        {
                            A -= 2;
                            ar += *A**X - *(A+1)**(X+1);
                            ai += *A**(X+1) + *(A+1)**X;
                        }
                        Atmp -= 2u*(P-1u); A -= 2;
                        ea = er*er + ei*ei;
                        *Y = (ar*er+ai*ei) / ea;
                        *(Y+1) = (ai*er-ar*ei) / ea;
                    }
                }
            }
        }
        free(A); free(Atmp);
    }
	
	return 0;
}


int ac2poly_levdurb_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ac2poly_levdurb_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u || Lx<2u) {}
    else
    {
        const size_t P = Lx - 1u;
        double *A, *Atmp, ar, ai, aa, aar, aai, er, ei, ea;
        if (!(A=(double *)malloc(2u*P*sizeof(double)))) { fprintf(stderr,"error in ac2poly_levdurb_z: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Atmp=(double *)malloc(2u*(P-1u)*sizeof(double)))) { fprintf(stderr,"error in ac2poly_levdurb_z: problem with malloc. "); perror("malloc"); return 1; }

        if (Lx==N)
        {
            aa = *X**X + *(X+1)**(X+1);
            ar = -(*(X+2)**X+*(X+3)**(X+1)) / aa;
            ai = -(*(X+3)**X-*(X+2)**(X+1)) / aa;
            *A++ = 1.0; *A++ = 0.0;
            *A++ = *Atmp++ = ar;
            *A++ = *Atmp++ = ai;
            *Y++ = -ar; *Y++ = -ai;
            er = *X++; ei = *X++;
            er += ar**X - ai**(X+1);
            ei += ar**(X+1) + ai**X;
            X += 2;
            for (size_t p=2u; p<P; ++p, X+=2)
            {
                ar = *X++; ai = *X--;
                X -= 2u*(p-1u);
                for (size_t q=1u; q<p; ++q, X+=2)
                {
                    A -= 2;
                    ar += *A**X - *(A+1)**(X+1);
                    ai += *A**(X+1) + *(A+1)**X;
                }
                ea = er*er + ei*ei;
                aa = (ar*er+ai*ei) / ea;
                ai = (ar*ei-ai*er) / ea;
                ar = -aa;
                *Y++ = -ar; *Y++ = -ai;
                *(A+2u*p-2u) = ar; *(A+2u*p-1u) = ai;
                for (size_t q=1u; q<p; ++q)
                {
                    Atmp -= 2;
                    *A++ += ar**Atmp - ai**(Atmp+1);
                    *A++ = ar**(Atmp+1) + ai**Atmp;
                }
                A -= 2u*(p-1u);
                for (size_t q=0u; q<2u*p; ++q, ++A, ++Atmp) { *Atmp = *A; }
                aar = 1.0 - ar*ar + ai*ai;
                aai = - ar*ai - ar*ai;
                ea = er*aar - ei*aai;
                ei = er*aai + ei*aar;
                er = ea;
            }
            ar = *X; ai = *(X+1);
            X -= 2u*(P-1u);
            for (size_t q=1u; q<P; ++q, X+=2)
            {
                A -= 2;
                ar += *A**X - *(A+1)**(X+1);
                ai += *A**(X+1) + *(A+1)**X;
            }
            Atmp -= 2u*(P-1u); A -= 2;
            ea = er*er + ei*ei;
            *Y++ = (ar*er+ai*ei) / ea;
            *Y++ = (ai*er-ar*ei) / ea;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, X+=2)
                {
                    aa = *X**X + *(X+1)**(X+1);
                    ar = -(*(X+2)**X+*(X+3)**(X+1)) / aa;
                    ai = -(*(X+3)**X-*(X+2)**(X+1)) / aa;
                    *A++ = 1.0; *A++ = 0.0;
                    *A++ = *Atmp++ = ar;
                    *A++ = *Atmp++ = ai;
                    *Y++ = -ar; *Y++ = -ai;
                    er = *X++; ei = *X++;
                    er += ar**X - ai**(X+1);
                    ei += ar**(X+1) + ai**X;
                    X += 2;
                    for (size_t p=2u; p<P; ++p, X+=2)
                    {
                        ar = *X++; ai = *X--;
                        X -= 2u*(p-1u);
                        for (size_t q=1u; q<p; ++q, X+=2)
                        {
                            A -= 2;
                            ar += *A**X - *(A+1)**(X+1);
                            ai += *A**(X+1) + *(A+1)**X;
                        }
                        ea = er*er + ei*ei;
                        aa = (ar*er+ai*ei) / ea;
                        ai = (ar*ei-ai*er) / ea;
                        ar = -aa;
                        *Y++ = -ar; *Y++ = -ai;
                        *(A+2u*p-2u) = ar; *(A+2u*p-1u) = ai;
                        for (size_t q=1u; q<p; ++q)
                        {
                            Atmp -= 2;
                            *A++ += ar**Atmp - ai**(Atmp+1);
                            *A++ = ar**(Atmp+1) + ai**Atmp;
                        }
                        A -= 2u*(p-1u);
                        for (size_t q=0u; q<2u*p; ++q, ++A, ++Atmp) { *Atmp = *A; }
                        aar = 1.0 - ar*ar + ai*ai;
                        aai = - ar*ai - ar*ai;
                        ea = er*aar - ei*aai;
                        ei = er*aai + ei*aar;
                        er = ea;
                    }
                    ar = *X; ai = *(X+1);
                    X -= 2u*(P-1u);
                    for (size_t q=1u; q<P; ++q, X+=2)
                    {
                        A -= 2;
                        ar += *A**X - *(A+1)**(X+1);
                        ai += *A**(X+1) + *(A+1)**X;
                    }
                    Atmp -= 2u*(P-1u); A -= 2;
                    ea = er*er + ei*ei;
                    *Y++ = (ar*er+ai*ei) / ea;
                    *Y++ = (ai*er-ar*ei) / ea;
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*(P-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=2u*(K*Lx-K-1u), Y-=2u*(K*P-K-1u))
                    {
                        aa = *X**X + *(X+1)**(X+1);
                        ar = -(*(X+2u*K)**X+*(X+2u*K+1u)**(X+1)) / aa;
                        ai = -(*(X+2u*K+1u)**X-*(X+2u*K)**(X+1)) / aa;
                        *A++ = 1.0; *A++ = 0.0;
                        *A++ = *Atmp++ = ar;
                        *A++ = *Atmp++ = ai;
                        *Y = -ar; *(Y+1) = -ai; Y += 2u*K;
                        er = *X; ei = *(X+1); X += 2u*K;
                        er += ar**X - ai**(X+1);
                        ei += ar**(X+1) + ai**X;
                        X += 2u*K;
                        for (size_t p=2u; p<P; ++p, X+=2u*K)
                        {
                            ar = *X++; ai = *X--;
                            X -= 2u*(p-1u)*K;
                            for (size_t q=1u; q<p; ++q, X+=2u*K)
                            {
                                A -= 2;
                                ar += *A**X - *(A+1)**(X+1);
                                ai += *A**(X+1) + *(A+1)**X;
                            }
                            ea = er*er + ei*ei;
                            aa = (ar*er+ai*ei) / ea;
                            ai = (ar*ei-ai*er) / ea;
                            ar = -aa;
                            *Y = -ar; *(Y+1) = -ai; Y += 2u*K;
                            *(A+2u*p-2u) = ar; *(A+2u*p-1u) = ai;
                            for (size_t q=1u; q<p; ++q)
                            {
                                Atmp -= 2;
                                *A++ += ar**Atmp - ai**(Atmp+1);
                                *A++ = ar**(Atmp+1) + ai**Atmp;
                            }
                            A -= 2u*(p-1u);
                            for (size_t q=0u; q<2u*p; ++q, ++A, ++Atmp) { *Atmp = *A; }
                            aar = 1.0 - ar*ar + ai*ai;
                            aai = - ar*ai - ar*ai;
                            ea = er*aar - ei*aai;
                            ei = er*aai + ei*aar;
                            er = ea;
                        }
                        ar = *X; ai = *(X+1);
                        X -= 2u*(P-1u)*K;
                        for (size_t q=1u; q<P; ++q, X+=2u*K)
                        {
                            A -= 2;
                            ar += *A**X - *(A+1)**(X+1);
                            ai += *A**(X+1) + *(A+1)**X;
                        }
                        Atmp -= 2u*(P-1u); A -= 2;
                        ea = er*er + ei*ei;
                        *Y = (ar*er+ai*ei) / ea;
                        *(Y+1) = (ai*er-ar*ei) / ea;
                    }
                }
            }
        }
        free(A); free(Atmp);
    }
	
	return 0;
}


#ifdef __cplusplus
}
}
#endif
