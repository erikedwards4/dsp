//These do a quick linear interpolation.
//These assume no extrapolation, and assume that Xi is sorted ascending.
//Basic extrapolation is supported just in case --
//the first/last members of Y are just repeated outside of the range of Xi.

//If decreasing, X and Xi must be monotonically decreasing.
//Otherwise, X and Xi must be monotonically increasing.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int interp1q_s (float *Yi, const float *X, const float *Y, const float *Xi, const size_t N, const size_t Ni, const char decreasing);
int interp1q_d (double *Yi, const double *X, const double *Y, const double *Xi, const size_t N, const size_t Ni, const char decreasing);
int interp1q_c (float *Yi, const float *X, const float *Y, const float *Xi, const size_t N, const size_t Ni, const char decreasing);
int interp1q_z (double *Yi, const double *X, const double *Y, const double *Xi, const size_t N, const size_t Ni, const char decreasing);


int interp1q_s (float *Yi, const float *X, const float *Y, const float *Xi, const size_t N, const size_t Ni, const char decreasing)
{
    if (N<1u) { fprintf(stderr,"error in interp1q_s: N must be positive\n"); return 1; }
    
    if (Ni==0u) {}
    else
    {
        float x1, y1;
        size_t n = 0u, ni = 0u;

        if (decreasing)
        {
            while (*Xi>=*X && ni<Ni) { *Yi = *Y; ++Xi; ++Yi; ++ni; }
            while (n<N && ni<Ni)
            {
                while (*X>=*Xi && n<N) { ++X; ++Y; ++n; }
                if (n<N)
                {
                    x1 = *(X-1); y1 = *(Y-1);
                    *Yi = y1 + (*Xi-x1)*(*Y-y1)/(*X-x1);
                    ++Xi; ++Yi; ++ni;
                }
            }
        }
        else
        {
            while (*Xi<=*X && ni<Ni) { *Yi = *Y; ++Xi; ++Yi; ++ni; }
            while (n<N && ni<Ni)
            {
                while (n<N && *X<=*Xi) { ++X; ++Y; ++n; }
                if (n<N)
                {
                    x1 = *(X-1); y1 = *(Y-1);
                    *Yi = y1 + (*Xi-x1)*(*Y-y1)/(*X-x1);
                    ++Xi; ++Yi; ++ni;
                }
            }
        }
        
        while (n<N) { ++X; ++Y; ++n; }
        y1 = *(Y-1);
        while (ni<Ni) { *Yi = y1; ++Yi; ++ni; }
    }

    return 0;
}


int interp1q_d (double *Yi, const double *X, const double *Y, const double *Xi, const size_t N, const size_t Ni, const char decreasing)
{
    if (N<1u) { fprintf(stderr,"error in interp1q_d: N must be positive\n"); return 1; }
    
    if (Ni==0u) {}
    else
    {
        double x1, y1;
        size_t n = 0u, ni = 0u;

        if (decreasing)
        {
            while (*Xi>=*X && ni<Ni) { *Yi = *Y; ++Xi; ++Yi; ++ni; }
            while (n<N && ni<Ni)
            {
                while (*X>=*Xi && n<N) { ++X; ++Y; ++n; }
                if (n<N)
                {
                    x1 = *(X-1); y1 = *(Y-1);
                    *Yi = y1 + (*Xi-x1)*(*Y-y1)/(*X-x1);
                    ++Xi; ++Yi; ++ni;
                }
            }
        }
        else
        {
            while (*Xi<=*X && ni<Ni) { *Yi = *Y; ++Xi; ++Yi; ++ni; }
            while (n<N && ni<Ni)
            {
                while (n<N && *X<=*Xi) { ++X; ++Y; ++n; }
                if (n<N)
                {
                    x1 = *(X-1); y1 = *(Y-1);
                    *Yi = y1 + (*Xi-x1)*(*Y-y1)/(*X-x1);
                    ++Xi; ++Yi; ++ni;
                }
            }
        }
        
        while (n<N) { ++X; ++Y; ++n; }
        y1 = *(Y-1);
        while (ni<Ni) { *Yi = y1; ++Yi; ++ni; }
    }

    return 0;
}


int interp1q_c (float *Yi, const float *X, const float *Y, const float *Xi, const size_t N, const size_t Ni, const char decreasing)
{
    if (N<1u) { fprintf(stderr,"error in interp1q_c: N must be positive\n"); return 1; }
    
    if (Ni==0u) {}
    else
    {
        float x1, y1r, y1i;
        size_t n = 0u, ni = 0u;

        if (decreasing)
        {
            while (*Xi>=*X && ni<Ni) { *Yi++ = *Y; *Yi++ = *Y; ++Xi; ++ni; }
            while (n<N && ni<Ni)
            {
                while (*X>=*Xi && n<N) { ++X; Y+=2; ++n; }
                if (n<N)
                {
                    x1 = *(X-1); y1r = *(Y-2); y1i = *(Y-1);
                    *Yi++ = y1r + (*Xi-x1)*(*Y-y1r)/(*X-x1);
                    *Yi++ = y1i + (*Xi-x1)*(*(Y+1)-y1i)/(*X-x1);
                    ++Xi; ++ni;
                }
            }
        }
        else
        {
            while (*Xi<=*X && ni<Ni) { *Yi++ = *Y; *Yi++ = *(Y+1); ++Xi; ++ni; }
            while (n<N && ni<Ni)
            {
                while (n<N && *X<=*Xi) { ++X; Y+=2; ++n; }
                if (n<N)
                {
                    x1 = *(X-1); y1r = *(Y-2); y1i = *(Y-1);
                    *Yi++ = y1r + (*Xi-x1)*(*Y-y1r)/(*X-x1);
                    *Yi++ = y1i + (*Xi-x1)*(*(Y+1)-y1i)/(*X-x1);
                    ++Xi; ++ni;
                }
            }
        }
        
        while (n<N) { ++X; Y+=2; ++n; }
        y1r = *(Y-2); y1i = *(Y-1);
        while (ni<Ni) { *Yi++ = y1r; *Yi++ = y1i; ++ni; }
    }

    return 0;
}


int interp1q_z (double *Yi, const double *X, const double *Y, const double *Xi, const size_t N, const size_t Ni, const char decreasing)
{
    if (N<1u) { fprintf(stderr,"error in interp1q_z: N must be positive\n"); return 1; }
    
    if (Ni==0u) {}
    else
    {
        double x1, y1r, y1i;
        size_t n = 0u, ni = 0u;

        if (decreasing)
        {
            while (*Xi>=*X && ni<Ni) { *Yi++ = *Y++; *Yi++ = *Y++; ++Xi; ++ni; }
            while (n<N && ni<Ni)
            {
                while (*X>=*Xi && n<N) { ++X; Y+=2; ++n; }
                if (n<N)
                {
                    x1 = *(X-1); y1r = *(Y-2); y1i = *(Y-1);
                    *Yi++ = y1r + (*Xi-x1)*(*Y-y1r)/(*X-x1);
                    *Yi++ = y1i + (*Xi-x1)*(*Y-y1i)/(*X-x1);
                    ++Xi; ++ni;
                }
            }
        }
        else
        {
            while (*Xi<=*X && ni<Ni) { *Yi++ = *Y++; *Yi++ = *Y++; ++Xi; ++ni; }
            while (n<N && ni<Ni)
            {
                while (n<N && *X<=*Xi) { ++X; Y+=2; ++n; }
                if (n<N)
                {
                    x1 = *(X-1); y1r = *(Y-2); y1i = *(Y-1);
                    *Yi++ = y1r + (*Xi-x1)*(*Y-y1r)/(*X-x1);
                    *Yi++ = y1i + (*Xi-x1)*(*Y-y1i)/(*X-x1);
                    ++Xi; ++ni;
                }
            }
        }
        
        while (n<N) { ++X; Y+=2; ++n; }
        y1r = *(Y-2); y1i = *(Y-1);
        while (ni<Ni) { *Yi++ = y1r; *Yi++ = y1i; ++ni; }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
