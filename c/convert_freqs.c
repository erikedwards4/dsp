//This converts between several basic frequency scales:
//Hz (linear), sqrt, cbrt, octave (log2).

#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int convert_freqs_s (float *frqs, const size_t F, const char in_scale[], const char out_scale[]);
int convert_freqs_d (double *frqs, const size_t F, const char in_scale[], const char out_scale[]);


int convert_freqs_s (float *frqs, const size_t F, const char in_scale[], const char out_scale[])
{
	if (strlen(in_scale)<2u) { fprintf(stderr,"error in convert_freqs_s: input freq scale must be string with length > 1\n"); return 1; }
	if (strlen(out_scale)<2u) { fprintf(stderr,"error in convert_freqs_s: output freq scale must be string with length > 1\n"); return 1; }

	//Convert values in frqs
	if (F==0u) {}
	else if (strncmp(in_scale,"hz",2u)==0) //Hz
	{
		if (strncmp(out_scale,"hz",2u)==0)
		{
		}
        else if (strncmp(out_scale,"sqrt",4u)==0)
		{
			for (size_t f=0; f<F; ++f) { frqs[f] = sqrtf(frqs[f]); }
		}
        else if (strncmp(out_scale,"cbrt",4u)==0)
		{
			for (size_t f=0; f<F; ++f) { frqs[f] = cbrtf(frqs[f]); }
		}
		else if (strncmp(out_scale,"octave",6u)==0)
		{
			for (size_t f=0; f<F; ++f) { frqs[f] = log2f(frqs[f]); }
		}
		else
		{
			fprintf(stderr,"error in convert_freqs_s: output frequency scale not recognized\n"); return 1;
		}
	}
    else if (strncmp(in_scale,"sqrt",2u)==0) //sqrt
	{
		if (strncmp(out_scale,"hz",2u)==0)
		{
            for (size_t f=0; f<F; ++f) { frqs[f] *= frqs[f]; }
		}
        else if (strncmp(out_scale,"sqrt",4u)==0)
		{
		}
        else if (strncmp(out_scale,"cbrt",4u)==0)
		{
			for (size_t f=0; f<F; ++f) { frqs[f] = cbrtf(frqs[f]*frqs[f]); }
		}
		else if (strncmp(out_scale,"octave",6u)==0)
		{
			for (size_t f=0; f<F; ++f) { frqs[f] = log2f(frqs[f]*frqs[f]); }
		}
		else
		{
			fprintf(stderr,"error in convert_freqs_s: output frequency scale not recognized\n"); return 1;
		}
	}
    else if (strncmp(in_scale,"cbrt",2u)==0) //cbrt
	{
		if (strncmp(out_scale,"hz",2u)==0)
		{
            for (size_t f=0; f<F; ++f) { frqs[f] *= frqs[f]*frqs[f]; }
		}
        else if (strncmp(out_scale,"sqrt",4u)==0)
		{
            for (size_t f=0; f<F; ++f) { frqs[f] = sqrtf(frqs[f]*frqs[f]*frqs[f]); }
		}
        else if (strncmp(out_scale,"cbrt",4u)==0)
		{
		}
		else if (strncmp(out_scale,"octave",6u)==0)
		{
			for (size_t f=0; f<F; ++f) { frqs[f] = log2f(frqs[f]*frqs[f]*frqs[f]); }
		}
		else
		{
			fprintf(stderr,"error in convert_freqs_s: output frequency scale not recognized\n"); return 1;
		}
	}
	else if (strncmp(in_scale,"octave",6u)==0) //Octaves
	{
		if (strncmp(out_scale,"hz",2u)==0)
		{
			for (size_t f=0; f<F; ++f) { frqs[f] = powf(2.0f,frqs[f]); }
		}
		else if (strncmp(out_scale,"sqrt",4u)==0)
		{
			for (size_t f=0; f<F; ++f) { frqs[f] = powf(2.0f,0.5f*frqs[f]); }
		}
        else if (strncmp(out_scale,"cbrt",4u)==0)
		{
			for (size_t f=0; f<F; ++f) { frqs[f] = cbrtf(powf(2.0f,frqs[f])); }
		}
		else if (strncmp(out_scale,"octave",6u)==0)
		{
		}
		else
		{
			fprintf(stderr,"error in convert_freqs_s: output frequency scale not recognized\n"); return 1;
		}
	}
	else
	{
		fprintf(stderr,"error in convert_freqs_s: input frequency scale not recognized\n"); return 1;
	}
	
	return 0;
}


int convert_freqs_d (double *frqs, const size_t F, const char in_scale[], const char out_scale[])
{
	if (strlen(in_scale)<2u) { fprintf(stderr,"error in convert_freqs_d: input freq scale must be string with length > 1\n"); return 1; }
	if (strlen(out_scale)<2u) { fprintf(stderr,"error in convert_freqs_d: output freq scale must be string with length > 1\n"); return 1; }

	//Convert values in frqs
	if (F==0u) {}
	else if (strncmp(in_scale,"hz",2u)==0) //Hz
	{
		if (strncmp(out_scale,"hz",2u)==0)
		{
		}
        else if (strncmp(out_scale,"sqrt",4u)==0)
		{
			for (size_t f=0; f<F; ++f) { frqs[f] = sqrt(frqs[f]); }
		}
        else if (strncmp(out_scale,"cbrt",4u)==0)
		{
			for (size_t f=0; f<F; ++f) { frqs[f] = cbrt(frqs[f]); }
		}
		else if (strncmp(out_scale,"octave",6u)==0)
		{
			for (size_t f=0; f<F; ++f) { frqs[f] = log2(frqs[f]); }
		}
		else
		{
			fprintf(stderr,"error in convert_freqs_d: output frequency scale not recognized\n"); return 1;
		}
	}
    else if (strncmp(in_scale,"sqrt",2u)==0) //sqrt
	{
		if (strncmp(out_scale,"hz",2u)==0)
		{
            for (size_t f=0; f<F; ++f) { frqs[f] *= frqs[f]; }
		}
        else if (strncmp(out_scale,"sqrt",4u)==0)
		{
		}
        else if (strncmp(out_scale,"cbrt",4u)==0)
		{
			for (size_t f=0; f<F; ++f) { frqs[f] = cbrt(frqs[f]*frqs[f]); }
		}
		else if (strncmp(out_scale,"octave",6u)==0)
		{
			for (size_t f=0; f<F; ++f) { frqs[f] = log2(frqs[f]*frqs[f]); }
		}
		else
		{
			fprintf(stderr,"error in convert_freqs_d: output frequency scale not recognized\n"); return 1;
		}
	}
    else if (strncmp(in_scale,"cbrt",2u)==0) //cbrt
	{
		if (strncmp(out_scale,"hz",2u)==0)
		{
            for (size_t f=0; f<F; ++f) { frqs[f] *= frqs[f]*frqs[f]; }
		}
        else if (strncmp(out_scale,"sqrt",4u)==0)
		{
            for (size_t f=0; f<F; ++f) { frqs[f] = sqrt(frqs[f]*frqs[f]*frqs[f]); }
		}
        else if (strncmp(out_scale,"cbrt",4u)==0)
		{
		}
		else if (strncmp(out_scale,"octave",6u)==0)
		{
			for (size_t f=0; f<F; ++f) { frqs[f] = log2(frqs[f]*frqs[f]*frqs[f]); }
		}
		else
		{
			fprintf(stderr,"error in convert_freqs_d: output frequency scale not recognized\n"); return 1;
		}
	}
	else if (strncmp(in_scale,"octave",6u)==0) //Octaves
	{
		if (strncmp(out_scale,"hz",2u)==0)
		{
			for (size_t f=0; f<F; ++f) { frqs[f] = pow(2.0,frqs[f]); }
		}
		else if (strncmp(out_scale,"sqrt",4u)==0)
		{
			for (size_t f=0; f<F; ++f) { frqs[f] = pow(2.0,0.5*frqs[f]); }
		}
        else if (strncmp(out_scale,"cbrt",4u)==0)
		{
			for (size_t f=0; f<F; ++f) { frqs[f] = cbrt(pow(2.0,frqs[f])); }
		}
		else if (strncmp(out_scale,"octave",6u)==0)
		{
		}
		else
		{
			fprintf(stderr,"error in convert_freqs_d: output frequency scale not recognized\n"); return 1;
		}
	}
	else
	{
		fprintf(stderr,"error in convert_freqs_d: input frequency scale not recognized\n"); return 1;
	}
	
	return 0;
}


#ifdef __cplusplus
}
}
#endif
