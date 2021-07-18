//Includes
#include <float.h>
#include "get_stft_freqs.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u};
const size_t I = 0u, O = 1u;
size_t nfft, F, dim;
double sr;

//Description
string descr;
descr += "Gets the first F frequencies for an STFT of order nfft.\n";
descr += "\n";
descr += "Use -n (--nfft) to give the STFT transform length.\n";
descr += "\n";
descr += "Use -f (--F) to give the number of freqs in the output [default=nfft/2+1].\n";
descr += "Note that F must be <= nfft.\n";
descr += "If F>nfft/2+1, then the STFT negative freqs will be included.\n";
descr += "\n";
descr += "Use -r (--srate) to give the sample rate in units of Hz [default=1].\n";
descr += "\n";
descr += "Use -d (--dim) to give the nonsingleton dim of the output vec.\n";
descr += "If d=0, then Y is a column vector [default].\n";
descr += "If d=1, then Y is a row vector.\n";
descr += "(d=2 and d=3 are also possible, but rarely used.)\n";
descr += "\n";
descr += "Since this is a generating function (no inputs), the output data type\n";
descr += "and file format can be specified by -t and -f, respectively. \n";
descr += "\n";
descr += "Examples:\n";
descr += "$ get_stft_freqs -n255 -r8000 -o Y \n";
descr += "$ get_stft_freqs -n512 -r16000 > Y \n";
descr += "$ get_stft_freqs -n127 -r1000 -d1 -t1 -f101 > Y \n";

//Argtable
struct arg_int    *a_n = arg_intn("n","nfft","<uint>",0,1,"STFT transform length [default=0]");
struct arg_int    *a_f = arg_intn("f","F","<uint>",0,1,"number of freqs to output [default=nfft/2+1]");
struct arg_dbl   *a_sr = arg_dbln("r","srate","<dbl>",0,1,"sample rate in Hz [default=1]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"nonsingleton dimension [default=0 -> col vec]");
struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=2 -> double]");
struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=102 -> colmajor]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get o1.F
if (a_ofmt->count==0) { o1.F = 147u; }
else if (a_ofmt->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be nonnegative" << endl; return 1; }
else if (a_ofmt->ival[0]>255) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be < 256" << endl; return 1; }
else { o1.F = size_t(a_ofmt->ival[0]); }

//Get o1.T
if (a_otyp->count==0) { o1.T = 1u; }
else if (a_otyp->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "data type must be positive int" << endl; return 1; }
else { o1.T = size_t(a_otyp->ival[0]); }
if ((o1.T==oktypes).sum()==0)
{
    cerr << progstr+": " << __LINE__ << errstr << "output data type must be in " << "{";
    for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1]) ? "}" : ","); }
    cerr << endl; return 1;
}

//Get nfft
if (a_n->count==0) { nfft = 0u; }
else if (a_n->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be a positive int" << endl; return 1; }
else { nfft = size_t(a_n->ival[0]); }
if (nfft<1u) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be a positive int" << endl; return 1; }

//Get F
if (a_f->count==0) { F = nfft/2u + 1u; }
else if (a_f->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "F must be a positive int" << endl; return 1; }
else { F = size_t(a_f->ival[0]); }
if (F<1u) { cerr << progstr+": " << __LINE__ << errstr << "F must be a positive int" << endl; return 1; }
if (F>nfft) { cerr << progstr+": " << __LINE__ << errstr << "F must be <= nfft>" << endl; return 1; }

//Get sr
if (a_sr->count==0) { sr = 1.0; }
else if (a_sr->dval[0]<DBL_EPSILON) { cerr << progstr+": " << __LINE__ << errstr << "sample rate must be positive" << endl; return 1; }
else { sr = a_sr->dval[0]; }

//Get dim
if (a_d->count==0) { dim = 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Checks

//Set output header info
o1.R = (dim==0u) ? F : 1u;
o1.C = (dim==1u) ? F : 1u;
o1.S = (dim==2u) ? F : 1u;
o1.H = (dim==3u) ? F : 1u;

//Other prep

//Process
if (o1.T==1u)
{
    float *Y;
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    if (codee::get_stft_freqs_s(Y,nfft,F,(float)sr))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] Y;
}

//Finish