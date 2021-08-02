//Includes
#include "ifft.rad2.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 1u, O = 1u;
size_t dim, nfft;
int sc, rl;

//Description
string descr;
descr += "1D IFFT (inverse fast Fourier transform) of each vector (1D signal) in X,\n";
descr += "using a variant of the radix-2 algorithm (only valid for power-of-2 nfft).\n";
descr += "\n";
descr += "This is meant for inverting fft (the 1D FFT in this namespace).\n";
descr += "Thus, X must be complex-valued and have appropriate length.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension along which to transform.\n";
descr += "Use -d0 to operate along cols, -d1 to operate along rows, etc.\n";
descr += "The default is 0 (along cols), unless X is a row vector.\n";
descr += "\n";
descr += "The transform nfft is automatically set to the length of X along dim.\n";
descr += "Again, this is mean for inverting fft (the 1D FFT in this namespace).\n";
descr += "Instead, use or don't use the -r option to control output length.\n";
descr += "\n";
descr += "The output (Y) is complex-valued with length nfft along dim,\n";
descr += "unless using the -r (--real) option.\n";
descr += "\n";
descr += "Use -r (--real) if X represents n/2+1 nonnegative freqs.\n";
descr += "In this case (and only this case), the output Y is real-valued.\n";
descr += "Note that this assumes that the FFT nfft was even.\n";
descr += "To allow for odd-length FFT, use the -n (--nfft) option.\n";
descr += "The output Y always has length nfft (but may be real or complex).\n";
descr += "\n";
descr += "Include -s (--scale) to scale, to invert with scaled fft.\n";
descr += "Otherwise, Y will be scaled by 1/nfft, to invert with unscaled fft.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ ifft.rad2 X -o Y \n";
descr += "$ ifft.rad2 -d1 -r X > Y \n";
descr += "$ cat X | ifft.rad2 -n9 -r > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to transform [default=0]");
struct arg_int    *a_n = arg_intn("n","nfft","<uint>",0,1,"transform length [default=L]");
struct arg_lit   *a_rl = arg_litn("r","real",0,1,"return real Y (i.e., X is only nonnegative freqs)");
struct arg_lit   *a_sc = arg_litn("s","scale",0,1,"include to scale by sqrt(2) [default is 1/nfft]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = i1.isrowvec() ? 1u : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Get rl
rl = (a_rl->count>0);

//Get sc
sc = (a_sc->count>0);

//Get nfft
if (a_n->count==0)
{
    nfft = (dim==0u) ? i1.R : (dim==1u) ? i1.C : (dim==2u) ? i1.S : i1.H;
    if (rl) { nfft = 2u*(nfft-1u); }
}
else if (a_n->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be positive" << endl; return 1; }
else { nfft = size_t(a_n->ival[0]); }

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
if (i1.isreal()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) must be complex" << endl; return 1; }
if (dim==0u && nfft<i1.R) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be >= nrows X for dim=0" << endl; return 1; }
if (dim==1u && nfft<i1.C) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be >= ncols X for dim=1" << endl; return 1; }
if (dim==2u && nfft<i1.S) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be >= nslices X for dim=2" << endl; return 1; }
if (dim==3u && nfft<i1.H) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be >= nhyperslices X for dim=3" << endl; return 1; }

//Set output header info
o1.F = i1.F;
o1.T = (rl) ? i1.T-100u : i1.T;
o1.R = (dim==0u) ? nfft : i1.R;
o1.C = (dim==1u) ? nfft : i1.C;
o1.S = (dim==2u) ? nfft : i1.S;
o1.H = (dim==3u) ? nfft : i1.H;

//Other prep

//Process
if (o1.T==1u)
{
    float *X, *Y;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    if (codee::ifft_rad2_s(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,nfft,sc))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}
else if (o1.T==101u)
{
    float *X, *Y;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { Y = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    if (codee::ifft_rad2_c(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,nfft,sc)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish
