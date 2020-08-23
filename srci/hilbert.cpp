//Includes
#include "hilbert.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u};
const size_t I = 1u, O = 1u;
size_t dim, nfft;

//Description
string descr;
descr += "1D Hilbert transform of each vector (1D signal) in X.\n";
descr += "The output (Y) is the analytic signal, with X in the real part,\n";
descr += "and the actual Hilbert transform in the imaginary part.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension along which to transform.\n";
descr += "Use -d0 to operate along cols, -d1 to operate along rows, etc.\n";
descr += "The default is 0 (along cols), unless X is a row vector.\n";
descr += "\n";
descr += "Use -n (--nfft) to specify transform length [default=L].\n";
descr += "The default (L) is the length of X along dim.\n";
descr += "X is zero-padded as necessary to match nfft.\n";
descr += "\n";
descr += "The output (Y) is complex-valued with the same size as X.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ hilbert -n256 X -o Y \n";
descr += "$ hilbert -n256 -d1 X > Y \n";
descr += "$ cat X | hilbert -n256 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to transform [default=0]");
struct arg_int    *a_n = arg_intn("n","nfft","<uint>",0,1,"transform length [default=L]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = (i1.isrowvec()==1u) ? 1u : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1u,2u,3}" << endl; return 1; }

//Get nfft
if (a_n->count==0) { nfft = (dim==0u) ? i1.R : (dim==1u) ? i1.C : (dim==2u) ? i1.S : i1.H; }
else if (a_n->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be positive" << endl; return 1; }
else { nfft = size_t(a_n->ival[0]); }

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
if (dim==0u && nfft<i1.R) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be >= nrows X for dim=0" << endl; return 1; }
if (dim==1u && nfft<i1.C) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be >= ncols X for dim=1" << endl; return 1; }
if (dim==2u && nfft<i1.S) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be >= nslices X for dim=2" << endl; return 1; }
if (dim==3u && nfft<i1.H) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be >= nhyperslices X for dim=3" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T+100u;
o1.R = i1.R; o1.C = i1.C;
o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1u)
{
    float *X, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { Y = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    if (codee::hilbert_s(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,nfft))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish
