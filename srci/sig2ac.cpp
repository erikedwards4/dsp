//Includes
#include "sig2ac.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 1u, O = 1u;
size_t dim, L;
char u;

//Description
string descr;
descr += "Gets autocovariance of each column or row of X.\n";
descr += "The means are NOT subtracted before taking dot products.\n";
descr += "This does not normalize lag 0 of Y (so just like Octave xcorr).\n";
descr += "\n";
descr += "Use -l (--L) to give the number of lags at which to compute.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension along which to operate.\n";
descr += "Default is 0 (along cols), unless X is a row vector.\n";
descr += "\n";
descr += "If dim==0, then Y has size L x C.\n";
descr += "If dim==1, then Y has size R x L.\n";
descr += "\n";
descr += "Include -u (--unbiased) to use unbiased calculation [default is biased].\n";
descr += "This uses N-l instead of N in the denominator (it is actually just less biased).\n";
descr += "This takes longer, doesn't match FFT estimate, and has larger MSE (see Wikipedia).\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ sig2ac -l127 X -o Y \n";
descr += "$ sig2ac -d1 -l127 X > Y \n";
descr += "$ cat X | sig2ac -l127 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_l = arg_intn("l","L","<uint>",0,1,"number of lags to compute [default=1]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to operate [default=0]");
struct arg_lit    *a_u = arg_litn("u","unbiased",0,1,"use unbiased (N-l) denominator [default=biased]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get L
if (a_l->count==0) { L = 1; }
else if (a_l->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "L (nlags) must be positive" << endl; return 1; }
else { L = a_l->ival[0]; }

//Get dim
if (a_d->count==0) { dim = i1.isrowvec() ? 1u : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim!=0u && dim!=1u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be 0 or 1" << endl; return 1; }

//Get u
u = (a_u->count>0);

//Checks
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input must be 1D or 2D" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
if (dim==0u && i1.R<2u) { cerr << progstr+": " << __LINE__ << errstr << "cannot work along a singleton dimension" << endl; return 1; }
if (dim==1u && i1.C<2u) { cerr << progstr+": " << __LINE__ << errstr << "cannot work along a singleton dimension" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0u) ? uint32_t(L) : i1.R;
o1.C = (dim==1u) ? uint32_t(L) : i1.C;
o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1u)
{
    float *X, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::sig2ac_s(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,L,u)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}
else if (i1.T==101u)
{
    float *X, *Y;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { Y = new float[2*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::sig2ac_c(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,L,u)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish
