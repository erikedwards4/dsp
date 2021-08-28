//Includes
#include "ac2psd.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 2u, O = 1u;
size_t dim;

//Description
string descr;
descr += "Gets power spectral density (PSD) from autocovariance (AC) for each vector in X.\n";
descr += "The 2nd input is W, which is a vector of F freqs (in radians).\n";
descr += "The PSD is obtained only at the F freqs.\n";
descr += "\n";
descr += "In a typical use case, X is a vector of length L,\n";
descr += "where P is the order of the polynomial,\n";
descr += "and Y is a vector of length F (same length as W).\n";
descr += "\n";
descr += "However, matrices and tensors are supported as follows:\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension along which to operate.\n";
descr += "Default is 0 (along cols), unless X is a row vector.\n";
descr += "\n";
descr += "If dim==0, then Y has size FxCxSxH.\n";
descr += "If dim==1, then Y has size RxFxSxH.\n";
descr += "If dim==2, then Y has size RxCxFxH.\n";
descr += "If dim==3, then Y has size RxCxSxF.\n";
descr += "\n";
descr += "Thus, this is a vec2vec operation,\n";
descr += "where the input vecs (ACs) have length L,\n";
descr += "and the output vecs (PSDs) have length F,\n";
descr += "and the vecs can extend along any dim.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ ac2psd X W -o Y \n";
descr += "$ ac2psd -d1 X W > Y \n";
descr += "$ cat W | ac2psd X > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X,W)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to operate [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = i1.isrowvec() ? 1u : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Checks
if (i2.iscomplex()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (W) must have real-valued data type" << endl; return 1; }
if (i1.isreal() && i1.T!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (i1.iscomplex() && (i1.T-100)!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have compatible data types" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (W) found to be empty" << endl; return 1; }

//Set output header info
o1.F = i1.F;
o1.T = (i1.isreal()) ? i1.T : i1.T-100u;
o1.R = (dim==0u) ? i2.N() : i1.R;
o1.C = (dim==1u) ? i2.N() : i1.C;
o1.S = (dim==2u) ? i2.N() : i1.S;
o1.H = (dim==3u) ? i2.N() : i1.H;

//Other prep

//Process
if (i1.T==1u)
{
    float *X, *W, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { W = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (W)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(W),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (W)" << endl; return 1; }
    if (codee::ac2psd_s(Y,X,W,i2.N(),i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] W; delete[] Y;
}
else if (i1.T==101u)
{
    float *X, *W, *Y;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { W = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (W)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(W),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (W)" << endl; return 1; }
    if (codee::ac2psd_c(Y,X,W,i2.N(),i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] W; delete[] Y;
}

//Finish
