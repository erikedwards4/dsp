//Includes
#include "ar2psd.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2,101,102};
const size_t I = 3, O = 1;
size_t dim;

//Description
string descr;
descr += "Gets power spectral density (PSD) from autoregressive (AR) params along rows or cols of X.\n";
descr += "Also takes a 2nd input V, which is a vector of variances (prediction errors).\n";
descr += "Also takes a 3rd input W, which is a vector of F freqs (in radians).\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension along which to operate.\n";
descr += "Default is 0 (along cols), unless X is a row vector.\n";
descr += "\n";
descr += "If dim==0, then V must have length C, and Y has size F x C.\n";
descr += "If dim==1, then V must have length R, and Y has size R x F,\n";
descr += "where X has size R x C.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ ar2psd X V W -o Y \n";
descr += "$ ar2psd X V W > Y \n";
descr += "$ cat W | ar2psd X V > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X,V,W)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to operate [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = (i1.R==1u) ? 1 : 0; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim!=0 && dim!=1) { cerr << progstr+": " << __LINE__ << errstr << "dim must be 0 or 1" << endl; return 1; }

//Checks
if (i2.T!=i3.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs 2 and 3 must have the same data type" << endl; return 1; }
if (i1.isreal() && i1.T!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (i1.iscomplex() && (i1.T-100)!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have compatible data types" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (V) found to be empty" << endl; return 1; }
if (i3.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 3 (W) found to be empty" << endl; return 1; }
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) must be a matrix" << endl; return 1; }
if (!i2.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (V) must be a vector" << endl; return 1; }
if (!i3.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 3 (W) must be a vector" << endl; return 1; }
if (dim==0 && i1.C!=i2.N()) { cerr << progstr+": " << __LINE__ << errstr << "inputs 1 and 2 must be size compatible" << endl; return 1; }
if (dim==1 && i1.R!=i2.N()) { cerr << progstr+": " << __LINE__ << errstr << "inputs 1 and 2 must be size compatible" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = (i1.isreal()) ? i1.T : i1.T-100;
o1.R = (dim==0) ? i3.N() : i1.R;
o1.C = (dim==1) ? i3.N() : i1.C;
o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1)
{
    float *X, *V, *W, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { V = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (V)" << endl; return 1; }
    try { W = new float[i3.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 3 (W)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(V),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (V)" << endl; return 1; }
    try { ifs3.read(reinterpret_cast<char*>(W),i3.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 3 (W)" << endl; return 1; }
    if (codee::ar2psd_s(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),V,W,int(i3.N()),dim))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] V; delete[] W; delete[] Y;
}
else if (i1.T==101)
{
    float *X, *V, *W, *Y;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { V = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (V)" << endl; return 1; }
    try { W = new float[i3.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 3 (W)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(V),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (V)" << endl; return 1; }
    try { ifs3.read(reinterpret_cast<char*>(W),i3.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 3 (W)" << endl; return 1; }
    if (codee::ar2psd_c(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),V,W,int(i3.N()),dim))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] V; delete[] W; delete[] Y;
}

//Finish
