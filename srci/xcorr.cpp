//Includes
#include "xcorr.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 2u, O = 1u;
size_t dim, L1, W;
string shape;

//Description
string descr;
descr += "1D cross-correlation of each vector in X1 by X2.\n";
descr += "This is not identical to the xcorr function of Octave.\n";
descr += "Rather, it is identical to conv, but without flipping X2.\n";
descr += "\n";
descr += "X2 is a vector that is NOT in flipped order.\n";
descr += "Note that some \"convolution\" functions actually do cross-corr.\n";
descr += "To use X2 in flipped order (convolution), use conv.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension (axis) along which to operate.\n";
descr += "Use -d0 to operate along cols, -d1 to operate along rows, etc.\n";
descr += "The default is 0 (along cols), unless X1 is a row vector.\n";
descr += "\n";
descr += "Use -s (--shape) to give the shape as 'full', 'same' or 'valid' [default='full'].\n";
descr += "For 'same', Y has length L1 along dim (same as X1).\n";
descr += "For 'full', Y has length L1+2*L2-2 along dim.\n";
descr += "For 'valid', Y has length L1-L2+1 along dim.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ xcorr X1 X2 -o Y \n";
descr += "$ xcorr -d1 -s'same' X1 X2 > Y \n";
descr += "$ cat X2 | xcorr -d1 -s'valid' X1 - > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X1,X2)");
struct arg_str   *a_sh = arg_strn("s","shape","<str>",0,1,"shape [default='full']");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to filter [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = i1.isvec() ? i1.nonsingleton1() : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Get shape
if (a_sh->count==0) { shape = "full"; }
else
{
	try { shape = string(a_sh->sval[0]); }
	catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem getting string for shape" << endl; return 1; }
}
for (string::size_type c=0u; c<shape.size(); ++c) { shape[c] = char(tolower(shape[c])); }
if (shape!="full" && shape!="same" && shape!="valid") { cerr << progstr+": " << __LINE__ << errstr << "shape string must be 'full', 'same' or 'valid'" << endl; return 1; }

//Checks
if (i1.T!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) found to be empty" << endl; return 1; }
if (!i2.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) must be a vector" << endl; return 1; }

//Set output header info
L1 = (dim==0u) ? i1.R : (dim==1u) ? i1.C : (dim==2u) ? i1.S : i1.H;
if (shape=="full") { W = L1 + i2.N() - 1u; }
else if (shape=="same") { W = L1; }
else { W = (L1>=i2.N()) ? L1-i2.N()+1u : 0u; }
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0u) ? W : i1.R;
o1.C = (dim==1u) ? W : i1.C;
o1.S = (dim==2u) ? W : i1.S;
o1.H = (dim==3u) ? W : i1.H;

//Other prep

//Process
if (i1.T==1u)
{
    float *X1, *X2, *Y;
    try { X1 = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
    try { X2 = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
    if (codee::xcorr_s(Y,X1,X2,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),i2.N(),shape.c_str(),dim))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; } 
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X1; delete[] X2; delete[] Y;
}
else if (i1.T==101u)
{
    float *X1, *X2, *Y;
    try { X1 = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
    try { X2 = new float[2u*i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
    try { Y = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
    if (codee::xcorr_c(Y,X1,X2,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),i2.N(),shape.c_str(),dim))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; } 
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X1; delete[] X2; delete[] Y;
}

//Finish
