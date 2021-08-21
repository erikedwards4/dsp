//Includes
#include "xcorr1d.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 2u, O = 1u;
size_t dim, L1, W, str, dil, pad;

//Description
string descr;
descr += "1D cross-correlation of each vector in X1 with X2.\n";
descr += "This version gives control over stride, dilation, etc.\n";
descr += "For a simpler interface with less flexibility, see xcorr.\n";
descr += "\n";
descr += "X2 is a vector that is NOT in flipped order.\n";
descr += "Note that some \"convolution\" functions actually do cross-corr.\n";
descr += "To use X2 in flipped order (convolution), use conv1d.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension (axis) along which to filter.\n";
descr += "Use -d0 to operate along cols, -d1 to operate along rows, etc.\n";
descr += "The default is 0 (along cols), unless X1 is a row vector.\n";
descr += "\n";
descr += "Use -s (--stride) to give the stride (step-size) in samples [default=1].\n";
descr += "\n";
descr += "Use -d (--dilation) to give the dilation factor [default=1].\n";
descr += "\n";
descr += "Use -p (--padding) to give the padding in samples [default=0]\n";
descr += "\n";
descr += "Y has the same size as X1, but the vecs along dim have length\n";
descr += "W = L_out = floor[1 + (N1-N2)/stride].\n";
descr += "          = floor[1 + (L1 + 2*pad - dil*(L2-1) - 1)/stride].\n";
descr += "where\n";
descr += "L1 = L_in is the length of vecs along dim in X1.\n";
descr += "L2 is the length of X2.\n";
descr += "N1 is the full length of vecs in X1 including padding.\n";
descr += "N2 is the full length of X2 including dilation.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ xcorr1d X1 X2 -o Y \n";
descr += "$ xcorr1d -d1 -p2 X1 X2 > Y \n";
descr += "$ cat X2 | xcorr1d -d1 -p5 -s3 -i2 X2 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X1,X2)");
struct arg_int  *a_str = arg_intn("s","step","<uint>",0,1,"step size in samps [default=1]");
struct arg_int  *a_dil = arg_intn("i","dilation","<uint>",0,1,"dilation factor [default=1]");
struct arg_int  *a_pad = arg_intn("p","padding","<uint>",0,1,"padding [default=0]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to filter [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = i1.isvec() ? i1.nonsingleton1() : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }
L1 = (dim==0u) ? i1.R : (dim==1u) ? i1.C : (dim==2u) ? i1.S : i1.H;

//Get padding
if (a_pad->count==0) { pad = 0u; }
else if (a_pad->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "pad must be nonnegative" << endl; return 1; }
else { pad = size_t(a_pad->ival[0]); }

//Get stride
if (a_str->count==0) { str = 1u; }
else if (a_str->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "stride must be positive" << endl; return 1; }
else { str = size_t(a_str->ival[0]); }

//Get dil
if (a_dil->count==0) { dil = 1u; }
else if (a_dil->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "dilation must be positive" << endl; return 1; }
else { dil = size_t(a_dil->ival[0]); }

//Checks
if (i1.T!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) found to be empty" << endl; return 1; }
if (!i2.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) must be a vector" << endl; return 1; }
if (L1+2u*pad<=dil*(i2.N()-1u)) { cerr << progstr+": " << __LINE__ << errstr << "L1+2*pad must be > dil*(L2-1)" << endl; return 1; }

//Set output header info
W = 1u + (L1 + 2u*pad - dil*(i2.N()-1u) - 1u) / str;
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
    if (codee::xcorr1d_s(Y,X1,X2,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),i2.N(),pad,str,dil,dim))
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
    if (codee::xcorr1d_c(Y,X1,X2,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),i2.N(),pad,str,dil,dim))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; } 
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X1; delete[] X2; delete[] Y;
}

//Finish
