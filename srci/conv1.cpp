//Includes
#include "conv1.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 2u, O = 1u;
size_t dim, L1, str, dil, Ly;
int es0;

//Description
string descr;
descr += "1D cross-correlation of each vector in X1 with X2.\n";
descr += "This version gives control over stride, start-samp, and output length.\n";
descr += "For a simpler interface with less flexibility, see conv.\n";
descr += "For a PyTorch-like interface with pad and dilation, see conv1d.\n";
descr += "\n";
descr += "X2 is a vector of length L2 that is in flipped order (see also xcorr1).\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension (axis) along which to operate.\n";
descr += "Use -d0 to operate along cols, -d1 to operate along rows, etc.\n";
descr += "The default is 0 (along cols), unless X1 is a row vector.\n";
descr += "\n";
descr += "Use -s (--stride) to give the stride (step-size) in samples [default=1].\n";
descr += "\n";
descr += "Use -i (--dilation) to give the dilation factor [default=1].\n";
descr += "\n";
descr += "Use -e (--es0) to give the end-samp of the initial frame [default=0].\n";
descr += "Set this to L2-1 if want first frame fully within X1 (e.g., 'valid').\n";
descr += "Set this to L2/2 if want first center-samp at 0 (e.g., 'same').\n";
descr += "Set this to 0 if want first frame to just overlap X1 (e.g., 'full').\n";
descr += "Set this to a negative int if want some initial zeros in the output.\n";
descr += "\n";
descr += "Use -l (--Ly) to give the length of output vecs in Y [default=L1/str]\n";
descr += "\n";
descr += "Y has the same size as X1, but the vecs along dim have length Ly.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ conv1 X1 X2 -o Y \n";
descr += "$ conv1 -s2 X1 X2 > Y \n";
descr += "$ cat X2 | conv1 -s2 -e15 -l400 X1 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X1,X2)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to operate [default=0]");
struct arg_int  *a_es0 = arg_intn("e","es0","<int>",0,1,"end sample of initial frame [default=0]");
struct arg_int  *a_str = arg_intn("s","stride","<uint>",0,1,"stride (step size) in samps [default=1]");
struct arg_int  *a_dil = arg_intn("i","dilation","<uint>",0,1,"dilation factor [default=1]");
struct arg_int   *a_ly = arg_intn("l","Ly","<uint>",0,1,"length of vecs in output Y [default=L1/str]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = i1.isvec() ? i1.nonsingleton1() : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Get es0
if (a_es0->count==0) { es0 = 0; }
else { es0 = a_es0->ival[0]; }

//Get stride
if (a_str->count==0) { str = 1u; }
else if (a_str->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "stride must be positive" << endl; return 1; }
else { str = size_t(a_str->ival[0]); }

//Get dil
if (a_dil->count==0) { dil = 1u; }
else if (a_dil->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "dilation must be positive" << endl; return 1; }
else { dil = size_t(a_dil->ival[0]); }

//Get Ly
L1 = (dim==0u) ? i1.R : (dim==1u) ? i1.C : (dim==2u) ? i1.S : i1.H;
if (a_ly->count==0) { Ly = L1/str; }
else if (a_ly->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "Ly (length of vecs in Y) must be positive" << endl; return 1; }
else { Ly = size_t(a_ly->ival[0]); }

//Checks
if (i1.T!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) found to be empty" << endl; return 1; }
if (!i2.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) must be a vector" << endl; return 1; }
if (L1<i2.N()) { cerr << progstr+": " << __LINE__ << errstr << "L1 must be >= L2)" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0u) ? Ly : i1.R;
o1.C = (dim==1u) ? Ly : i1.C;
o1.S = (dim==2u) ? Ly : i1.S;
o1.H = (dim==3u) ? Ly : i1.H;

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
    if (codee::conv1_s(Y,X1,X2,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),i2.N(),es0,str,dil,Ly,dim))
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
    if (codee::conv1_c(Y,X1,X2,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),i2.N(),es0,str,dil,Ly,dim))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; } 
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X1; delete[] X2; delete[] Y;
}

//Finish
