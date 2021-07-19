//Includes
#include "apply_win.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 2u, O = 1u;
size_t dim;

//Description
string descr;
descr += "Elementwise function for 2 inputs with broadcasting.\n";
descr += "Elementwise multiplication: Y = X1 .* X2,\n";
descr += "where X1 is a tensor (up to 4-D) and X2 is a vector (1-D).\n";
descr += "\n";
descr += "X2 is a real-valued window of length L,\n";
descr += "as made by one of the window functions.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dim along which to apply X2.\n";
descr += "X1 must have length L along dim.\n";
descr += "\n";
descr += "Output (Y) has the same size as X1.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ apply_win X1 X2 -o Y \n";
descr += "$ apply_win X1 X2 > Y \n";
descr += "$ apply_win X1 <(hamming -l401) > Y \n";
descr += "$ hamming -l401 | apply_win X1 - > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X1,X2)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to apply X2");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Checks
if (i1.T!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) found to be empty" << endl; return 1; }
if (!i2.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) must be a vector" << endl; return 1; }
if (dim==0u && i1.R!=i2.N()) { cerr << progstr+": " << __LINE__ << errstr << "winlength must equal R (nrows X1) for dim=0" << endl; return 1; }
if (dim==1u && i1.C!=i2.N()) { cerr << progstr+": " << __LINE__ << errstr << "winlength must equal C (ncols X1) for dim=1" << endl; return 1; }
if (dim==2u && i1.S!=i2.N()) { cerr << progstr+": " << __LINE__ << errstr << "winlength must equal S (nslices X1) for dim=2" << endl; return 1; }
if (dim==3u && i1.H!=i2.N()) { cerr << progstr+": " << __LINE__ << errstr << "winlength must equal H (nhyperslices X1) for dim=3" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = i1.R; o1.C = i1.C; o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1u)
{
    float *X1, *X2;
    try { X1 = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
    try { X2 = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
    if (codee::apply_win_inplace_s(X1,X2,i1.R,i1.C,i1.S,i1.H,i2.N(),dim,o1.iscolmajor()))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(X1),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X1; delete[] X2;
}
else if (i1.T==101u)
{
    float *X1, *X2;
    try { X1 = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
    try { X2 = new float[2u*i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
    if (codee::apply_win_inplace_c(X1,X2,i1.R,i1.C,i1.S,i1.H,i2.N(),dim,o1.iscolmajor()))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(X1),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X1; delete[] X2;
}

//Finish
