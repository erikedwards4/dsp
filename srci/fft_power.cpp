//Includes
#include "fft_power.c"

//Declarations
const valarray<size_t> oktypes = {101u,102u};
const size_t I = 1u, O = 1u;

//Description
string descr;
descr += "Gets real-valued power for complex-valued FFT output.\n";
descr += "\n";
descr += "This is simply element-wise xr^2 + xi^2, \n";
descr += "where xr, xi are the real, imaginary parts of one element of X.\n";
descr += "\n";
descr += "The output Y has the same size as X, but is real-valued.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ fft_power X -o Y \n";
descr += "$ fft_powerX > Y \n";
descr += "$ cat X | fft_power > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Checks
if (i1.isreal()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) must be complex-valued" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T-100u;
o1.R = i1.R; o1.C = i1.C; o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==101u)
{
    float *X, *Y;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    if (codee::fft_power_c(Y,X,i1.N()))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish
