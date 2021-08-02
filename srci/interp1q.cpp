//Includes
#include "interp1q.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 3u, O = 1u;
int decreasing;

//Description
string descr;
descr += "Vecs2vec function for 3 inputs (X,Y,Xi) and 1 output (Yi).\n";
descr += "Does quick (linear) interpolation from coords (X,Y) to (Xi,Yi).\n";
descr += "\n";
descr += "Both X and Xi must be monotonically increasing.\n";
descr += "Use -d (--decreasing) if X and Xi are monotonically decreasing\n";
descr += "\n";
descr += "For complex inputs, this interpolates real and imag parts separately.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ interp1q X Y Xi -o Yi \n";
descr += "$ interp1q X Y Xi > Yi \n";
descr += "$ cat Xi | interp1q -d X Y - > Yi \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X,Y,Xi)");
struct arg_lit  *a_dec = arg_litn("d","decreasing",0,1,"include if X, Xi are monotonically decreasing");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Yi)");

//Get options

//Get decreasing
decreasing = (a_dec->count>0);

//Checks
if (i1.iscomplex()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) must be real-valued" << endl; return 1; }
if (i3.iscomplex()) { cerr << progstr+": " << __LINE__ << errstr << "input 3 (Xi) must be real-valued" << endl; return 1; }
if (i1.T!=i2.T && i1.T!=i2.T-100u) { cerr << progstr+": " << __LINE__ << errstr << "inputs 1 and 2 must have compatible data types" << endl; return 1; }
if (i1.T!=i3.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs 1 (X) and 3 (Xi) must have the same data type" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (Y) found to be empty" << endl; return 1; }
if (!i1.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) must be a vector" << endl; return 1; }
if (!i2.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (Y) must be a vector" << endl; return 1; }
if (!i3.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 3 (Xi) must be a vector" << endl; return 1; }
if (i1.N()!=i2.N()) { cerr << progstr+": " << __LINE__ << errstr << "inputs 1 and 2 must have the same length" << endl; return 1; }

//Set output header info
o1.F = i2.F; o1.T = i2.T;
o1.R = i3.R; o1.C = i3.C; o1.S = i3.S; o1.H = i3.H;

//Other prep

//Process
if (o1.T==1u)
{
    float *X, *Y, *Xi, *Yi;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { Y = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (Y)" << endl; return 1; }
    try { Xi = new float[i3.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 3 (Xi)" << endl; return 1; }
    try { Yi = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Yi)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(Y),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (Y)" << endl; return 1; }
    try { ifs3.read(reinterpret_cast<char*>(Xi),i3.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 3 (Xi)" << endl; return 1; }
    if (codee::interp1q_s(Yi,X,Y,Xi,i1.N(),o1.N(),decreasing))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Yi),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Yi)" << endl; return 1; }
    }
    delete[] X; delete[] Y; delete[] Xi; delete[] Yi;
}
else if (o1.T==101u)
{
    float *X, *Y, *Xi, *Yi;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { Y = new float[2u*i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (Y)" << endl; return 1; }
    try { Xi = new float[i3.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 3 (Xi)" << endl; return 1; }
    try { Yi = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Yi)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(Y),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (Y)" << endl; return 1; }
    try { ifs3.read(reinterpret_cast<char*>(Xi),i3.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 3 (Xi)" << endl; return 1; }
    if (codee::interp1q_c(Yi,X,Y,Xi,i1.N(),o1.N(),decreasing))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Yi),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Yi)" << endl; return 1; }
    }
    delete[] X; delete[] Y; delete[] Xi; delete[] Yi;
}

//Finish
