//Includes
#include <cfloat>
#include "blue.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 0u, O = 1u;
size_t N, dim;
char zmn;
double std;

//Description
string descr;
descr += "Zero-mean, Gaussian blue noise (1-D).\n";
descr += "Makes vector of Gaussian white noise with specified stddev, and then\n";
descr += "takes FFT, multiplies by the f^1 characteristic, and takes IFFT.\n";
descr += "\n";
descr += "This uses modified code from PCG random, but does not require it to be installed.\n";
descr += "\n";
descr += "Use -n (--N) to give the number of samples in the output vec.\n";
descr += "\n";
descr += "Use -d (--dim) to give the nonsingleton dim of the output vec.\n";
descr += "If d=0, then Y is a column vector [default].\n";
descr += "If d=1, then Y is a row vector.\n";
descr += "(d=2 and d=3 are also possible, but rarely used.)\n";
descr += "\n";
descr += "Use -z (--zero_mean) to zero the mean before output.\n";
descr += "Although the expected value of the mean is 0, the generated mean is not.\n";
descr += "\n";
descr += "For complex output, real/imag parts are separately set using the same params.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ blue -n16 -o Y \n";
descr += "$ blue -d1 -n16 -z -t1 > Y \n";
descr += "$ blue -d1 -n16 -t102 > Y \n";

//Argtable
struct arg_dbl  *a_std = arg_dbln("s","std","<dbl>",0,1,"std dev parameter [default=1.0]");
struct arg_lit  *a_zmn = arg_litn("z","zero_mean",0,1,"include to zero mean before output");
struct arg_int    *a_n = arg_intn("n","N","<uint>",0,1,"num samples in output [default=1]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"nonsingleton dimension [default=0 -> col vec]");
struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=1]");
struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=147]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get o1.F
if (a_ofmt->count==0) { o1.F = 147u; }
else if (a_ofmt->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be nonnegative" << endl; return 1; }
else if (a_ofmt->ival[0]>255) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be < 256" << endl; return 1; }
else { o1.F = size_t(a_ofmt->ival[0]); }

//Get o1.T
if (a_otyp->count==0) { o1.T = 1u; }
else if (a_otyp->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "output data type must be positive" << endl; return 1; }
else { o1.T = size_t(a_otyp->ival[0]); }
if ((o1.T==oktypes).sum()==0)
{
    cerr << progstr+": " << __LINE__ << errstr << "output data type must be in " << "{";
    for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1u]) ? "}" : ","); }
    cerr << endl; return 1;
}

//Get N
if (a_n->count==0) { N = 1u; }
else if (a_n->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "N must be a positive int" << endl; return 1; }
else { N = size_t(a_n->ival[0]); }

//Get dim
if (a_d->count==0) { dim = 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Get std
std = (a_std->count>0) ? a_std->dval[0] : 1.0;
if (std<0.0) { cerr << progstr+": " << __LINE__ << errstr << "s (stddev) must be nonnegative " << endl; return 1; }
if (o1.T==1u && std>=double(FLT_MAX)) { cerr << progstr+": " << __LINE__ << errstr << "s (stddev) must be < " << double(FLT_MAX) << endl; return 1; }

//Get zmn
zmn = (a_zmn->count>0);

//Checks

//Set output header
o1.R = (dim==0u) ? N : 1u;
o1.C = (dim==1u) ? N : 1u;
o1.S = (dim==2u) ? N : 1u;
o1.H = (dim==3u) ? N : 1u;

//Other prep

//Process
if (o1.T==1u)
{
    float *Y;
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    if (codee::blue_s(Y,o1.N(),(float)std,zmn))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] Y;
}
else if (o1.T==101u)
{
    float *Y;
    try { Y = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    if (codee::blue_c(Y,o1.N(),(float)std,zmn))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] Y;
}

//Finish
