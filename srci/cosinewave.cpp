//Includes
#include <cfloat>
#include "cosinewave.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 0u, O = 1u;
size_t N, dim;
double amp, frq, phs;

//Description
string descr;
descr += "Generates 1-D cosinewave signal.\n";
descr += "This is same as sinewave, but uses different phase convention:\n";
descr += "For cosinewave the signal is at its maximum at phase 0,\n";
descr += "whereas for sinewave the signal is at 0 at phase 0.\n";
descr += "\n";
descr += "Use -l (--length) to give the output vector length in sample points.\n";
descr += "\n";
descr += "Use -d (--dim) to give the nonsingleton dim of the output vec.\n";
descr += "If d=0, then Y is a column vector [default].\n";
descr += "If d=1, then Y is a row vector.\n";
descr += "(d=2 and d=3 are also possible, but rarely used.)\n";
descr += "\n";
descr += "Use -a (--amplitude) to give the amplitude (max) of the output.\n";
descr += "Note that this is half of the peak-to-peak amplitude.\n";
descr += "\n";
descr += "Use -f (--freq) to give the frequency in units of cycles/sample.\n";
descr += "This is the frequency in Hz divided by the sample rate in Hz.\n";
descr += "Thus, f should be 0<f<=0.5, where 0 is DC and 0.5 is Nyquist.\n";
descr += "\n";
descr += "Use -p (--phase) to give the phase in radians.\n";
descr += "This is usually in [0 2*pi), but the value will be taken modulo 2*pi.\n";
descr += "\n";
descr += "Since this is a generating function (no inputs), the output data type\n";
descr += "and file format can be specified by -t and -f, respectively. \n";
descr += "\n";
descr += "Examples:\n";
descr += "$ cosinewave -n32 -a2.5 -f0.2 -o Y \n";
descr += "$ cosinewave -n32 -f0.01 -p1.57079632679 > Y \n";
descr += "$ cosinewave -n32 -f0.01 -d1 -t1 -f101 > Y \n";

//Argtable
struct arg_dbl  *a_amp = arg_dbln("a","amp","<dbl>",0,1,"amplitude [default=1.0]");
struct arg_dbl  *a_frq = arg_dbln("f","freq","<dbl>",0,1,"frequency in cycles/sample [default=0.5]");
struct arg_dbl  *a_phs = arg_dbln("p","phase","<dbl>",0,1,"phase in radians [default=0.0]");
struct arg_int    *a_n = arg_intn("n","N","<uint>",0,1,"num samples in output [default=1]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"nonsingleton dimension [default=0 -> col vec]");
struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=2 -> double]");
struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=102 -> colmajor]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get o1.F
if (a_ofmt->count==0) { o1.F = 147u; }
else if (a_ofmt->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be nonnegative" << endl; return 1; }
else if (a_ofmt->ival[0]>255) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be < 256" << endl; return 1; }
else { o1.F = size_t(a_ofmt->ival[0]); }

//Get o1.T
if (a_otyp->count==0) { o1.T = 1u; }
else if (a_otyp->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "data type must be positive int" << endl; return 1; }
else { o1.T = size_t(a_otyp->ival[0]); }
if ((o1.T==oktypes).sum()==0)
{
    cerr << progstr+": " << __LINE__ << errstr << "output data type must be in " << "{";
    for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1]) ? "}" : ","); }
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

//Get amp
amp = (a_amp->count>0) ? a_amp->dval[0] : 1.0;
if (amp<0.0) { cerr << progstr+": " << __LINE__ << errstr << "amplitude must be nonnegative" << endl; return 1; }
if (o1.T==1 && amp>=double(FLT_MAX)) { cerr << progstr+": " << __LINE__ << errstr << "amplitude must be < " << double(FLT_MAX) << endl; return 1; }

//Get frq
frq = (a_frq->count>0) ? a_frq->dval[0] : 0.5;
if (frq<double(FLT_EPSILON)) { cerr << progstr+": " << __LINE__ << errstr << "frequency must be positive" << endl; return 1; }
if (frq>0.5) { cerr << progstr+": " << __LINE__ << errstr << "frequency must be <= 0.5" << endl; return 1; }

//Get phs
phs = (a_phs->count>0) ? a_phs->dval[0] : 0.0;

//Checks

//Set output header info
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
    if (codee::cosinewave_s(Y,o1.N(),(float)amp,(float)frq,(float)phs))
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
    if (codee::cosinewave_c(Y,N,(float)amp,(float)frq,(float)phs))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] Y;
}

//Finish
