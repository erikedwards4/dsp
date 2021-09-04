//Includes
#include "smooth_diff.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u};
const size_t I = 0u, O = 1u;
size_t dim, n, N;

//Description
string descr;
descr += "Gets filter coefficients for the smooth differentiators of Holoborodko [2008].\n";
descr += "This is a generating function (no inputs except opts and params).\n";
descr += "These are MA filters; the resulting coeffs (B) can be used with fir.\n";
descr += "\n";
descr += "These are anti-symmetric filters that are designed to allow exactness of\n";
descr += "diff for polynomials of degree n, but smoothness for noise robustness.\n";
descr += "\n";
descr += "Use -n (--degree) to give the polynomial degree in {2,4} [default=2]\n";
descr += "\n";
descr += "Use -l (--N) to give the filter length in {5,7,9,11} [default=7]\n";
descr += "For n=4, N must be in {7,9,11}.\n";
descr += "\n";
descr += "Use -d (--dim) to give the nonsingleton dim of the output vec.\n";
descr += "If d=0, then Y is a column vector [default].\n";
descr += "If d=1, then Y is a row vector.\n";
descr += "(d=2 and d=3 are also possible, but rarely used.)\n";
descr += "\n";
descr += "Since this is a generating function (no inputs), the output data type\n";
descr += "and file format can be specified by -t and -f, respectively. \n";
descr += "\n";
descr += "Examples:\n";
descr += "$ smooth_diff -l9 -o B \n";
descr += "$ smooth_diff -n4 -l11 > B \n";
descr += "$ smooth_diff -l7 -d1 -t1 -f101 > B \n";

//Argtable
struct arg_int    *a_n = arg_intn("n","degree","<uint>",0,1,"polynomial degree [default=2]");
struct arg_int    *a_l = arg_intn("l","N","<uint>",0,1,"output (filter coeffs) length [default=7]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"nonsingleton dimension [default=0 -> col vec]");
struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=1]");
struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=147]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (B)");

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
    for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1u]) ? "}" : ","); }
    cerr << endl; return 1;
}

//Get dim
if (a_d->count==0) { dim = 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Get n
if (a_n->count==0) { n = 2u; }
else if (a_n->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "n must be nonnegative" << endl; return 1; }
else { n = size_t(a_n->ival[0]); }
if (n!=2u && n!=4u) { cerr << progstr+": " << __LINE__ << errstr << "n (polynomial degree) must be in {2,4}" << endl; return 1; }

//Get N
if (a_l->count==0) { N = 7u; }
else if (a_l->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "N must be nonnegative" << endl; return 1; }
else { N = size_t(a_l->ival[0]); }
if (N%2u==0u) { cerr << progstr+": " << __LINE__ << errstr << "N (filter length) must be odd" << endl; return 1; }
if (N<5u || N>11u) { cerr << progstr+": " << __LINE__ << errstr << "N (filter length) must be in {5,7,9,11}" << endl; return 1; }

//Checks
if (n==4u && N==5u) { cerr << progstr+": " << __LINE__ << errstr << "N (filter length) must be in {7,9,11} for n=4" << endl; return 1; }

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
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (B)" << endl; return 1; }
    if (codee::smooth_diff_s(Y,N,n))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (B)" << endl; return 1; }
    }
    delete[] Y;
}

//Finish
