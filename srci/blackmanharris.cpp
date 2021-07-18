//Includes
#include "blackmanharris.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 0u, O = 1u;
size_t L, dim, norm;

//Description
string descr;
descr += "Gets 1D Blackman-Harris window.\n";
descr += "\n";
descr += "Use -l (--winlength) to give the window length in number of taps (sample points).\n";
descr += "\n";
descr += "Use -d (--dim) to give the nonsingleton dim of the output vec.\n";
descr += "If d=0, then Y is a column vector [default].\n";
descr += "If d=1, then Y is a row vector.\n";
descr += "(d=2 and d=3 are also possible, but rarely used.)\n";
descr += "\n";
descr += "Use -n (--norm) to normalize (divide) the output by a norm:\n";
descr += "Use -n0 to do no normalization [default].\n";
descr += "Use -n1 to normalize by the L1-norm (sum abs values).\n";
descr += "Use -n2 to normalize by the L2-norm (root-sum-of-squares).\n";
descr += "Use -n3 to normalize by the Inf-norm (maximum).\n";
descr += "\n";
descr += "Since this is a generating function (no inputs), the output data type\n";
descr += "and file format can be specified by -t and -f, respectively. \n";
descr += "\n";
descr += "Examples:\n";
descr += "$ blackmanharris -l255 -o Y \n";
descr += "$ blackmanharris -l255 > Y \n";
descr += "$ blackmanharris -l127 -d1 -t1 -f101 > Y \n";

//Argtable
struct arg_int   *a_wl = arg_intn("l","winlength","<uint>",0,1,"window length [default=7]");
struct arg_int  *a_nrm = arg_intn("n","norm","<uint>",0,1,"normalize output by norm [default=0]");
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
else if (a_otyp->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "data type must be positive int" << endl; return 1; }
else { o1.T = size_t(a_otyp->ival[0]); }
if ((o1.T==oktypes).sum()==0)
{
    cerr << progstr+": " << __LINE__ << errstr << "output data type must be in " << "{";
    for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1u]) ? "}" : ","); }
    cerr << endl; return 1;
}

//Get L
if (a_wl->count==0) { L = 7u; }
else if (a_wl->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "winlength must be a positive int" << endl; return 1; }
else { L = size_t(a_wl->ival[0]); }

//Get dim
if (a_d->count==0) { dim = 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Get norm
if (a_nrm->count==0) { norm = 0u; }
else if (a_nrm->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "norm must be nonnegative" << endl; return 1; }
else { norm = size_t(a_nrm->ival[0]); }
if (norm>3u) { cerr << progstr+": " << __LINE__ << errstr << "norm must be in {0,1,2,3}" << endl; return 1; }

//Checks

//Set output header info
o1.R = (dim==0u) ? L : 1u;
o1.C = (dim==1u) ? L : 1u;
o1.S = (dim==2u) ? L : 1u;
o1.H = (dim==3u) ? L : 1u;

//Other prep

//Process
if (o1.T==1u)
{
    float *Y;
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    if (codee::blackmanharris_s(Y,L,norm))
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
    if (codee::blackmanharris_c(Y,L,norm))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] Y;
}

//Finish
