//Includes
#include <cfloat>
#include "ac2cc.c"

//Declarations
const valarray<size_t> oktypes = {1u};
const size_t I = 1u, O = 1u;
size_t dim, Lx, K;
double preg;

//Description
string descr;
descr += "Gets cepstral coeffs (CCs) starting from the autocovariance (AC).\n";
descr += "Does Levinson-Durbin recursion of each vector in X,\n";
descr += "where each vector in X is one AC function.\n";
descr += "\n";
descr += "Use -k (--K) to specify the number of CCs to compute [default=Lx-1].\n";
descr += "X must have at least K+1 lags (i.e., Lx>K).\n";
descr += "\n";
descr += "A small regularizing power is added to the raw power before log compression.\n";
descr += "Use -p (--preg) to specify this constant [default=FLT_EPS].\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension along which to operate.\n";
descr += "Default is 0 (along cols), unless X is a row vector.\n";
descr += "\n";
descr += "If dim==0, then Y has size K x C x S x H\n";
descr += "If dim==1, then Y has size R x K x S x H.\n";
descr += "If dim==2, then Y has size R x C x K x H.\n";
descr += "If dim==3, then Y has size R x C x S x K.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ ac2cc -k13 X -o Y -o V \n";
descr += "$ ac2cc -d1 -k40 X -o Y -o V \n";
descr += "$ cat X | ac2cc -d1 -p1e-5 -k13 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to operate [default=0]");
struct arg_int    *a_k = arg_intn("k","K","<uint>",0,1,"num CCs to compute [default=Lx-1]");
struct arg_dbl *a_preg = arg_dbln("p","preg","<dbl>",0,1,"power regularization constant [default=FLT_EPS]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = i1.isrowvec() ? 1u : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Get K
Lx = (dim==0u) ? i1.R : (dim==1u) ? i1.C : (dim==2u) ? i1.S : i1.H;
if (a_k->count==0) { K = Lx - 1u; }
else if (a_d->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "K must be positve" << endl; return 1; }
else { K = size_t(a_k->ival[0]); }

//Get preg
preg = (a_preg->count>0) ? a_preg->dval[0] : double(FLT_EPSILON);
if (preg<0.0) { cerr << progstr+": " << __LINE__ << errstr << "preg must be nonnegative" << endl; return 1; }

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
if (Lx<2u) { cerr << progstr+": " << __LINE__ << errstr << "Lx (length of vecs in X) must be > 1" << endl; return 1; }
if (Lx<=K) { cerr << progstr+": " << __LINE__ << errstr << "K (num CCs) must be < Lx (length of vecs in X)" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0u) ? K : i1.R;
o1.C = (dim==1u) ? K : i1.C;
o1.S = (dim==2u) ? K : i1.S;
o1.H = (dim==3u) ? K : i1.H;

//Other prep

//Process
if (i1.T==1u)
{
    float *X, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::ac2cc_s(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,K,preg)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish
