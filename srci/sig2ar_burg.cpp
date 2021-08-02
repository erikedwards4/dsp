//Includes
#include "sig2ar_burg.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u};
const size_t I = 1u, O = 2u;
size_t dim, P;
int m0;

//Description
string descr;
descr += "Does linear prediction for each row or col of X,\n";
descr += "using the Burg autoregressive (AR) method,\n";
descr += "which uses time-domain forward and backward predictions.\n";
descr += "\n";
descr += "Use -p (--P) to specify the number of LP coefficients.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension along which to operate.\n";
descr += "Default is 0 (along cols), unless X is a row vector.\n";
descr += "\n";
descr += "Include -m (--mean0) to zero the mean of each row or col of X [default=no].\n";
descr += "This is done before computing AR coeffs (works along same dim as above).\n";
descr += "\n";
descr += "If dim==0, then Y has size P x C.\n";
descr += "If dim==1, then Y has size R x P.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ sig2ar_burg -p3 X -o Y -o V \n";
descr += "$ sig2ar_burg -d1 -p5 X -o Y -o V \n";
descr += "$ cat X | sig2ar_burg -p7 -m > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_p = arg_intn("p","P","<uint>",0,1,"number of LP coeffs [default=1]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to operate [default=0]");
struct arg_lit   *a_m0 = arg_litn("m","mean0",0,1,"zero mean before computing [default=no]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output files (Y,V)");

//Get options

//Get dim
if (a_d->count==0) { dim = i1.isrowvec() ? 1u : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim!=0u && dim!=1u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be 0 or 1" << endl; return 1; }

//Get L
if (a_p->count==0) { P = 1; }
else if (a_p->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "P must be positive" << endl; return 1; }
else { P = a_p->ival[0]; }

//Get m0
m0 = (a_m0->count>0);

//Checks
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input must be 1D or 2D" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }

//Set output header info
o1.F = o2.F = i1.F;
o1.T = o2.T = i1.T;
o1.R = (dim==0u) ? uint32_t(P) : i1.R;
o1.C = (dim==1u) ? uint32_t(P) : i1.C;
o2.R = (dim==0u) ? 1u : i1.R;
o2.C = (dim==1u) ? 1u : i1.C;
o1.S = o2.S = i1.S;
o1.H = o2.H = i1.H;

//Other prep

//Process
if (i1.T==1u)
{
    float *X, *Y, *V;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 1 (Y)" << endl; return 1; }
    try { V = new float[o2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 2 (V)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::sig2ar_burg_s(Y,V,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,P,m0)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    if (wo2)
    {
        try { ofs2.write(reinterpret_cast<char*>(V),o2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (V)" << endl; return 1; }
    }
    delete[] X; delete[] Y; delete[] V;
}

//Finish
