//Includes
#include "lcr.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u};
const size_t I = 1u, O = 1u;
size_t dim, Lx, Lw;
int g, causal;
double lvl;

//Description
string descr;
descr += "Gets the level crossing rate (LCR) of X along dim.\n";
descr += "This uses the usual definition (i.e., rectangular window of length Lw).\n";
descr += "\n";
descr += "Use -l (--Lw) to give the length of the moving-average window.\n";
descr += "\n";
descr += "Use -v (--level) to give the level of X to test for [default=0].\n";
descr += "For -v0 [default], this is identical to zero crossings (ZCs).\n";
descr += "\n";
descr += "Use -g (--going) to specify positive- or negative-going LCs.\n";
descr += "Use -g0 to detect positive- and negative-going LCs [default].\n";
descr += "Use -g1 to detect only positive-going LCs.\n";
descr += "Use -g-1 to detect only negative-going LCs.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension along which to operate.\n";
descr += "The default is 0 (along cols), unless X is a row vector.\n";
descr += "\n";
descr += "Include -c (--causal) to give causal output [default=false].\n";
descr += "The default (noncausal) uses a centered moving-average window.\n";
descr += "The causal option averages only over the current and past samps.\n";
descr += "\n";
descr += "Y is real-valued and has the same size and file format as X.\n";
descr += "Y is in units of LCs/sample (so multiply by fs to get LCs/sec).\n";
descr += "Thus, each element of Y is in [0.0 1.0].\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ lcr -l100 -v0.1 X -o Y \n";
descr += "$ lcr -d1 -l7 -v-8 X > Y \n";
descr += "$ cat X | lcr -l127 -v1e-5 -g1 -c > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_dbl  *a_lvl = arg_dbln("v","level","<dbl>",0,1,"level to test for [default=0.0]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to operate [default=0]");
struct arg_int    *a_l = arg_intn("l","Lw","<uint>",0,1,"winlength over which to average LCs [default=1]");
struct arg_int    *a_g = arg_intn("g","going","<uint>",0,1,"if using positive- or negative-going LCs [default=0]");
struct arg_lit    *a_c = arg_litn("c","causal",0,1,"include to give causal output [default=false]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = i1.isrowvec() ? 1u : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Get g
g = (a_g->count>0) ? a_g->ival[0] : 0;
if (g!=0 && g!=1 && g!=-1) { cerr << progstr+": " << __LINE__ << errstr << "g must be in {-1,0,1}" << endl; return 1; }

//Get Lw
if (a_l->count==0) { Lw = 1u; }
else if (a_l->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "Lw must be positive" << endl; return 1; }
else { Lw = size_t(a_l->ival[0]); }

//Get g
g = (a_g->count>0) ? a_g->ival[0] : 0;
if (g!=0 && g!=1 && g!=-1) { cerr << progstr+": " << __LINE__ << errstr << "g must be in {-1,0,1}" << endl; return 1; }

//Get lvl
lvl = (a_lvl->count>0) ? a_lvl->dval[0] : 0.0;

//Get causal
causal = (a_c->count>0);

//Checks
Lx = (dim==0u) ? i1.R : (dim==1u) ? i1.C : (dim==2u) ? i1.S : i1.H;
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
if (Lx<2u) { cerr << progstr+": " << __LINE__ << errstr << "cannot work along a singleton dimension" << endl; return 1; }
if (Lx<Lw) { cerr << progstr+": " << __LINE__ << errstr << "Lw (winlength) must be <= Lx (length of vecs in X)" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = i1.R; o1.C = i1.C; o1.S = i1.S; o1.H = i1.H;

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
    if (codee::lcr_s(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,Lw,g,float(lvl),causal))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish
