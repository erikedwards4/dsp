//Includes
#include "pow_compress.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u};
const size_t I = 1u, O = 1u;
double p, preg;

//Description
string descr;
descr += "Elementwise function for 1 input.\n";
descr += "Applies power compression to each element in X.\n";
descr += "This is intended for nonnegative X only (such as power).\n";
descr += "\n";
descr += "A small regularization (preg) is added before compression.\n";
descr += "\n";
descr += "Use -r (--preg) to give the small additive value [default=FLT_EPS].\n";
descr += "preg must be nonnegative.\n";
descr += "\n";
descr += "Use -p (--pow) to give the power exponent [default=1/3].\n";
descr += "p must be in [0 1].\n";
descr += "\n";
descr += "If p=1, then Y = X + preg.\n";
descr += "If p=0, then Y = log(X+preg).\n";
descr += "If p=1/2, then Y = sqrt(X+preg).\n";
descr += "If p=1/3, then Y = cbrt(X+preg).\n";
descr += "Otherwise, then Y = (X+preg)^p.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ pow_compress -p0.5 X -o Y \n";
descr += "$ pow_compress -p0 -r1e-5 X > Y \n";
descr += "$ cat X | pow_compress -p0.1 - > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_dbl  *a_pow = arg_dbln("p","pow","<dbl>",0,1,"power exponent (0 <= pow <= 1) [default=1/3]");
struct arg_dbl  *a_reg = arg_dbln("r","preg","<dbl>",0,1,"power regularizer [default=FLT_EPS]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get p
p = (a_pow->count>0) ? a_pow->dval[0] : 1.0/3.0;
if (p<0.0 || p>1.0) { cerr << progstr+": " << __LINE__ << errstr << "pow must be in [0.0 1.0]" << endl; return 1; }

//Get preg
preg = (a_reg->count>0) ? a_reg->dval[0] : double(FLT_EPSILON);
if (preg<0.0) { cerr << progstr+": " << __LINE__ << errstr << "preg must be nonnegative" << endl; return 1; }

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = i1.R; o1.C = i1.C;
o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1u)
{
    float *X;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::pow_compress_inplace_s(X,i1.N(),float(p),float(preg)))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(X),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X;
}

//Finish
