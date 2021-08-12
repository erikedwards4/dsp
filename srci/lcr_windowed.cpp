//Includes
#include "lcr_windowed.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u};
const size_t I = 2u, O = 1u;
size_t dim, W, Lx;
double stp, c0, lvl;
int g;

//Description
string descr;
descr += "Gets level crossing rate (LCR) of X1 along dim.\n";
descr += "Output (Y) has a moving average of LCs using the window X2.\n";
descr += "\n";
descr += "X2 is the output of a window function, e.g. hamming.\n";
descr += "X2 should sum to 1 if the output Y is interpreted as a moving average.\n";
descr += "\n";
descr += "Use -v (--level) to give the level of X to test for [default=0].\n";
descr += "For -v0 [default], this is identical to zero crossings (ZCs).\n";
descr += "\n";
descr += "This uses the framing convention of window_univar_flt:\n";
descr += "\n";
descr += "Use -s (--step) to give the step-size (frame-shift) in samples [default=160].\n";
descr += "This is a positive floating-point value.\n";
descr += "\n";
descr += "Use -c (--c0) to give the center-sample of the first frame [default=0].\n";
descr += "This is a positive floating-point value.\n";
descr += "\n";
descr += "Use -w (--nframes) to give W, the number of frames [default=(N-1)/stp].\n";
descr += "This is a positive int (use less than default to use only part of X).\n";
descr += "\n";
descr += "Note: -c and -w do not usually need to be set (defaults are recommended).\n";
descr += "\n";
descr += "Use -g (--going) to specify positive- or negative-going LCs.\n";
descr += "Use -g0 to detect positive- and negative-going LCs [default].\n";
descr += "Use -g1 to detect only positive-going LCs.\n";
descr += "Use -g-1 to detect only negative-going LCs.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension along which to operate.\n";
descr += "The default is 0 (along cols), unless X1 is a row vector.\n";
descr += "\n";
descr += "This is a vec-to-vec operation: \n";
descr += "each vec in X1 is replaced by a vec of length W in Y.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ lcr_windowed -v0.5 X1 X2 -o Y \n";
descr += "$ lcr_windowed -d1 -v-1 -g1 -s5 X1 X2 > Y \n";
descr += "$ hamming -l227 | lcr_windowed -v1e-5 -g-1 -s3 X1 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X1,X2)");
struct arg_dbl  *a_lvl = arg_dbln("v","level","<dbl>",0,1,"level to test for [default=0.0]");
struct arg_dbl   *a_c0 = arg_dbln("c","c0","<dbl>",0,1,"center of first frame in samps [default=0.0]");
struct arg_dbl  *a_stp = arg_dbln("s","stp","<dbl>",0,1,"step size btwn frames [default=160.0]");
struct arg_int    *a_w = arg_intn("w","nframes","<uint>",0,1,"number of frames [default=(N-1)/stp]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to operate [default=0]");
struct arg_int    *a_g = arg_intn("g","going","<uint>",0,1,"if using positive- or negative-going LCs [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = i1.isrowvec() ? 1u : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Get c0
c0 = (a_c0->count>0) ? a_c0->dval[0] : 0.0;
if (c0>double(i1.N()-1u)) { cerr << progstr+": " << __LINE__ << errstr << "c0 (center of first frame) must be <= N-1" << endl; return 1; }

//Get stp
stp = (a_stp->count>0) ? a_stp->dval[0] : 160.0;
if (stp<double(FLT_EPSILON)) { cerr << progstr+": " << __LINE__ << errstr << "stp (step size) must be positive" << endl; return 1; }

//Get W
if (a_w->count==0) { W = size_t((double(i1.N()-1u)-c0)/stp); }
else if (a_w->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "W (nframes) must be positive" << endl; return 1; }
else { W = size_t(a_w->ival[0]); }

//Get g
g = (a_g->count>0) ? a_g->ival[0] : 0;
if (g!=0 && g!=1 && g!=-1) { cerr << progstr+": " << __LINE__ << errstr << "g must be in {-1,0,1}" << endl; return 1; }

//Get lvl
lvl = (a_lvl->count>0) ? a_lvl->dval[0] : 0.0;

//Checks
Lx = (dim==0u) ? i1.R : (dim==1u) ? i1.C : (dim==2u) ? i1.S : i1.H;
if (i1.T!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) found to be empty" << endl; return 1; }
if (!i2.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) must be a vector" << endl; return 1; }
if (i2.N()>=Lx) { cerr << progstr+": " << __LINE__ << errstr << "X2 (win) length must be < length of vecs in X" << endl; return 1; }
if (Lx<2u) { cerr << progstr+": " << __LINE__ << errstr << "cannot work along a singleton dimension" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0u) ? W : i1.R;
o1.C = (dim==1u) ? W : i1.C;
o1.S = (dim==2u) ? W : i1.S;
o1.H = (dim==3u) ? W : i1.H;

//Other prep

//Process
if (i1.T==1u)
{
    float *X1, *X2, *Y;
    try { X1 = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
    try { X2 = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
    if (codee::lcr_windowed_s(Y,X1,X2,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),i2.N(),W,dim,float(c0),float(stp),g,float(lvl)))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X1; delete[] X2; delete[] Y;
}

//Finish
