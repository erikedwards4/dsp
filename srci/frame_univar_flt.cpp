//Includes
#include <cfloat>
#include "frame_univar_flt.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 1u, O = 1u;
size_t L, W;
double stp, c0;

//Description
string descr;
descr += "Takes univariate X and produces a series of (overlapping) frames.\n";
descr += "The output Y has size LxW or WxL, where L is the length of each frame, \n";
descr += "and W is the number of frames (a.k.a. windows).\n";
descr += "\n";
descr += "This _float version has different options and conventions;\n";
descr += "and allows float (non-integer) values for tep size and start samp.\n";
descr += "\n";
descr += "Use -s (--step) to give the step-size (frame-shift) in samples [default=160].\n";
descr += "This is a positive floating-point value.\n";
descr += "\n";
descr += "Use -c (--c0) to give the center-sample of the first frame [default=0].\n";
descr += "This is a positive floating-point value.\n";
descr += "\n";
descr += "Use -l (--winlength) to give L, the length of each frame [default=401].\n";
descr += "This is a positive int (because a window is a vector of int length).\n";
descr += "\n";
descr += "Use -w (--nframes) to give W, the number of frames [default=(N-1)/stp].\n";
descr += "This is a positive int (can use less than default to use only part of X).\n";
descr += "\n";
descr += "Only after the (floating-point) centers of each frame are set,\n";
descr += "then the center of each frame is rounded to the nearest integer sample.\n";
descr += "\n";
descr += "X is extrapolated with zeros if the first/last frames overlap the edge.\n";
descr += "\n";
descr += "The following framing convention is used here:\n";
descr += "Samples from one frame are contiguous in memory, for row- and col-major.\n";
descr += "So, if Y is row-major, then it has size W x L; \n";
descr += "but if Y is col-major, then it has size L x W. \n";
descr += "\n";
descr += "Examples:\n";
descr += "$ frame_univar_flt -l255 -s65 -c0 X -o Y \n";
descr += "$ frame_univar_flt -l255 -w20 -c60.25 -s60.25 X > Y \n";
descr += "$ cat X | frame_univar_flt -l127 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_dbl   *a_c0 = arg_dbln("c","c0","<dbl>",0,1,"center of first frame in samps [default=0.0]");
struct arg_dbl  *a_stp = arg_dbln("s","stp","<dbl>",0,1,"step size btwn frames [default=160.0]");
struct arg_int    *a_l = arg_intn("l","winlength","<uint>",0,1,"length in samps of each frame [default=401]");
struct arg_int    *a_w = arg_intn("w","nframes","<uint>",0,1,"number of frames [default=(N-1)/stp]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get c0
c0 = (a_c0->count>0) ? a_c0->dval[0] : 0.0;
if (c0>double(i1.N()-1u)) { cerr << progstr+": " << __LINE__ << errstr << "c0 (center of first frame) must be <= N-1" << endl; return 1; }

//Get stp
stp = (a_stp->count>0) ? a_stp->dval[0] : 160.0;
if (stp<double(FLT_EPSILON)) { cerr << progstr+": " << __LINE__ << errstr << "stp (step size) must be positive" << endl; return 1; }

//Get L
if (a_l->count==0) { L = 401u; }
else if (a_l->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "L (winlength) must be positive" << endl; return 1; }
else { L = size_t(a_l->ival[0]); }

//Get W
if (a_w->count==0) { W = size_t((double(i1.N()-1u)-c0)/stp); }
else if (a_w->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "W (nframes) must be positive" << endl; return 1; }
else { W = size_t(a_w->ival[0]); }

//Checks
if (!i1.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) must be a vector" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }

//Set output header
o1.F = i1.F; o1.T = i1.T;
o1.R = (i1.isrowmajor()) ? W : L;
o1.C = (i1.isrowmajor()) ? L : W;
o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (o1.T==1u)
{
    float *X, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::frame_univar_flt_s(Y,X,i1.N(),L,W,float(c0),float(stp)))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}
else if (o1.T==101u)
{
    float *X, *Y;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { Y = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::frame_univar_flt_c(Y,X,i1.N(),L,W,float(c0),float(stp)))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish
