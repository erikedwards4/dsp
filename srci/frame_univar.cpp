//Includes
#include <cfloat>
#include "frame_univar.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 1u, O = 1u;
size_t L, stp, W;
int snip_edges;

//Description
string descr;
descr += "Takes univariate X and produces a series of (overlapping) frames.\n";
descr += "The output Y has size LxW or WxL, where L is the length of each frame, \n";
descr += "and W is the number of frames (a.k.a. windows).\n";
descr += "\n";
descr += "Use -l (--winlength) to give L, the length of each frame [default=1].\n";
descr += "\n";
descr += "Use -s (--step) to give the step-size (frame-shift) in samples [default=1].\n";
descr += "\n";
descr += "Use -e (--snip-edges) to set snip-edges to true [default=false].\n";
descr += "This is a setting from HTK, Kaldi, Librosa, etc., which controls\n";
descr += "the placement of the first/last frames w.r.t. the start/end of X.\n";
descr += "This is used here for compatibility.\n";
descr += "\n";
descr += "The number of output frames (W) is set as in Kaldi:\n";
descr += "If snip-edges=true:  W = 1u + (N-L)/stp   \n";
descr += "If snip-edges=false: W = (N+stp/2u) / stp \n";
descr += "\n";
descr += "If snip-edges=true, the first frame starts at samp 0,\n";
descr += "and the last frame fits entirely within the length of X..\n";
descr += "If snip-edges=false, the first frame is centered at samp stp/2,\n";
descr += "and the last frame can overlap the end of X.\n";
descr += "\n";
descr += "Also following Kaldi for compatibility, X is extrapolated \n";
descr += "by reversing the edge samples of X, if snip-edges=false. \n";
descr += "\n";
descr += "The following framing convention is used here:\n";
descr += "Samples from one frame are contiguous in memory, for row- and col-major.\n";
descr += "So, if Y is row-major, then it has size W x L; \n";
descr += "but if Y is col-major, then it has size L x W. \n";
descr += "\n";
descr += "Examples:\n";
descr += "$ frame_univar -l255 -s65 X -o Y \n";
descr += "$ frame_univar -l255 -e X > Y \n";
descr += "$ cat X | frame_univar -l127 -e > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int   *a_wl = arg_intn("l","winlength","<uint>",0,1,"length in samps of each frame [default=1]");
struct arg_int  *a_stp = arg_intn("s","step","<uint>",0,1,"step in samps between each frame [default=1]");
struct arg_lit  *a_sne = arg_litn("e","snip-edges",0,1,"include to snip edges [default=false]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get L
if (a_wl->count==0) { L = 1u; }
else if (a_wl->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "L must be positive" << endl; return 1; }
else { L = size_t(a_wl->ival[0]); }

//Get stp
if (a_stp->count==0) { stp = 1u; }
else if (a_stp->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "stp must be positive" << endl; return 1; }
else { stp = size_t(a_stp->ival[0]); }

//Get snip_edges
snip_edges = (a_sne->count>0);

//Checks
if (!i1.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) must be a vector" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }

//Set output header
W = (snip_edges) ? 1u+(i1.N()-L)/stp : (i1.N()+stp/2u)/stp;
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
    if (codee::frame_univar_s(Y,X,i1.N(),L,stp,snip_edges))
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
    if (codee::frame_univar_c(Y,X,i1.N(),L,stp,snip_edges))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish
