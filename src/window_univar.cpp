//@author Erik Edwards
//@date 2018-present
//@license BSD 3-clause


#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <cstring>
#include <valarray>
#include <unordered_map>
#include <argtable2.h>
#include "../util/cmli.hpp"
#include <cfloat>
#include "window_univar.c"

#ifdef I
#undef I
#endif


int main(int argc, char *argv[])
{
    using namespace std;


    //Declarations
    int ret = 0;
    const string errstr = ": \033[1;31merror:\033[0m ";
    const string warstr = ": \033[1;35mwarning:\033[0m ";
    const string progstr(__FILE__,string(__FILE__).find_last_of("/")+1,strlen(__FILE__)-string(__FILE__).find_last_of("/")-5);
    const valarray<size_t> oktypes = {1u,2u,101u,102u};
    const size_t I = 2u, O = 1u;
    ifstream ifs1, ifs2; ofstream ofs1;
    int8_t stdi1, stdi2, stdo1, wo1;
    ioinfo i1, i2, o1;
    size_t L, stp, W;
    int snip_edges;


    //Description
    string descr;
    descr += "Windows univariate X1 to produces a series of (overlapping) frames.\n";
    descr += "\n";
    descr += "Each frame of X1 is windowed (element-wise multiplied) by X2,\n";
    descr += "which is a window vector of length L (generated previously).\n";
    descr += "\n";
    descr += "The output Y has size LxW or WxL, where L is the length of each frame, \n";
    descr += "and W is the number of frames (a.k.a. windows).\n";
    descr += "\n";
    descr += "The window (X2) is made by a generating function (hamming, hann, etc.).\n";
    descr += "The window (X2) must be real-valued.\n";
    descr += "\n";
    descr += "Other than the window input (X2), this is identical to frame_univar.\n";
    descr += "\n";
    descr += "Use -s (--step) to give the step-size (frame-shift) in samples [default=1].\n";
    descr += "\n";
    descr += "Use -e (--snip-edges) to set snip-edges to true [default=false].\n";
    descr += "This is a setting from HTK, Kaldi, Librosa, etc., which controls\n";
    descr += "the placement of the first/last frames w.r.t. the start/end of X1.\n";
    descr += "This is used here for compatibility.\n";
    descr += "\n";
    descr += "The number of output frames (W) is set as in Kaldi:\n";
    descr += "If snip-edges=true:  W = 1u + (N-L)/stp   \n";
    descr += "If snip-edges=false: W = (N+stp/2u) / stp \n";
    descr += "\n";
    descr += "If snip-edges=true, the first frame starts at samp 0,\n";
    descr += "and the last frame fits entirely within the length of X.\n";
    descr += "If snip-edges=false, the first frame is centered at samp stp/2,\n";
    descr += "and the last frame can overlap the end of X.\n";
    descr += "\n";
    descr += "Also following Kaldi for compatibility, X1 is extrapolated \n";
    descr += "by reversing the edge samples of X1, if snip-edges=false. \n";
    descr += "\n";
    descr += "The following framing convention is used here:\n";
    descr += "Samples from one frame are contiguous in memory, for row- and col-major.\n";
    descr += "So, if Y is row-major, then it has size W x L; \n";
    descr += "but if Y is col-major, then it has size L x W. \n";
    descr += "\n";
    descr += "Examples:\n";
    descr += "$ window_univar -s65 X1 X2 -o Y \n";
    descr += "$ window_univar -e X1 X2 > Y \n";
    descr += "$ cat X1 | window_univar -e - X2 > Y \n";
    descr += "$ hamming -l401 | window_univar -s160 X1 - > Y \n";
    descr += "$ window_univar -s160 X1 <(hamming -l401) > Y \n";


    //Argtable
    int nerrs;
    struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X1,X2)");
    struct arg_int  *a_stp = arg_intn("s","step","<uint>",0,1,"step in samps between each frame [default=1]");
    struct arg_lit  *a_sne = arg_litn("e","snip-edges",0,1,"include to snip edges [default=false]");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_fi, a_stp, a_sne, a_fo, a_help, a_end};
    if (arg_nullcheck(argtable)!=0) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating argtable" << endl; return 1; }
    nerrs = arg_parse(argc, argv, argtable);
    if (a_help->count>0)
    {
        cout << "Usage: " << progstr; arg_print_syntax(stdout, argtable, "\n");
        cout << endl; arg_print_glossary(stdout, argtable, "  %-25s %s\n");
        cout << endl << descr; return 1;
    }
    if (nerrs>0) { arg_print_errors(stderr,a_end,(progstr+": "+to_string(__LINE__)+errstr).c_str()); return 1; }


    //Check stdin
    stdi1 = (a_fi->count==0 || strlen(a_fi->filename[0])==0u || strcmp(a_fi->filename[0],"-")==0);
    stdi2 = (a_fi->count<=1 || strlen(a_fi->filename[1])==0u || strcmp(a_fi->filename[1],"-")==0);
    if (stdi1+stdi2>1) { cerr << progstr+": " << __LINE__ << errstr << "can only use stdin for one input" << endl; return 1; }
    if (stdi1+stdi2>0 && isatty(fileno(stdin))) { cerr << progstr+": " << __LINE__ << errstr << "no stdin detected" << endl; return 1; }


    //Check stdout
    if (a_fo->count>0) { stdo1 = (strlen(a_fo->filename[0])==0u || strcmp(a_fo->filename[0],"-")==0); }
    else { stdo1 = (!isatty(fileno(stdout))); }
    wo1 = (stdo1 || a_fo->count>0);


    //Open inputs
    if (stdi1) { ifs1.copyfmt(cin); ifs1.basic_ios<char>::rdbuf(cin.rdbuf()); } else { ifs1.open(a_fi->filename[0]); }
    if (!ifs1) { cerr << progstr+": " << __LINE__ << errstr << "problem opening input file 1" << endl; return 1; }
    if (stdi2) { ifs2.copyfmt(cin); ifs2.basic_ios<char>::rdbuf(cin.rdbuf()); } else { ifs2.open(a_fi->filename[1]); }
    if (!ifs2) { cerr << progstr+": " << __LINE__ << errstr << "problem opening input file 2" << endl; return 1; }


    //Read input headers
    if (!read_input_header(ifs1,i1)) { cerr << progstr+": " << __LINE__ << errstr << "problem reading header for input file 1" << endl; return 1; }
    if (!read_input_header(ifs2,i2)) { cerr << progstr+": " << __LINE__ << errstr << "problem reading header for input file 2" << endl; return 1; }
    if ((i1.T==oktypes).sum()==0 || (i2.T==oktypes).sum()==0)
    {
        cerr << progstr+": " << __LINE__ << errstr << "input data type must be in " << "{";
        for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1u]) ? "}" : ","); }
        cerr << endl; return 1;
    }


    //Get options

    //Get stp
    if (a_stp->count==0) { stp = 1u; }
    else if (a_stp->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "stp must be positive" << endl; return 1; }
    else { stp = size_t(a_stp->ival[0]); }

    //Get snip_edges
    snip_edges = (a_sne->count>0);


    //Checks
    if (i1.T!=i2.T && i1.T-100u!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have compatible data types" << endl; return 1; }
    if (!i1.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) must be a vector" << endl; return 1; }
    if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) found to be empty" << endl; return 1; }
    if (!i2.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) must be a vector" << endl; return 1; }
    if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) found to be empty" << endl; return 1; }
    if (i2.iscomplex()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) must be real-valued" << endl; return 1; }


    //Set output header info
    L = i2.N();
    W = (snip_edges) ? 1u+(i1.N()-L)/stp : (i1.N()+stp/2u)/stp;
    o1.F = i1.F; o1.T = i1.T;
    o1.R = (i1.isrowmajor()) ? W : L;
    o1.C = (i1.isrowmajor()) ? L : W;
    o1.S = i1.S; o1.H = i1.H;


    //Open output
    if (wo1)
    {
        if (stdo1) { ofs1.copyfmt(cout); ofs1.basic_ios<char>::rdbuf(cout.rdbuf()); } else { ofs1.open(a_fo->filename[0]); }
        if (!ofs1) { cerr << progstr+": " << __LINE__ << errstr << "problem opening output file 1" << endl; return 1; }
    }


    //Write output header
    if (wo1 && !write_output_header(ofs1,o1)) { cerr << progstr+": " << __LINE__ << errstr << "problem writing header for output file 1" << endl; return 1; }


    //Other prep


    //Process
    if (o1.T==1u)
    {
        float *X1, *X2, *Y;
        try { X1 = new float[i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
        try { X2 = new float[i2.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
        try { Y = new float[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
        if (codee::window_univar_s(Y,X1,X2,i1.N(),L,stp,snip_edges))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] X1; delete[] X2; delete[] Y;
    }
    else if (o1.T==2u)
    {
        double *X1, *X2, *Y;
        try { X1 = new double[i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
        try { X2 = new double[i2.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
        try { Y = new double[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
        if (codee::window_univar_d(Y,X1,X2,i1.N(),L,stp,snip_edges))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] X1; delete[] X2; delete[] Y;
    }
    else if (o1.T==101u)
    {
        float *X1, *X2, *Y;
        try { X1 = new float[2u*i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
        try { X2 = new float[2u*i2.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
        try { Y = new float[2u*o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
        try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
        if (codee::window_univar_c(Y,X1,X2,i1.N(),L,stp,snip_edges))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] X1; delete[] X2; delete[] Y;
    }
    else if (o1.T==102u)
    {
        double *X1, *X2, *Y;
        try { X1 = new double[2u*i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
        try { X2 = new double[2u*i2.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
        try { Y = new double[2u*o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
        try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
        if (codee::window_univar_z(Y,X1,X2,i1.N(),L,stp,snip_edges))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] X1; delete[] X2; delete[] Y;
    }
    else
    {
        cerr << progstr+": " << __LINE__ << errstr << "data type not supported" << endl; return 1;
    }
    
    //Close fstreams
    ifs1.close(); ifs2.close();

    ofs1.close();

    //Exit
    return ret;
}
