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
#include "cmli.hpp"
#include <cfloat>
#include "frame_univar.c"

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
    const size_t I = 1u, O = 1u;
    ifstream ifs1; ofstream ofs1;
    int8_t stdi1, stdo1, wo1;
    ioinfo i1, o1;
    size_t L, stp, W;
    int snip_edges;


    //Description
    string descr;
    descr += "Takes univariate X and produces a series of (overlapping) frames.\n";
    descr += "The output Y has size LxW or WxL, where L is the length of each frame, \n";
    descr += "and W is the number of frames (a.k.a. windows).\n";
    descr += "\n";
    descr += "Use -l (--winlength) to give L, the length of each frame [default=401].\n";
    descr += "\n";
    descr += "Use -s (--step) to give the step-size (frame-shift) in samples [default=160].\n";
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
    int nerrs;
    struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
    struct arg_int   *a_wl = arg_intn("l","winlength","<uint>",0,1,"length in samps of each frame [default=401]");
    struct arg_int  *a_stp = arg_intn("s","step","<uint>",0,1,"step in samps between each frame [default=160]");
    struct arg_lit  *a_sne = arg_litn("e","snip-edges",0,1,"include to snip edges [default=false]");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_fi, a_wl, a_stp, a_sne, a_fo, a_help, a_end};
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
    stdi1 = (a_fi->count==0 || strlen(a_fi->filename[0])==0 || strcmp(a_fi->filename[0],"-")==0);
    if (stdi1>0 && isatty(fileno(stdin))) { cerr << progstr+": " << __LINE__ << errstr << "no stdin detected" << endl; return 1; }


    //Check stdout
    if (a_fo->count>0) { stdo1 = (strlen(a_fo->filename[0])==0 || strcmp(a_fo->filename[0],"-")==0); }
    else { stdo1 = (!isatty(fileno(stdout))); }
    wo1 = (stdo1 || a_fo->count>0);


    //Open input
    if (stdi1) { ifs1.copyfmt(cin); ifs1.basic_ios<char>::rdbuf(cin.rdbuf()); } else { ifs1.open(a_fi->filename[0]); }
    if (!ifs1) { cerr << progstr+": " << __LINE__ << errstr << "problem opening input file" << endl; return 1; }


    //Read input header
    if (!read_input_header(ifs1,i1)) { cerr << progstr+": " << __LINE__ << errstr << "problem reading header for input file" << endl; return 1; }
    if ((i1.T==oktypes).sum()==0)
    {
        cerr << progstr+": " << __LINE__ << errstr << "input data type must be in " << "{";
        for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1u]) ? "}" : ","); }
        cerr << endl; return 1;
    }


    //Get options

    //Get L
    if (a_wl->count==0) { L = 401u; }
    else if (a_wl->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "L must be positive" << endl; return 1; }
    else { L = size_t(a_wl->ival[0]); }

    //Get stp
    if (a_stp->count==0) { stp = 160u; }
    else if (a_stp->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "stp must be positive" << endl; return 1; }
    else { stp = size_t(a_stp->ival[0]); }

    //Get snip_edges
    snip_edges = (a_sne->count>0);


    //Checks
    if (!i1.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) must be a vector" << endl; return 1; }
    if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }


    //Set output header info
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
    else if (o1.T==2)
    {
        double *X, *Y;
        try { X = new double[i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
        try { Y = new double[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        if (codee::frame_univar_d(Y,X,i1.N(),L,stp,snip_edges))
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
    else if (o1.T==102u)
    {
        double *X, *Y;
        try { X = new double[2u*i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
        try { Y = new double[2u*o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        if (codee::frame_univar_z(Y,X,i1.N(),L,stp,snip_edges))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] X; delete[] Y;
    }
    else
    {
        cerr << progstr+": " << __LINE__ << errstr << "data type not supported" << endl; return 1;
    }
    

    //Exit
    return ret;
}

