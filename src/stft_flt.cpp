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
#include "stft_flt.c"

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
    const valarray<size_t> oktypes = {1u,2u};
    const size_t I = 2u, O = 1u;
    ifstream ifs1, ifs2; ofstream ofs1;
    int8_t stdi1, stdi2, stdo1, wo1;
    ioinfo i1, i2, o1;
    size_t L, W, nfft, F;
    double stp, c0;
    int mn0, amp, lg;


    //Description
    string descr;
    descr += "Does STFT (short-term Fourier transform) of univariate X1 using window X2.\n";
    descr += "\n";
    descr += "The window (X2) is made by a generating function (hamming, hann, etc.).\n";
    descr += "The signal (X1) and the window (X2) must be real-valued.\n";
    descr += "\n";
    descr += "Each frame of X1 is windowed (element-wise multiplied) by X2,\n";
    descr += "which is a window vector of length L (generated previously).\n";
    descr += "The FFT is done on each windowed frame, and\n";
    descr += "the real-valued power at each positive FFT freq is returned as the STFT.\n";
    descr += "\n";
    descr += "The output Y has size FxW or WxF, where F is nfft/2+1, \n";
    descr += "and nfft is the next-pow-2 of L, \n";
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
    descr += "Use -w (--nframes) to give W, the number of frames [default=(N-1)/stp].\n";
    descr += "This is a positive int (use less than default to use only part of X).\n";
    descr += "\n";
    descr += "Only after the (floating-point) centers of each frame are set,\n";
    descr += "then the center of each frame is rounded to the nearest integer sample.\n";
    descr += "\n";
    descr += "X is extrapolated with zeros if the first/last frames overlap the edge.\n";
    descr += "\n";
    descr += "The following framing convention is used here:\n";
    descr += "Samples from one frame are contiguous in memory, for row- and col-major.\n";
    descr += "So, if Y is row-major, then it has size W x F; \n";
    descr += "but if Y is col-major, then it has size F x W. \n";
    descr += "\n";
    descr += "Include -z (--zero-mean) to subtract the mean from each frame [default=false].\n";
    descr += "This is applied just after windowing.\n";
    descr += "\n";
    descr += "Include -a (--amplitude) to output amplitude rather than power [default=false].\n";
    descr += "This simply takes the sqrt of each element of Y before output.\n";
    descr += "\n";
    descr += "Include -l (--log) to output log amplitude or power [default=false].\n";
    descr += "This simply takes the log of each element of Y before output.\n";
    descr += "\n";
    descr += "Examples:\n";
    descr += "$ stft_flt -s65 X1 X2 -o Y \n";
    descr += "$ stft_flt X1 X2 > Y \n";
    descr += "$ cat X1 | stft_flt - X2 > Y \n";
    descr += "$ hamming -l401 | stft_flt -s160 X1 > Y \n";
    descr += "$ stft_flt -s160 X1 <(hamming -l401) > Y \n";


    //Argtable
    int nerrs;
    struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X1,X2)");
    struct arg_dbl   *a_c0 = arg_dbln("c","c0","<dbl>",0,1,"center of first frame in samps [default=0.0]");
    struct arg_dbl  *a_stp = arg_dbln("s","stp","<dbl>",0,1,"step size btwn frames [default=160.0]");
    struct arg_int    *a_w = arg_intn("w","nframes","<uint>",0,1,"number of frames [default=(N-1)/stp]");
    struct arg_lit  *a_mnz = arg_litn("z","zero-mean",0,1,"include to zero the mean of each frame [default=false]");
    struct arg_lit  *a_amp = arg_litn("a","amplitude",0,1,"include to output amplitude (sqrt of power) [default=false]");
    struct arg_lit  *a_log = arg_litn("l","log",0,1,"include to output log of amplitude or power [default=false]");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_fi, a_c0, a_stp, a_w, a_mnz, a_amp, a_log, a_fo, a_help, a_end};
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
    stdi2 = (a_fi->count<=1 || strlen(a_fi->filename[1])==0 || strcmp(a_fi->filename[1],"-")==0);
    if (stdi1+stdi2>1) { cerr << progstr+": " << __LINE__ << errstr << "can only use stdin for one input" << endl; return 1; }
    if (stdi1+stdi2>0 && isatty(fileno(stdin))) { cerr << progstr+": " << __LINE__ << errstr << "no stdin detected" << endl; return 1; }


    //Check stdout
    if (a_fo->count>0) { stdo1 = (strlen(a_fo->filename[0])==0 || strcmp(a_fo->filename[0],"-")==0); }
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

    //Get mn0
    mn0 = (a_mnz->count>0);

    //Get amp
    amp = (a_amp->count>0);

    //Get lg
    lg = (a_log->count>0);


    //Checks
    if (i1.iscomplex() || i2.iscomplex()) { cerr << progstr+": " << __LINE__ << errstr << "inputs must be real-valued" << endl; return 1; }
    if (i1.T!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
    if (!i1.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) must be a vector" << endl; return 1; }
    if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) found to be empty" << endl; return 1; }
    if (!i2.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) must be a vector" << endl; return 1; }
    if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) found to be empty" << endl; return 1; }


    //Set output header info
    L = i2.N();
    nfft = 1u;
    while (nfft<L) { nfft *= 2u; }
    F = nfft/2u + 1u;
    o1.F = i1.F; o1.T = i1.T;
    o1.R = (i1.isrowmajor()) ? W : F;
    o1.C = (i1.isrowmajor()) ? F : W;
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
        if (codee::stft_flt_s(Y,X1,X2,i1.N(),L,W,nfft,float(c0),float(stp),mn0,amp,lg))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] X1; delete[] X2; delete[] Y;
    }
    else if (o1.T==2)
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
        if (codee::stft_flt_d(Y,X1,X2,i1.N(),L,W,nfft,double(c0),double(stp),mn0,amp,lg))
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
    

    //Exit
    return ret;
}

