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
#include "fft.fftw.c"

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
    size_t dim, nfft, Ly;
    char sc;


    //Description
    string descr;
    descr += "1D FFT (fast Fourier transform) of each vector (1D signal) in X.\n";
    descr += "\n";
    descr += "Use -d (--dim) to give the dimension along which to transform.\n";
    descr += "Use -d0 to operate along cols, -d1 to operate along rows, etc.\n";
    descr += "The default is 0 (along cols), unless X is a row vector.\n";
    descr += "\n";
    descr += "Use -n (--nfft) to specify transform length [default=L].\n";
    descr += "The default (L) is the length of X along dim.\n";
    descr += "X is zero-padded as necessary to match nfft.\n";
    descr += "\n";
    descr += "The output (Y) is complex-valued with length nfft along dim\n";
    descr += "for complex X, and length nfrqs along dim for real X, \n";
    descr += "where nfrqs = floor(nfft/2)+1 = num nonnegative FFT frequencies.\n";
    descr += "\n";
    descr += "Note: to get same result + negative freqs, just convert X to complex.\n";
    descr += "Alternately, use fft.rad2, which outputs all nfft freqs.\n";
    descr += "\n";
    descr += "Include -s (--scale) to scale by sqrt(0.5/L), for formal definition.\n";
    descr += "\n";
    descr += "Examples:\n";
    descr += "$ fft.fftw -n256 X -o Y \n";
    descr += "$ fft.fftw -n256 -d1 X > Y \n";
    descr += "$ cat X | fft.fftw -n256 > Y \n";


    //Argtable
    int nerrs;
    struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
    struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to transform [default=0]");
    struct arg_int    *a_n = arg_intn("n","nfft","<uint>",0,1,"transform length [default=L]");
    struct arg_lit   *a_sc = arg_litn("s","scale",0,1,"include to scale by sqrt(0.5/n) [default=no scale]");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_fi, a_d, a_n, a_sc, a_fo, a_help, a_end};
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

    //Get dim
    if (a_d->count==0) { dim = i1.isrowvec() ? 1u : 0u; }
    else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
    else { dim = size_t(a_d->ival[0]); }
    if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

    //Get nfft
    if (a_n->count==0) { nfft = (dim==0u) ? i1.R : (dim==1u) ? i1.C : (dim==2u) ? i1.S : i1.H; }
    else if (a_n->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be positive" << endl; return 1; }
    else { nfft = size_t(a_n->ival[0]); }

    //Get sc
    sc = (a_sc->count>0);


    //Checks
    if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
    if (dim==0u && nfft<i1.R) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be >= nrows X for dim=0" << endl; return 1; }
    if (dim==1u && nfft<i1.C) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be >= ncols X for dim=1" << endl; return 1; }
    if (dim==2u && nfft<i1.S) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be >= nslices X for dim=2" << endl; return 1; }
    if (dim==3u && nfft<i1.H) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be >= nhyperslices X for dim=3" << endl; return 1; }


    //Set output header info
    o1.F = i1.F;
    o1.T = i1.isreal() ? i1.T+100u : i1.T;
    Ly = i1.isreal() ? nfft/2u+1u : nfft;
    //Ly = i1.isreal() ? nfft : nfft;
    o1.R = (dim==0u) ? Ly : i1.R;
    o1.C = (dim==1u) ? Ly : i1.C;
    o1.S = (dim==2u) ? Ly : i1.S;
    o1.H = (dim==3u) ? Ly : i1.H;


    //Open output
    if (wo1)
    {
        if (stdo1) { ofs1.copyfmt(cout); ofs1.basic_ios<char>::rdbuf(cout.rdbuf()); } else { ofs1.open(a_fo->filename[0]); }
        if (!ofs1) { cerr << progstr+": " << __LINE__ << errstr << "problem opening output file 1" << endl; return 1; }
    }


    //Write output header
    if (wo1 && !write_output_header(ofs1,o1)) { cerr << progstr+": " << __LINE__ << errstr << "problem writing header for output file 1" << endl; return 1; }


    //Other prep
    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    

    //Process
    if (i1.T==1u)
    {
        float *X, *Y;
        try { X = new float[i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
        try { Y = new float[2u*o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
        if (codee::fft_fftw_s(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,nfft,sc))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] X; delete[] Y;
    }
    else if (i1.T==2)
    {
        double *X, *Y;
        try { X = new double[i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
        try { Y = new double[2u*o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
        if (codee::fft_fftw_d(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,nfft,sc))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] X; delete[] Y;
    }
    else if (i1.T==101u)
    {
        float *X, *Y;
        try { X = new float[2u*i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
        try { Y = new float[2u*o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
        if (codee::fft_fftw_c(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,nfft,sc))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] X; delete[] Y;
    }
    else if (i1.T==102u)
    {
        double *X, *Y;
        try { X = new double[2u*i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
        try { Y = new double[2u*o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
        if (codee::fft_fftw_z(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,nfft,sc))
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
    

    //Finish
    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);


    //Exit
    return ret;
}

