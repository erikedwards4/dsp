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
#include "smooth_diffdiff.c"

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
    const size_t O = 1u;
    ofstream ofs1;
    int8_t stdo1, wo1;
    ioinfo o1;
    size_t dim, n, N;


    //Description
    string descr;
    descr += "Gets filter coefficients for the 2nd-order smooth differentiators of Holoborodko [2008].\n";
    descr += "This is a generating function (no inputs except opts and params).\n";
    descr += "These are MA filters; the resulting coeffs (B) can be used with fir.\n";
    descr += "\n";
    descr += "These are symmetric filters that are designed to allow exactness of\n";
    descr += "2nd-diff for polynomials of degree n, but smoothness for noise robustness.\n";
    descr += "\n";
    descr += "Use -n (--degree) to give the polynomial degree in {3,5} [default=3]\n";
    descr += "\n";
    descr += "Use -l (--N) to give the filter length in {5,7,9} [default=7]\n";
    descr += "For n=5, N must be in {7,9}.\n";
    descr += "\n";
    descr += "Use -d (--dim) to give the nonsingleton dim of the output vec.\n";
    descr += "If d=0, then Y is a column vector [default].\n";
    descr += "If d=1, then Y is a row vector.\n";
    descr += "(d=2 and d=3 are also possible, but rarely used.)\n";
    descr += "\n";
    descr += "Since this is a generating function (no inputs), the output data type\n";
    descr += "and file format can be specified by -t and -f, respectively. \n";
    descr += "\n";
    descr += "Examples:\n";
    descr += "$ smooth_diffdiff -l9 -o B \n";
    descr += "$ smooth_diffdiff -n5 -l7 > B \n";
    descr += "$ smooth_diffdiff -l7 -d1 -t1 -f101 > B \n";


    //Argtable
    int nerrs;
    struct arg_int    *a_n = arg_intn("n","degree","<uint>",0,1,"polynomial degree [default=3]");
    struct arg_int    *a_l = arg_intn("l","N","<uint>",0,1,"output (filter coeffs) length [default=7]");
    struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"nonsingleton dimension [default=0 -> col vec]");
    struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=1]");
    struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=147]");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (B)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_n, a_l, a_d, a_otyp, a_ofmt, a_fo, a_help, a_end};
    if (arg_nullcheck(argtable)!=0) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating argtable" << endl; return 1; }
    nerrs = arg_parse(argc, argv, argtable);
    if (a_help->count>0)
    {
        cout << "Usage: " << progstr; arg_print_syntax(stdout, argtable, "\n");
        cout << endl; arg_print_glossary(stdout, argtable, "  %-25s %s\n");
        cout << endl << descr; return 1;
    }
    if (nerrs>0) { arg_print_errors(stderr,a_end,(progstr+": "+to_string(__LINE__)+errstr).c_str()); return 1; }


    //Check stdout
    if (a_fo->count>0) { stdo1 = (strlen(a_fo->filename[0])==0u || strcmp(a_fo->filename[0],"-")==0); }
    else { stdo1 = (!isatty(fileno(stdout))); }
    wo1 = (stdo1 || a_fo->count>0);


    //Get options

    //Get o1.F
    if (a_ofmt->count==0) { o1.F = 147u; }
    else if (a_ofmt->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be nonnegative" << endl; return 1; }
    else if (a_ofmt->ival[0]>255) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be < 256" << endl; return 1; }
    else { o1.F = size_t(a_ofmt->ival[0]); }

    //Get o1.T
    if (a_otyp->count==0) { o1.T = 1u; }
    else if (a_otyp->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "data type must be positive int" << endl; return 1; }
    else { o1.T = size_t(a_otyp->ival[0]); }
    if ((o1.T==oktypes).sum()==0)
    {
        cerr << progstr+": " << __LINE__ << errstr << "output data type must be in " << "{";
        for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1u]) ? "}" : ","); }
        cerr << endl; return 1;
    }

    //Get dim
    if (a_d->count==0) { dim = 0u; }
    else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
    else { dim = size_t(a_d->ival[0]); }
    if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

    //Get n
    if (a_n->count==0) { n = 3u; }
    else if (a_n->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "n must be nonnegative" << endl; return 1; }
    else { n = size_t(a_n->ival[0]); }
    if (n!=3u && n!=5u) { cerr << progstr+": " << __LINE__ << errstr << "n (polynomial degree) must be in {3,5}" << endl; return 1; }

    //Get N
    if (a_l->count==0) { N = 7u; }
    else if (a_l->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "N must be nonnegative" << endl; return 1; }
    else { N = size_t(a_l->ival[0]); }
    if (N%2u==0u) { cerr << progstr+": " << __LINE__ << errstr << "N (filter length) must be odd" << endl; return 1; }
    if (N<5u || N>9u) { cerr << progstr+": " << __LINE__ << errstr << "N (filter length) must be in {5,7,9}" << endl; return 1; }


    //Checks
    if (n==5u && N==5u) { cerr << progstr+": " << __LINE__ << errstr << "N (filter length) must be in {7,9} for n=5" << endl; return 1; }


    //Set output header info
    o1.R = (dim==0u) ? N : 1u;
    o1.C = (dim==1u) ? N : 1u;
    o1.S = (dim==2u) ? N : 1u;
    o1.H = (dim==3u) ? N : 1u;


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
        float *Y;
        try { Y = new float[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (B)" << endl; return 1; }
        if (codee::smooth_diffdiff_s(Y,N,n))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (B)" << endl; return 1; }
        }
        delete[] Y;
    }
    else if (o1.T==2u)
    {
        double *Y;
        try { Y = new double[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (B)" << endl; return 1; }
        if (codee::smooth_diffdiff_d(Y,N,n))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (B)" << endl; return 1; }
        }
        delete[] Y;
    }
    else
    {
        cerr << progstr+": " << __LINE__ << errstr << "data type not supported" << endl; return 1;
    }
    
    //Close fstreams
    ofs1.close();

    //Exit
    return ret;
}
