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
#include "gauss.c"

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
    const size_t O = 1u;
    ofstream ofs1;
    int8_t stdo1, wo1;
    ioinfo o1;
    size_t L, dim, norm;
    double a;


    //Description
    string descr;
    descr += "Gets 1D Gauss window with stdev parameter a.\n";
    descr += "\n";
    descr += "Use -l (--winlength) to give the window length in number of taps (sample points).\n";
    descr += "\n";
    descr += "Use -a (--alpha) to give the (positive) stdev parameter [default=2.5].\n";
    descr += "The relation is: stdev = (L-1)/(2*alpha). \n";
    descr += "                 alpha = (L-1)/(2*stdev). \n";
    descr += "where L is the winlength. \n";
    descr += "\n";
    descr += "Use -d (--dim) to give the nonsingleton dim of the output vec.\n";
    descr += "If d=0, then Y is a column vector [default].\n";
    descr += "If d=1, then Y is a row vector.\n";
    descr += "\n";
    descr += "Use -n (--norm) to normalize (divide) the output by a norm:\n";
    descr += "Use -n0 to do no normalization [default].\n";
    descr += "Use -n1 to normalize by the L1-norm (sum abs values).\n";
    descr += "Use -n2 to normalize by the L2-norm (root-sum-of-squares).\n";
    descr += "Use -n3 to normalize by the Inf-norm (maximum).\n";
    descr += "\n";
    descr += "Since this is a generating function (no inputs),\n";
    descr += "the output data type and file format can be specified by\n";
    descr += "-t and -f, respectively (these are the usual CMLI opts).\n";
    descr += "\n";
    descr += "Examples:\n";
    descr += "$ gauss -l255 -o Y \n";
    descr += "$ gauss -l255 -a3 > Y \n";
    descr += "$ gauss -l127 -a3 -d1 -t1 -f101 > Y \n";


    //Argtable
    int nerrs;
    struct arg_int   *a_wl = arg_intn("l","winlength","<uint>",0,1,"window length [default=7]");
    struct arg_dbl    *a_a = arg_dbln("a","alpha","<dbl>",0,1,"alpha stdev param [default=2.5]");
    struct arg_int  *a_nrm = arg_intn("n","norm","<uint>",0,1,"normalize output by norm [default=0]");
    struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"nonsingleton dimension [default=0 -> col vec]");
    struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=2 -> double]");
    struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=102 -> colmajor]");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_wl, a_a, a_nrm, a_d, a_otyp, a_ofmt, a_fo, a_help, a_end};
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
    if (a_fo->count>0) { stdo1 = (strlen(a_fo->filename[0])==0 || strcmp(a_fo->filename[0],"-")==0); }
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
        for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1]) ? "}" : ","); }
        cerr << endl; return 1;
    }

    //Get L
    if (a_wl->count==0) { L = 7u; }
    else if (a_wl->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "winlength must be a positive int" << endl; return 1; }
    else { L = size_t(a_wl->ival[0]); }

    //Get a
    if (a_a->count==0) { a = 2.5; }
    else if (a_a->dval[0]<=0.0) { cerr << progstr+": " << __LINE__ << errstr << "stdev param (a) must be positive" << endl; return 1; }
    else { a = a_a->dval[0]; }

    //Get dim
    if (a_d->count==0) { dim = 0u; }
    else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
    else { dim = size_t(a_d->ival[0]); }
    if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

    //Get norm
    if (a_nrm->count==0) { norm = 0u; }
    else if (a_nrm->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "norm must be nonnegative" << endl; return 1; }
    else { norm = size_t(a_nrm->ival[0]); }
    if (norm>3) { cerr << progstr+": " << __LINE__ << errstr << "norm must be in {0,1,2,3}" << endl; return 1; }


    //Set output header info
    o1.R = (dim==0u) ? L : 1u;
    o1.C = (dim==1u) ? L : 1u;
    o1.S = (dim==2u) ? L : 1u;
    o1.H = (dim==3u) ? L : 1u;


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
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        if (codee::gauss_s(Y,L,float(a),norm))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] Y;
    }
    else if (o1.T==2)
    {
        double *Y;
        try { Y = new double[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        if (codee::gauss_d(Y,L,double(a),norm))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] Y;
    }
    else if (o1.T==101u)
    {
        float *Y;
        try { Y = new float[2u*o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        if (codee::gauss_c(Y,L,float(a),norm))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] Y;
    }
    else if (o1.T==102u)
    {
        double *Y;
        try { Y = new double[2u*o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        if (codee::gauss_z(Y,L,double(a),norm))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] Y;
    }
    else
    {
        cerr << progstr+": " << __LINE__ << errstr << "data type not supported" << endl; return 1;
    }
    

    //Exit
    return ret;
}

