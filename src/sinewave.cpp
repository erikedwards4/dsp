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
#include "sinewave.c"

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
    size_t N, dim;
    double amp, frq, phs;


    //Description
    string descr;
    descr += "Generates 1-D sinewave signal.\n";
    descr += "\n";
    descr += "Use -n (--N) to give the output vector length in sample points.\n";
    descr += "\n";
    descr += "Use -d (--dim) to give the nonsingleton dim of the output vec.\n";
    descr += "If d=0, then Y is a column vector [default].\n";
    descr += "If d=1, then Y is a row vector.\n";
    descr += "(d=2 and d=3 are also possible, but rarely used.)\n";
    descr += "\n";
    descr += "Use -a (--amplitude) to give the amplitude (max) of the output.\n";
    descr += "Note that this is half of the peak-to-peak amplitude.\n";
    descr += "\n";
    descr += "Use -f (--freq) to give the frequency in units of cycles/sample.\n";
    descr += "This is the frequency in Hz divided by the sample rate in Hz.\n";
    descr += "Thus, f should be in [0 0.5], where 0 is DC and 0.5 is Nyquist.\n";
    descr += "\n";
    descr += "Use -p (--phase) to give the phase in radians.\n";
    descr += "This is usually in [0 2*pi), but the value will be taken modulo 2*pi.\n";
    descr += "\n";
    descr += "Since this is a generating function (no inputs), the output data type\n";
    descr += "and file format can be specified by -t and -f, respectively. \n";
    descr += "\n";
    descr += "Examples:\n";
    descr += "$ sinewave -n32 -a2.5 -f0.2 -o Y \n";
    descr += "$ sinewave -n32 -f0.01 -p1.57079632679 > Y \n";
    descr += "$ sinewave -n32 -f0.01 -d1 -t1 -f101 > Y \n";


    //Argtable
    int nerrs;
    struct arg_dbl  *a_amp = arg_dbln("a","amp","<dbl>",0,1,"amplitude [default=1.0]");
    struct arg_dbl  *a_frq = arg_dbln("f","freq","<dbl>",0,1,"frequency in cycles/sample [default=0.5]");
    struct arg_dbl  *a_phs = arg_dbln("p","phase","<dbl>",0,1,"phase in radians [default=0.0]");
    struct arg_int    *a_n = arg_intn("n","N","<uint>",0,1,"num samples in output [default=1]");
    struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"nonsingleton dimension [default=0 -> col vec]");
    struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=1]");
    struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=147]");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_amp, a_frq, a_phs, a_n, a_d, a_otyp, a_ofmt, a_fo, a_help, a_end};
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
        for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1u]) ? "}" : ","); }
        cerr << endl; return 1;
    }

    //Get N
    if (a_n->count==0) { N = 1u; }
    else if (a_n->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "N must be a positive int" << endl; return 1; }
    else { N = size_t(a_n->ival[0]); }

    //Get dim
    if (a_d->count==0) { dim = 0u; }
    else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
    else { dim = size_t(a_d->ival[0]); }
    if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

    //Get amp
    amp = (a_amp->count>0) ? a_amp->dval[0] : 1.0;
    if (amp<0.0) { cerr << progstr+": " << __LINE__ << errstr << "amplitude must be nonnegative" << endl; return 1; }
    if (o1.T==1u && amp>=double(FLT_MAX)) { cerr << progstr+": " << __LINE__ << errstr << "amplitude must be < " << double(FLT_MAX) << endl; return 1; }

    //Get frq
    frq = (a_frq->count>0) ? a_frq->dval[0] : 0.5;
    if (frq<double(FLT_EPSILON)) { cerr << progstr+": " << __LINE__ << errstr << "frequency must be positive" << endl; return 1; }
    if (frq>0.5) { cerr << progstr+": " << __LINE__ << errstr << "frequency must be <= 0.5" << endl; return 1; }

    //Get phs
    phs = (a_phs->count>0) ? a_phs->dval[0] : 0.0;


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
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        if (codee::sinewave_s(Y,o1.N(),(float)amp,(float)frq,(float)phs))
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
        if (codee::sinewave_d(Y,o1.N(),(double)amp,(double)frq,(double)phs))
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
        if (codee::sinewave_c(Y,N,(float)amp,(float)frq,(float)phs))
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
        if (codee::sinewave_z(Y,N,(double)amp,(double)frq,(double)phs))
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

