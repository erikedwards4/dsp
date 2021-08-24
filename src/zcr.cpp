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
#include "zcr.c"

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
    size_t dim, Lx, Lw;
    int g, causal;


    //Description
    string descr;
    descr += "Gets the zero crossing rate (ZCR) of X along dim.\n";
    descr += "This uses the usual definition (i.e., rectangular window of length Lw).\n";
    descr += "\n";
    descr += "Use -l (--Lw) to give the length of the moving-average window.\n";
    descr += "\n";
    descr += "Use -g (--going) to specify positive- or negative-going ZCs.\n";
    descr += "Use -g0 to detect positive- and negative-going ZCs [default].\n";
    descr += "Use -g1 to detect only positive-going ZCs.\n";
    descr += "Use -g-1 to detect only negative-going ZCs.\n";
    descr += "\n";
    descr += "Use -d (--dim) to give the dimension along which to operate.\n";
    descr += "The default is 0 (along cols), unless X is a row vector.\n";
    descr += "\n";
    descr += "Include -c (--causal) to give causal output [default=false].\n";
    descr += "The default (noncausal) uses a centered moving-average window.\n";
    descr += "The causal option averages only over the current and past samps.\n";
    descr += "\n";
    descr += "Y is real-valued and has the same size and file format as X.\n";
    descr += "Y is in units of ZCs/sample (so multiply by fs to get ZCs/sec).\n";
    descr += "Thus, each element of Y is in [0.0 1.0].\n";
    descr += "\n";
    descr += "Examples:\n";
    descr += "$ zcr -l100 X -o Y \n";
    descr += "$ zcr -d1 -l7 X > Y \n";
    descr += "$ cat X | zcr -l127 -g1 -c > Y \n";


    //Argtable
    int nerrs;
    struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
    struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to operate [default=0]");
    struct arg_int    *a_l = arg_intn("l","Lw","<uint>",0,1,"winlength over which to average ZCs [default=1]");
    struct arg_int    *a_g = arg_intn("g","going","<uint>",0,1,"if using positive- or negative-going ZCs [default=0]");
    struct arg_lit    *a_c = arg_litn("c","causal",0,1,"include to give causal output [default=false]");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_fi, a_d, a_l, a_g, a_c, a_fo, a_help, a_end};
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

    //Get g
    g = (a_g->count>0) ? a_g->ival[0] : 0;
    if (g!=0 && g!=1 && g!=-1) { cerr << progstr+": " << __LINE__ << errstr << "g must be in {-1,0,1}" << endl; return 1; }

    //Get Lw
    if (a_l->count==0) { Lw = 1u; }
    else if (a_l->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "Lw must be positive" << endl; return 1; }
    else { Lw = size_t(a_l->ival[0]); }

    //Get g
    g = (a_g->count>0) ? a_g->ival[0] : 0;
    if (g!=0 && g!=1 && g!=-1) { cerr << progstr+": " << __LINE__ << errstr << "g must be in {-1,0,1}" << endl; return 1; }

    //Get causal
    causal = (a_c->count>0);


    //Checks
    Lx = (dim==0u) ? i1.R : (dim==1u) ? i1.C : (dim==2u) ? i1.S : i1.H;
    if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
    if (Lx<2u) { cerr << progstr+": " << __LINE__ << errstr << "cannot work along a singleton dimension" << endl; return 1; }
    if (Lx<Lw) { cerr << progstr+": " << __LINE__ << errstr << "Lw (winlength) must be <= Lx (length of vecs in X)" << endl; return 1; }


    //Set output header info
    o1.F = i1.F;
    o1.T = i1.iscomplex() ? i1.T-100u : i1.T;
    o1.R = i1.R; o1.C = i1.C; o1.S = i1.S; o1.H = i1.H;


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
    if (i1.T==1u)
    {
        float *X, *Y;
        try { X = new float[i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
        try { Y = new float[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        if (codee::zcr_s(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,Lw,g,causal))
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
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
        try { Y = new double[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        if (codee::zcr_d(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,Lw,g,causal))
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
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
        try { Y = new float[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        if (codee::zcr_c(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,Lw,g,causal))
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
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
        try { Y = new double[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        if (codee::zcr_z(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,Lw,g,causal))
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

