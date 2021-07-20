//Includes
#include "convert_freqs.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u};
const size_t I = 1u, O = 1u;
string iscale, oscale;

//Description
string descr;
descr += "Converts frequencies from one scale to another.\n";
descr += "Input X has frequencies on the original frequency scale.\n";
descr += "\n";
descr += "Output Y has frequencies on the new frequency scale.\n";
descr += "Y has the same size, data type and file format as X.\n";
descr += "\n";
descr += "Use -f (--from) to give a string for the input frequency scale.\n";
descr += "Use -t (--to) to give a string for the output frequency scale.\n";
descr += "Both input strings are lower-cased before use.\n";
descr += "\n";
descr += "The available frequency scales are: \n";
descr += "'bark', 'cochlea', 'erb', 'hz', 'mel', 'midi', 'octave', 'piano', \n";
descr += "'ihcs' (inner hair cells), and 'sgcs' (spiral ganglion cells).\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ convert_freqs -f'hz' -t'mel' X -o Y \n";
descr += "$ convert_freqs -f'hz' -t'sgcs' X > Y \n";
descr += "$ cat X | convert_freqs -f'bark' -t'midi' > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_str   *a_is = arg_strn("f","from","<str>",0,1,"input frequency scale [default='hz']");
struct arg_str   *a_os = arg_strn("t","to","<str>",0,1,"output frequency scale [default='hz']");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get iscale
if (a_is->count==0) { iscale = "hz"; }
else
{
	try { iscale = string(a_is->sval[0]); }
	catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem getting string for input frequency scale" << endl; return 1; }
}
for (string::size_type c=0u; c<iscale.size(); ++c) { iscale[c] = char(tolower(iscale[c])); }

//Get oscale
if (a_os->count==0) { oscale = "hz"; }
else
{
	try { oscale = string(a_os->sval[0]); }
	catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem getting string for output frequency scale" << endl; return 1; }
}
for (string::size_type c=0u; c<oscale.size(); ++c) { oscale[c] = char(tolower(oscale[c])); }

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = i1.R; o1.C = i1.C; o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1u)
{
    float *X;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::convert_freqs_s(X,i1.N(),iscale.c_str(),oscale.c_str())) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(X),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X;
}

//Finish
