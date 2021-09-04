//Includes
#include "tkeo_smooth.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 1u, O = 1u;
size_t dim, Lx, Ly, str, dil;
int cs0;

//Description
string descr;
descr += "Does smooth TKEO along dim for each vector in X.\n";
descr += "Teager-Kaiser Energy Operator (TKEO) according to de Matos [2018],\n";
descr += "which is computed using outputs of smooth_diff and smooth_diffdiff:\n";
descr += "y[n] = x'[n]*x'[n] - x[n]*x''[n]\n";
descr += "\n";
descr += "For complex X, the TKEO is real-valued and is the sum of the TKEO\n";
descr += "for real and imaginary parts [Hamila et al. 1999; de Matos 2018].\n";
descr += "\n";
descr += "For real X, the exact method of de Matos [2018] would be:\n";
descr += "analytic_sig X | tkeo_smooth > Y\n";
descr += "where Y is the \"complex TK energy\" (but is real-valued).\n";
descr += "A further Hilbert transform of the complex TK energy gives the TKV.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension along which to operate.\n";
descr += "Default is 0 (along cols), unless X is a row vector.\n";
descr += "\n";
descr += "Use -s (--stride) to give the stride (step-size) in samples [default=1].\n";
descr += "\n";
descr += "Use -i (--dilation) to give the dilation factor [default=1].\n";
descr += "\n";
descr += "Use -c (--cs0) to give the center-samp of the initial frame [default=0].\n";
descr += "\n";
descr += "Use -l (--Ly) to give the length of output vecs in Y [default=Lx/str]\n";
descr += "\n";
descr += "Y has the same size as X1, but the vecs along dim have length Ly.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ tkeo_smooth X -o Y \n";
descr += "$ tkeo_smooth -d1 X -o Y \n";
descr += "$ cat X | tkeo_smooth -d1 -s2 -i2 -p0 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to operate [default=0]");
struct arg_int  *a_cs0 = arg_intn("c","cs0","<int>",0,1,"center samp of initial frame [default=0]");
struct arg_int  *a_str = arg_intn("s","stride","<uint>",0,1,"stride (step size) in samps [default=1]");
struct arg_int  *a_dil = arg_intn("i","dilation","<uint>",0,1,"dilation factor [default=1]");
struct arg_int   *a_ly = arg_intn("l","Ly","<uint>",0,1,"length of vecs in output Y [default=Lx/str]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = i1.isrowvec() ? 1u : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Get cs0
if (a_cs0->count==0) { cs0 = 0; }
else { cs0 = a_cs0->ival[0]; }

//Get stride
if (a_str->count==0) { str = 1u; }
else if (a_str->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "stride must be positive" << endl; return 1; }
else { str = size_t(a_str->ival[0]); }

//Get dil
if (a_dil->count==0) { dil = 1u; }
else if (a_dil->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "dilation must be positive" << endl; return 1; }
else { dil = size_t(a_dil->ival[0]); }

//Get Ly
Lx = (dim==0u) ? i1.R : (dim==1u) ? i1.C : (dim==2u) ? i1.S : i1.H;
if (a_ly->count==0) { Ly = Lx/str; }
else if (a_ly->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "Ly (length of vecs in Y) must be positive" << endl; return 1; }
else { Ly = size_t(a_ly->ival[0]); }

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }

//Set output header info
o1.F = i1.F;
o1.T = i1.isreal() ? i1.T : i1.T-100u;
o1.R = (dim==0u) ? Ly : i1.R;
o1.C = (dim==1u) ? Ly : i1.C;
o1.S = (dim==2u) ? Ly : i1.S;
o1.H = (dim==3u) ? Ly : i1.H;

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
    if (codee::tkeo_smooth_s(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,cs0,str,dil,Ly))
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
    if (codee::tkeo_smooth_c(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,cs0,str,dil,Ly))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish
