//Includes
#include "fft.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2,101,102};
const size_t I = 1, O = 1;
size_t dim, nfft;

//Description
string descr;
descr += "Does 1D FFT along rows or cols of RxC input matrix X.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension along which to transform.\n";
descr += "Use -d0 to operate along cols, and -d1 to operate along rows.\n";
descr += "The default is 0 (along cols), unless X is a row vector.\n";
descr += "\n";
descr += "Use -n (--nfft) to specify transform length [default is R or C].\n";
descr += "X is zero-padded as necessary to match nfft.\n";
descr += "\n";
descr += "The output (Y) is complex-valued with size: \n";
descr += "d=0, X real   :   nfrqs x C \n";
descr += "d=1, X real   :   R x nfrqs \n";
descr += "d=0, X complex:   nfft x C \n";
descr += "d=1, X complex:   R x nfft \n";
descr += "where nfrqs = floor(nfft/2)+1 = num nonnegative FFT frequencies.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ fft -n256 X -o Y \n";
descr += "$ fft -n256 -d1 X > Y \n";
descr += "$ cat X | fft -n256 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to take FFT [default=0]");
struct arg_int    *a_n = arg_intn("n","nfft","<uint>",0,1,"transform length [default is R or C]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = 0; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim!=0 && dim!=1) { cerr << progstr+": " << __LINE__ << errstr << "dim must be 0 or 1" << endl; return 1; }

//Get nfft
if (a_n->count==0) { nfft = (dim==0) ? i1.R : i1.C; }
else if (a_n->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be positive" << endl; return 1; }
else { nfft = a_n->ival[0]; }

//Checks
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) must be 1D or 2D" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
if (dim==0 && nfft<i1.R) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be > nrows of X" << endl; return 1; }
if (dim==1 && nfft<i1.C) { cerr << progstr+": " << __LINE__ << errstr << "nfft must be > ncols of X" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.isreal() ? i1.T+100 : i1.T;
o1.R = (dim==1) ? i1.R : i1.isreal() ? uint32_t(nfft)/2+1 : uint32_t(nfft);
o1.C = (dim==0) ? i1.C : i1.isreal() ? uint32_t(nfft)/2+1 : uint32_t(nfft);
o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1)
{
    float *X, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { Y = new float[2*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    if (codee::fft_s(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,nfft)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}
else if (i1.T==101)
{
    float *X, *Y;
    try { X = new float[2*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { Y = new float[2*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    if (codee::fft_c(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,nfft)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish
