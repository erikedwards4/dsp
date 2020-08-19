//Includes
#include "fir.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2,101,102};
const size_t I = 2, O = 1;
size_t dim, N, T, L;

//Description
string descr;
descr += "Neural soma stage.\n";
descr += "Does FIR filtering of each row or col of X.\n";
descr += "using FIR filter coefficients in B.\n";
descr += "\n";
descr += "X has size NxT or TxN, where N is the number of neurons\n";
descr += "and T is the number of observations (e.g. time points).\n";
descr += "\n";
descr += "Use -d (--dim) to specify the dimension (axis) of length N [default=0].\n";
descr += "\n";
descr += "For dim=0, Y[n,t] = B[n,0]*X[n,t-L+1] + B[n,1]*X[n,t-L+2] + ... + B[n,L-1]*X[n,t]\n";
descr += "with sizes X:  N x T \n";
descr += "           B:  N x L \n";
descr += "           Y:  N x T \n";
descr += "\n";
descr += "For dim=1, Y[t,n] = B[0,n]*X[t-L+1,n] + B[1,n]*X[t-L+2,n] + ... + B[L-1,n]*X[t,n]\n";
descr += "with sizes X:  T x N \n";
descr += "           B:  L x N \n";
descr += "           Y:  T x N \n";
descr += "\n";
descr += "Note that B is a matrix (1 row or col per neuron).\n";
descr += "Note that each row or col of B is flipped before convolving (usual convention),\n";
descr += "so each row or col is in reverse chronological order (past towards right or bottom).\n";
descr += "\n";
descr += "If L=1, this is a simple scaling of each row or col of X by B[n].\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ fir X B -o Y \n";
descr += "$ fir X B > Y \n";
descr += "$ cat X | fir - B > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X,B)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension (0 or 1) [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = (i1.C==1u) ? 1 : 0; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>1) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1}" << endl; return 1; }

//Checks
if (i1.T!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) found to be empty" << endl; return 1; }
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) must be a matrix" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (B) found to be empty" << endl; return 1; }
if (!i2.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (B) must be a matrix" << endl; return 1; }
if (i2.N()!=1u)
{
    if (dim==0 && i2.R!=i1.R) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) and 2 (B) must have same num rows for dim=0" << endl; return 1; }
    if (dim==1 && i2.C!=i1.C) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) and 2 (B) must have same num cols for dim=1" << endl; return 1; }
}

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = i1.R; o1.C = i1.C; o1.S = i1.S; o1.H = i1.H;

//Other prep
N = (dim==0) ? int(o1.R) : int(o1.C);
T = (dim==0) ? int(o1.C) : int(o1.R);
L = (dim==0) ? int(i2.R) : int(i2.C);
if (T<2) { cerr << progstr+": " << __LINE__ << errstr << "num time points must be > 1" << endl; return 1; }

//Process
if (i1.T==1)
{
    float *X, *B, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { B = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (B)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(B),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (B)" << endl; return 1; }
    if (openn::fir_s(Y,X,B,N,T,L,dim,i1.iscolmajor()))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; } 
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] B; delete[] Y;
}
else if (i1.T==101)
{
    float *X, *B, *Y;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { B = new float[2u*i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (B)" << endl; return 1; }
    try { Y = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(B),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (B)" << endl; return 1; }
    if (openn::fir_c(Y,X,B,N,T,L,dim,i1.iscolmajor()))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; } 
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] B; delete[] Y;
}

//Finish
