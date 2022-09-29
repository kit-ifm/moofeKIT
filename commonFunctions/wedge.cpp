/*
 * function C = fastWedge(A,B)
 *
 * mex -v fastWedge.cpp
 *
 */

/* $Revision: 0.1 R.Pfefferkorn 22.7.2018 $ */
/* $Revision: 0.2 R.Pfefferkorn 28.8.2018 $ - added complex number support */
#include <mex.h>
#include <matrix.h>
#include <iostream>

using namespace std;

/* Input Arguments */
#define	A_IN prhs[0]
#define	B_IN prhs[1]

/* Output Arguments */
#define	C_OUT plhs[0]

/* Indices */
enum indMat:size_t {XX=0,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ};
enum indVec:size_t {X=0,Y,Z};


// /* ###################################################################### */
// /* Output to MATLAB console */
// /* ###################################################################### */
// class mystream : public std::streambuf
// {
// protected:
// virtual std::streamsize xsputn(const char *s, std::streamsize n) { mexPrintf("%.*s", n, s); return n; }
// virtual int overflow(int c=EOF) { if (c != EOF) { mexPrintf("%.1s", &c); } return 1; }
// };
// class scoped_redirect_cout
// {
// public:
// 	scoped_redirect_cout() { old_buf = std::cout.rdbuf(); std::cout.rdbuf(&mout); }
// 	~scoped_redirect_cout() { std::cout.rdbuf(old_buf); }
// private:
// 	mystream mout;
// 	std::streambuf *old_buf;
// };
// static scoped_redirect_cout mycout_redirect;

/* ###################################################################### */
/* Computational functions COMPLEX input */
/* ###################################################################### */
static void complexMult(double *zr,double *zi, const double ar, const double ai, const double br, const double bi){
    *(zr) += ar*br-ai*bi;
    *(zi) += ar*bi+ai*br;
    return;
}
static void complexCrossHelp(double *Cr, double *Ci, const double *Ar, const double *Ai, const double *Br, const double *Bi, 
size_t a, size_t x, size_t y, size_t u, size_t v){
    Cr[a] = 0;
    Ci[a] = 0;
    complexMult(&Cr[a],&Ci[a], Ar[x], Ai[x], Br[u], Bi[u]);
    complexMult(&Cr[a],&Ci[a],-Ar[y],-Ai[y], Br[v], Bi[v]);
    complexMult(&Cr[a],&Ci[a], Ar[u], Ai[u], Br[x], Bi[x]);
    complexMult(&Cr[a],&Ci[a],-Ar[v],-Ai[v], Br[y], Bi[y]);
    return;
}
static void calcCrossMatMatComplex(mxArray *Cmx, const mxArray *Amx, const mxArray *Bmx){
    //get pointers to complex arrays
    double *Ar,*Br,*Cr, *Ai,*Bi,*Ci;
    Ar = mxGetPr(Amx);
    Br = mxGetPr(Bmx);
    Cr = mxGetPr(Cmx);
    Ai = mxGetPi(Amx);
    Bi = mxGetPi(Bmx);
    Ci = mxGetPi(Cmx);
    
    //create zero arrays for imaginary part if pointers are empty
    bool delAi=false, delBi=false;
    if (Ai==NULL){
        Ai = new double[9];
        for (int i;i<9;i++){Ai[i]=0;}
        delAi = true;
    }
    if (Bi==NULL){
        Bi = new double[9];
        for (int i;i<9;i++){Bi[i]=0;}
        delBi = true;
    }
    
    double *A,*B,*C;
    A = mxGetPr(Amx);
    B = mxGetPr(Bmx);
    C = mxGetPr(Cmx);
    
    //bonet,gil,ortigosa (2015)
    complexCrossHelp(Cr,Ci,Ar,Ai,Br,Bi, XX, YY,YZ,ZZ,ZY);
    complexCrossHelp(Cr,Ci,Ar,Ai,Br,Bi, XY, YZ,YX,ZX,ZZ);
    complexCrossHelp(Cr,Ci,Ar,Ai,Br,Bi, XZ, YX,YY,ZY,ZX);
    
    complexCrossHelp(Cr,Ci,Ar,Ai,Br,Bi, YX, XZ,XY,ZY,ZZ);
    complexCrossHelp(Cr,Ci,Ar,Ai,Br,Bi, YY, ZZ,ZX,XX,XZ);
    complexCrossHelp(Cr,Ci,Ar,Ai,Br,Bi, YZ, ZX,ZY,XY,XX);
    
    complexCrossHelp(Cr,Ci,Ar,Ai,Br,Bi, ZX, XY,XZ,YZ,YY);
    complexCrossHelp(Cr,Ci,Ar,Ai,Br,Bi, ZY, XZ,XX,YX,YZ);
    complexCrossHelp(Cr,Ci,Ar,Ai,Br,Bi, ZZ, XX,XY,YY,YX);
    
    //delete memory
    if (delAi) delete [] Ai;
    if (delBi) delete [] Bi;
    
    return;
}

/* ###################################################################### */
/* Computational functions REAL input */
/* ###################################################################### */
static void calcCrossMatMat(double *C, const double *A, const double *B){
    //bonet,gil,ortigosa (2015)
    C[XX] = A[YY]*B[ZZ] - A[YZ]*B[ZY] + A[ZZ]*B[YY] - A[ZY]*B[YZ];
    C[XY] = A[YZ]*B[ZX] - A[YX]*B[ZZ] + A[ZX]*B[YZ] - A[ZZ]*B[YX];
    C[XZ] = A[YX]*B[ZY] - A[YY]*B[ZX] + A[ZY]*B[YX] - A[ZX]*B[YY];
    
    C[YX] = A[XZ]*B[ZY] - A[XY]*B[ZZ] + A[ZY]*B[XZ] - A[ZZ]*B[XY];
    C[YY] = A[ZZ]*B[XX] - A[ZX]*B[XZ] + A[XX]*B[ZZ] - A[XZ]*B[ZX];
    C[YZ] = A[ZX]*B[XY] - A[ZY]*B[XX] + A[XY]*B[ZX] - A[XX]*B[ZY];
    
    C[ZX] = A[XY]*B[YZ] - A[XZ]*B[YY] + A[YZ]*B[XY] - A[YY]*B[XZ];
    C[ZY] = A[XZ]*B[YX] - A[XX]*B[YZ] + A[YX]*B[XZ] - A[YZ]*B[XX];
    C[ZZ] = A[XX]*B[YY] - A[XY]*B[YX] + A[YY]*B[XX] - A[YX]*B[XY];
    
    return;
}
static void calcCrossMatMatSame(double *C, const double *A){
    //wedge function from esra (speed up for the same matrix twice)    
    C[XX] = 2*(A[YY]*A[ZZ] - A[YZ]*A[ZY]);
    C[XY] = 2*(A[YZ]*A[ZX] - A[YX]*A[ZZ]);
    C[XZ] = 2*(A[YX]*A[ZY] - A[YY]*A[ZX]);
    
    C[YX] = 2*(A[XZ]*A[ZY] - A[XY]*A[ZZ]);
    C[YY] = 2*(A[ZZ]*A[XX] - A[ZX]*A[XZ]);
    C[YZ] = 2*(A[ZX]*A[XY] - A[ZY]*A[XX]);
    
    C[ZX] = 2*(A[XY]*A[YZ] - A[XZ]*A[YY]);
    C[ZY] = 2*(A[XZ]*A[YX] - A[XX]*A[YZ]);
    C[ZZ] = 2*(A[XX]*A[YY] - A[XY]*A[YX]);
    
    return;
}
static void calcCrossVecMat(double *C, const double *v, const double *A){
    //bonet,gil,ortigosa (2015)
    C[XX] = v[Y]*A[ZX]-v[Z]*A[YX];
    C[YX] = v[Z]*A[XX]-v[X]*A[ZX];
    C[ZX] = v[X]*A[YX]-v[Y]*A[XX];
    
    C[XY] = v[Y]*A[ZY]-v[Z]*A[YY];
    C[YY] = v[Z]*A[XY]-v[X]*A[ZY];
    C[ZY] = v[X]*A[YY]-v[Y]*A[XY];
    
    C[XZ] = v[Y]*A[ZZ]-v[Z]*A[YZ];
    C[YZ] = v[Z]*A[XZ]-v[X]*A[ZZ];
    C[ZZ] = v[X]*A[YZ]-v[Y]*A[XZ];
    
    return;
}
static void calcCrossMatVec(double *C, const double *A, const double *v){
    //bonet,gil,ortigosa (2015)
    C[XX] = v[Z]*A[XY]-v[Y]*A[XZ];
    C[YX] = v[Z]*A[YY]-v[Y]*A[YZ];
    C[ZX] = v[Z]*A[ZY]-v[Y]*A[ZZ];
    
    C[XY] = v[X]*A[XZ]-v[Z]*A[XX];
    C[YY] = v[X]*A[YZ]-v[Z]*A[YX];
    C[ZY] = v[X]*A[ZZ]-v[Z]*A[ZX];
    
    C[XZ] = v[Y]*A[XX]-v[X]*A[XY];
    C[YZ] = v[Y]*A[YX]-v[X]*A[YY];
    C[ZZ] = v[Y]*A[ZX]-v[X]*A[ZY];
    
    return;
}

/* ###################################################################### */
/* Interface function */
/* ###################################################################### */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
    double *A,*B;
    size_t Am,An,Bm,Bn;
    
    /* Check for proper number of arguments */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:wedge:invalidNumInputs","Two input arguments required.");
    } else if (nlhs > 1) {
        mexErrMsgIdAndTxt("MATLAB:wedge:maxlhs","Too many output arguments.");
    }
    
    /* Check data type of input argument */
    if (!(mxIsDouble(A_IN)) || !(mxIsDouble(B_IN))) {
        mexErrMsgIdAndTxt( "MATLAB:wedge:invalidInputType","Input arrays must be of type double.");
    }
    /* Check number of dimensions of input */
    if (mxGetNumberOfDimensions(A_IN)!=2 || mxGetNumberOfDimensions(B_IN)!=2) {
        mexErrMsgIdAndTxt("MATLAB:wedge:invalidInputType","Input arrays must have two dimensions, e.g. 3x3.");
    } 
    
    /* dimensions of input parameters */
    A = mxGetPr(A_IN);
    B = mxGetPr(B_IN);
    Am = mxGetM(A_IN);
    An = mxGetN(A_IN);
    Bm = mxGetM(A_IN);
    Bn = mxGetN(B_IN);
    
    bool isvecA, isvecB, ismatA, ismatB;
    ismatA = Am==3 && An==3;
    ismatB = Bm==3 && Bn==3;
    isvecA = (Am==3 && An==1) || (Am==1 && An==3);
    isvecB = (Bm==3 && Bn==1) || (Bm==1 && Bn==3);
    
    bool isComplex;
    isComplex = (mxIsComplex(A_IN) || mxIsComplex(B_IN));
                
    //REAL INPUT
    if (!isComplex) {
        /* input pointers*/
        A = mxGetPr(A_IN);
        B = mxGetPr(B_IN);
        /* Output */
        C_OUT = mxCreateDoubleMatrix(3, 3, mxREAL);
        
        /*computations */
        if (ismatA && ismatB) {
            if (A==B) {
                calcCrossMatMatSame(mxGetPr(C_OUT),A);}
            else {
                calcCrossMatMat(mxGetPr(C_OUT),A,B);}}
        else if (isvecA && ismatB) {
            calcCrossVecMat(mxGetPr(C_OUT),A,B);}
        else if (ismatA && isvecB) {
            calcCrossMatVec(mxGetPr(C_OUT),A,B);}
        else {
            mexErrMsgIdAndTxt("MATLAB:fastWedge:invalid","Check dimension of both inputs.");}
    }//COMPLEX INPUT
    else {
        /* Output */
        C_OUT = mxCreateDoubleMatrix(3, 3, mxCOMPLEX);
        if (ismatA && ismatB) {
           calcCrossMatMatComplex(C_OUT,A_IN,B_IN);}
        else {
            mexErrMsgIdAndTxt("MATLAB:fastWedge:invalid","Check dimension of both inputs.");}
    }
    
    //end function
    return;
}