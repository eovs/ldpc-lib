#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "modulation.h"

#ifndef SKIP_MEX
#include <mex.h>
#endif  //SKIP_MEX



#if 0
typedef struct
{
    int   Lorg; 
    int   Lfact;
    short Q;
    short m;
    short ns;
    short *x;
    short **y;
    short *z1;
    short *z2;
    short p[4];
} QAM_MODULATOR_STATE;

#define NOT_FOUND   -1
#endif  //0

static int Qav[4] = { 4, 16, 64, 256 };
#ifndef SKIP_MEX
static int logQav[4] = { 2, 4, 6, 8 };
#endif

int findQ( int q )
{
    for( int i = 0; i < 4; i++ )
    {
        if( q == Qav[i] )
            return i;
    }
    return NOT_FOUND;
                
}

void unpackRow_double2short(double m[], int width, short result[] )
{
    int j;
	for( j = 0; j < width; ++j)
	{
		result[j] = (short)m[j];
	}
}

void packMatrix_short2double( short *m[], int nrows, int ncols, double result[])
{
	for ( int i = 0; i < nrows; ++i )
	{
		for ( int j = 0; j < ncols; ++j)
		{
			result[j * nrows + i] = m[i][j];
//			mexPrintf("i = %d j = %d m[i][j] = %d\n", i, j, m[i][j] );
        }
	}
}


QAM_MODULATOR_STATE* QAM_modulator_open( int Q, int L, int m )
{
    QAM_MODULATOR_STATE* st;
    int i,j;
    
    st = (QAM_MODULATOR_STATE*)calloc(sizeof(QAM_MODULATOR_STATE),1);
    if( !st ) return NULL;
    
    st->Lorg = L;
    st->Q = Q;
    st->m = m;
    
    int r = L % m;
    if( r )
        st->Lfact = L+m-r;
    else
        st->Lfact = L;
    st->ns = st->Lfact / m;
    
    st->x = (short*)calloc( st->Lfact, sizeof(st->x[0]) );
    if( !st->x ) return NULL;
    
    st->y = (short**)calloc(2,sizeof( short* ));
    if( !st->y ) return NULL;
    st->y[0] = (short*)calloc(st->ns, sizeof(st->y[0][0]) );
    st->y[1] = (short*)calloc(st->ns, sizeof(st->y[1][0]) );
    if( !st->y[0] || !st->y[1] ) return NULL;
    
    st->z1 = (short*)calloc( st->ns, sizeof(st->z1[0]));
    if( !st->z1 ) return NULL;
    st->z2 = (short*)calloc( st->ns, sizeof(st->z2[0]));
    if( !st->z2 ) return NULL;
    
    for( i = m/2-1, j = 0; i>= 0; i--, j++ )
        st->p[j] = 1<<i;
    
    st->dx = (double*)calloc( st->Lfact, sizeof(st->dx[0]) );
    return st;
    
}

void QAM_modulator_close(QAM_MODULATOR_STATE* st)
{
    st->Lorg = 0;
    st->Q = 0;
    st->Lfact = 0;
    free( st->x );
    free( st->y[0] );
    free( st->y[1] );
    free( st->y );
    free( st->z1 );
    free( st->z2 );
	free( st->dx );
    
    free( st );
}


static short gray[]={0, 1, 3, 2, 7, 6, 4, 5, 15, 14, 12, 13,  8,  9, 11, 10};  //% anti-gray
static short s[] = {0,1,3,7, 15};
static void GrayPAM( short x[], int size, short y[], int order )
{
    //a=gray(x+1); % gray-coding
    //s=[1,3,7, 15];
    //y=2*a-s(order);
    for( int i = 0; i < size; i++ )
    {
        y[i] = 2*gray[x[i]]-s[order];
        //mexPrintf("i = %d x[i] = %d gray[x[i]] = %d s[order] = %d\n",i,x[i],gray[x[i]],s[order]);
    }
}


void QAM_modulator( QAM_MODULATOR_STATE *qam_mod_state, double *in, double *out  )
{
   int n, i, j, l;
   int m = qam_mod_state->m;
   
    memset( qam_mod_state->y[0], 0, qam_mod_state->ns * sizeof(qam_mod_state->y[0][0]) );
    memset( qam_mod_state->y[1], 0, qam_mod_state->ns * sizeof(qam_mod_state->y[1][0]) );
    
    //z1=p*z(1:m/2,:);
    //z2=p*z(m/2+1:m,:);
    int mdiv2 = qam_mod_state->m/2;
    for( n = 0, j = 0; n < qam_mod_state->Lfact; n += qam_mod_state->m, j++ )
    {
        int sum = 0;
        for( i = 0; i < mdiv2; i++ )
            sum += (int)(qam_mod_state->p[i] * in[n+i]);
        #ifndef SKIP_MEX
        if( sum < 0 || sum >= (1<<mdiv2) )
            mexErrMsgTxt("input out of range in QAM modulator");
        #endif
        qam_mod_state->z1[j] = sum;
        
        sum = 0;
        for( i = mdiv2, l = 0; i < m; i++, l++ )
            sum += (int)(qam_mod_state->p[l] * in[n+i]);
        #ifndef SKIP_MEX
        if( sum < 0 || sum >= (1<<mdiv2) )
            mexErrMsgTxt("input out of range in QAM modulator");
        #endif
        qam_mod_state->z2[j] = sum;
        
    }
//    mexPrintf("z1[0] = %d z1[1] = %d z2[0] = %d z2[1] = %d\n", qam_mod_state->z1[0], qam_mod_state->z1[1],qam_mod_state->z2[0],qam_mod_state->z2[1]);
//    mexPrintf("p[0] = %d p[1] = %d p[2] = %d p[3] = %d\n", qam_mod_state->p[0], qam_mod_state->p[1],qam_mod_state->p[2],qam_mod_state->p[3]);
    //y(1,:)=GrayPAM(z1,(m/2));
    //y(2,:)=GrayPAM(z2,(m/2));
    #ifndef SKIP_MEX
    if( mdiv2 > 4 )  mexErrMsgTxt("too high order in QAM modulator");
    #endif
    GrayPAM( qam_mod_state->z1, qam_mod_state->ns, qam_mod_state->y[0], mdiv2 );
    GrayPAM( qam_mod_state->z2, qam_mod_state->ns, qam_mod_state->y[1], mdiv2 );
    
    //pOut[0]=mxCreateDoubleMatrix(2,qam_mod_state->ns,mxREAL);
    //double *p = mxGetPr(pOut[0]);
    //packMatrix_short2double( qam_mod_state->y, 2, qam_mod_state->ns, p );
    for( i = j = 0; i < qam_mod_state->ns; i++, j+= 2 )
    {
        out[j]   = qam_mod_state->y[0][i];
        out[j+1] = qam_mod_state->y[1][i];
    }
        
    
}

#if 0
function y=GrayPAM(x, order)
% transforms x from range 0...2^order-1 to PAM-code

if order > 4,  error('too high order'), end;
if any(x<0) || any(x>=2^order), error('input out of range'), end;
%gray=[0 1 3 2 6 7 5 4 12 13 15 14 10 11  9  8]; % gray 
gray=[0 1 3 2 7 6 4 5 15 14 12 13  8  9 11 10];  % anti-gray

%  000 001 011 010 110 111 101 100
%  0000 0001 0011 0010 0110 0111 0101 0100
%  1100 1101 1111 1110 1010 1011 1001 1000
a=gray(x+1); % gray-coding
s=[1,3,7, 15];
y=2*a-s(order);
#endif  //0


#ifndef SKIP_MEX
//function y = GrayQAMmodulator( x, Q)
void mexFunction(int nOut, mxArray *pOut[], int nInp, const mxArray *pInp[])
{
    static QAM_MODULATOR_STATE *qam_mod_state = NULL;
    static int Q = 0;
    static int L = 0;
    int i,j,n,l;
    int k, r;
    int m;
    int indQ;
   // bool newmem = false;
    //int p[4];
    int ns; // number of signals
    
    if( nInp != 2 )
    {
         mexErrMsgTxt("Only 2 input argument allowed for GrayQAMmodulatorC");
    }
    if( nOut != 1 )
    {
        mexErrMsgTxt("Only one output argument allowed for GrayQAMmodulatorC");
    }
    
    k = (int)mxGetM(pInp[0]);
    L = (int)mxGetN(pInp[0]);
    if( k != 1 )
    {
        mexErrMsgTxt("Only 1xL input vector allowed");
    }
    
    Q = (int)mxGetPr(pInp[1])[0];
    indQ = findQ( Q );
    if( indQ == NOT_FOUND )
    {
        mexErrMsgTxt("Only QAM-4, 16, 64, 256 are allowed");
    }
    m = logQav[indQ];
    
    
    
    if( qam_mod_state == NULL  )
    {   // need to be opened
        qam_mod_state = QAM_modulator_open( Q, L, m );
        if( !qam_mod_state )
             mexErrMsgTxt("Alloacation error in QAM_modulator_open" );
    }
    else
    {
        if( (Q != qam_mod_state->Q ) || (L != qam_mod_state->Lorg ) )
        {
            QAM_modulator_close( qam_mod_state );
            qam_mod_state = QAM_modulator_open( Q, L, m );  // reopen for another parameters
            if( !qam_mod_state )
                 mexErrMsgTxt("Alloacation error in QAM_modulator_open" );
        }
    }

    unpackRow_double2short(mxGetPr(pInp[0]), L, qam_mod_state->x );
    
    #if 0
    //y=zeros(2,ns);
    memset( qam_mod_state->y[0], 0, qam_mod_state->ns * sizeof(qam_mod_state->y[0][0]) );
    memset( qam_mod_state->y[1], 0, qam_mod_state->ns * sizeof(qam_mod_state->y[1][0]) );
    
    //z1=p*z(1:m/2,:);
    //z2=p*z(m/2+1:m,:);
    int mdiv2 = m/2;
    for( n = 0, j = 0; n < qam_mod_state->Lfact; n += qam_mod_state->m, j++ )
    {
        int sum = 0;
        for( i = 0; i < mdiv2; i++ )
            sum += qam_mod_state->p[i] * qam_mod_state->x[n+i];
        if( sum < 0 || sum >= (1<<mdiv2) )
            mexErrMsgTxt("input out of range in QAM modulator");
        qam_mod_state->z1[j] = sum;
        
        sum = 0;
        for( i = mdiv2, l = 0; i < m; i++, l++ )
            sum += qam_mod_state->p[l] *qam_mod_state->x[n+i];
        if( sum < 0 || sum >= (1<<mdiv2) )
            mexErrMsgTxt("input out of range in QAM modulator");
        qam_mod_state->z2[j] = sum;
        
    }
//    mexPrintf("z1[0] = %d z1[1] = %d z2[0] = %d z2[1] = %d\n", qam_mod_state->z1[0], qam_mod_state->z1[1],qam_mod_state->z2[0],qam_mod_state->z2[1]);
//    mexPrintf("p[0] = %d p[1] = %d p[2] = %d p[3] = %d\n", qam_mod_state->p[0], qam_mod_state->p[1],qam_mod_state->p[2],qam_mod_state->p[3]);
    //y(1,:)=GrayPAM(z1,(m/2));
    //y(2,:)=GrayPAM(z2,(m/2));
    if( mdiv2 > 4 )  mexErrMsgTxt("too high order in QAM modulator");
    GrayPAM( qam_mod_state->z1, qam_mod_state->ns, qam_mod_state->y[0], mdiv2 );
    GrayPAM( qam_mod_state->z2, qam_mod_state->ns, qam_mod_state->y[1], mdiv2 );
    
    pOut[0]=mxCreateDoubleMatrix(2,qam_mod_state->ns,mxREAL);
   
    double *p = mxGetPr(pOut[0]);
    packMatrix_short2double( qam_mod_state->y, 2, qam_mod_state->ns, p );
    #endif  //0

    pOut[0]=mxCreateDoubleMatrix(2,qam_mod_state->ns,mxREAL);
    double *p = mxGetPr(pOut[0]);

	for( i = 0; i < qam_mod_state->Lfact; i++ )
		qam_mod_state->dx[i] = qam_mod_state->x[i];

    QAM_modulator( qam_mod_state, qam_mod_state->dx, p );
    

}
#endif	//SKIP_MEX