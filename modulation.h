#ifndef _MODULATION_H_
#define _MODULATION_H_

enum MODULATION_TYPE 
{
	MODULATION_SKIP = 0, 
	MODULATION_QAM4, 
	MODULATION_QAM16, 
	MODULATION_QAM64, 
	MODULATION_QAM256
};



// Only one definition from below must be defined
#define DET_PERM_V1
//#define DET_PERM_V2
//#define DET_PERM_V3

typedef struct
{
    int val;
    int ind;
} PERM_FOR_ORDER;

typedef struct
{
    int **HC;       // base matrix
    int b;          // number of rows
    int c;          // number of columns
    int M;
    int QAM;
    int mlog;
    int halfmlog;
    int Ncode;
    int Rcode;
    int n0;
    int c0;
    int *cw;
    int *V;
    int *p;
    int *q;
    int *s;
    int *perm;      // direct permutation
	int *perm_short;
	int *invperm_short;
    PERM_FOR_ORDER *temp;
    int *invperm;   // inverse permutation
    double *tmpbuf;
    int perm_mode;
    int extra_bits;
    int block_size;
    int step_size;
    int Nblocks;
	int ShortBlockSize;
	double *buffer;
	double *perm_codeword;
} PERMSTATE;


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
	double *dx;
} QAM_MODULATOR_STATE;

typedef struct
{
    
    double  T;
    double  sigma;
    short   Q;      //QAM
    int     n;      // code length
    int     m;      // log2( QAM )
    int     ns;     // number of signals
    int     DemodOutType;   
    double **x;
    //double *L;
    //double *P1;
} QAM_DEMODULATOR_STATE;

#define NOT_FOUND   -1



PERMSTATE* Permutations_Open( int b, int c, int M, int QAM, int halfmlog, int perm_mode, int block_size, int step_size );
void Permutations_Close( PERMSTATE *st );
void Permutation( PERMSTATE* state, int direction, double pinp[], double pout[] );
void Permutation_Init(  PERMSTATE* state, short **hc );

QAM_MODULATOR_STATE* QAM_modulator_open( int Q, int L, int m );
void QAM_modulator_close(QAM_MODULATOR_STATE* st);
void QAM_modulator( QAM_MODULATOR_STATE *qam_mod_state, double *in, double *out  );

QAM_DEMODULATOR_STATE* QAM_demodulator_open( double T, double sigma, short Q, int n, int m, int ns, int out_type );
void QAM_demodulator_close(QAM_DEMODULATOR_STATE* st );
void Demodulate( QAM_DEMODULATOR_STATE* st, double x[], double pRes[] );

#endif //_MODULATION_H_