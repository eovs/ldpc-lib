#include "graph.h"

//Globals for code analysis
int CNodes[MAX_EQ];     // current nodes
int f[MAX_EQ];          // auxilary storage used in find routines
int ACTIVE[MAX_EQ];     // auxilary 
int DEPTH[MAX_EQ];      // 
int PBranch[MAX_EQ];    // previous branch
int WT[MAX_EQ];         //,WC[MAX_EQ];

int d_prev;
int path[NMAX];


// For fast equation test
int DTABLE[DMAX][NDMAX];
int FLAGS[MAX_EQ],
	FLAGSS[MAX_EQ];
int T[CMAX][CMAX];

int H[CMAX][CMAX*WMAX];
int	HH[CMAX][CMAX*WMAX];
//int	HTB[RMAX][CMAX*WMAX];
int	HM[RMAX][CMAX]; //Mask 
int HD[CMAX][CMAX];


// int    HQ[CMAX][CMAX];
// Variables for equations
int	EQ[MAX_EQ][CMAX],
	ACE[MAX_EQ],
//Depth[MAX_EQ],
	CYCLES[MAX_EQ][2],
	CYCLESS[MAX_EQ][2],
//  Index[MAX_EQ],
	TREE[MAX_EQ][2],
	TREES[MAX_EQ][2];
	//tree[MAX_EQ][2];
int   RESULTSG[CMAX];
//    RESULTSA[NMAX][CMAX];
int Graph[NMAX][CMAX];
int Branches[NMAX][CMAX];
int Numbers[NMAX]; 
int Path[MAX_EQ][CMAX*WMAX]; // Path to node (part of equation)
int tree[CMAX*MAX_EQ][2];    // temporal tree 

// DECODING

int C[NMAX][WMAX];       // columns
int B[NMAX][WMAX],       // hard decisions
	V[NMAX][WMAX];       //rows
int VI[NMAX][WMAX];      //indices of rows in columns
int codeword[NMAX];
float received[NMAX],
	noise[NMAX];
//float  L[NMAX][NMAX];
double  Z[NMAX][WMAX];

// BASE MATRIX GENERATION
int ncombs[CMAX];
int combs3[300][3],
	combs4[12000][4],
	combs5[50000][5],
	combs6[150000][6],
	combs7[350000][7],
	combs8[750000][8],
	combs9[750000][9],
    combs10[75000][10],
	combs11[75000][11],
	combs12[75000][12],
    combs13[75000][13],
    combs14[75000][14],
    combs15[75000][15],
    combs16[75000][16];
int col_set[200000]; 
int LCONF[4000][12];

