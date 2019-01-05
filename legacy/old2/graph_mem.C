#include "graph.h"

//Globals for code analysis
int CNodes[MAX_EQ];     // current nodes
int f[MAX_EQ];          // auxilary storage used in find routines
int ACTIVE[MAX_EQ];     // auxilary 
int DEPTH[MAX_EQ];      // 
int PBranch[MAX_EQ];    // previous branch
int WT[MAX_EQ];         //,WC[MAX_EQ];

// For fast equation test

int DTABLE[DMAX][NDMAX];
int FLAGS[MAX_EQ];
//int FLAGS_ARRAY[MAX_EQ][GMAX];

//int T[MAX_EQ][CMAX];
int T[CMAX][CMAX];

int H[CMAX][CMAX];
int	HH[CMAX][CMAX];
int	HTB[RMAX][NMAX];
int	HM[RMAX][NMAX]; //Mask 
int HD[CMAX][CMAX];
// int    HQ[CMAX][CMAX];
// Variables for equations
int	EQ[MAX_EQ][CMAX],
//	Length[MAX_EQ],
//Depth[MAX_EQ],
	CYCLES[MAX_EQ][2],
//  Index[MAX_EQ],
	TREE[MAX_EQ][2];
	//tree[MAX_EQ][2];
int   RESULTSG[CMAX];
//    RESULTSA[NMAX][CMAX];
int Graph[NMAX][CMAX];
int Branches[NMAX][CMAX];
int Numbers[NMAX]; 
int Path[MAX_EQ][CMAX]; // Path to node (part of equation)
int tree[CMAX*MAX_EQ][2];    // temporal tree 

int d_prev;
int path[NMAX];

