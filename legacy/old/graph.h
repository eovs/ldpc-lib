#define NMAX 150
#define KMAX NMAX     
#define RMAX NMAX
#define CMAX 200
#define DMAX 1000
#define NDMAX 64
#define GMAX 10

//#define LARGE

#define MAX_EQ 600000  // maximum number of equations

#define MAX_N 200   // maximum length of base code

//Globals for code analysis
extern int CNodes[MAX_EQ];     // current nodes
extern int f[MAX_EQ];          // auxilary storage used in find routines
extern int ACTIVE[MAX_EQ];     // auxilary 
extern int DEPTH[MAX_EQ];      // 
extern int PBranch[MAX_EQ];    // previous branch
extern int WT[MAX_EQ];         //,WC[MAX_EQ];

// For fast equation test

extern int DTABLE[DMAX][NDMAX];
extern int FLAGS[MAX_EQ];
//int FLAGS_ARRAY[MAX_EQ][GMAX];

//int T[MAX_EQ][CMAX];
extern int T[CMAX][CMAX];

extern int H[CMAX][CMAX];
extern int	HH[CMAX][CMAX];
extern int	HTB[RMAX][NMAX];
extern int	HM[RMAX][NMAX]; //Mask 
extern int HD[CMAX][CMAX];
// int    HQ[CMAX][CMAX];
// Variables for equations
extern int	EQ[MAX_EQ][CMAX],
//	Length[MAX_EQ],
//Depth[MAX_EQ],
	CYCLES[MAX_EQ][2],
//  Index[MAX_EQ],
	TREE[MAX_EQ][2];
	//tree[MAX_EQ][2];
extern int   RESULTSG[CMAX];
//    RESULTSA[NMAX][CMAX];
extern int Graph[NMAX][CMAX];
extern int Branches[NMAX][CMAX];
extern int Numbers[NMAX]; 
extern int Path[MAX_EQ][CMAX]; // Path to node (part of equation)
extern int tree[CMAX*MAX_EQ][2];    // temporal tree 

extern int d_prev;
extern int path[NMAX];

void  bin_matr_print(int n, int k, int f);
int tanner(int n,int r);
//int   h2incidence(int n, int *r);
//void  inc2graph (int n, int r);
void  H2graph (int n, int r);
//int   dijkstra(int n);
int   graph_girth(int n, int r, int b);
void  hpoly2tb(int c, int b, int L);
//int   g2h(int n, int k, int inform_set[], int check_set[]);
//int   h2g(int n, int r, int inform_set[], int check_set[]);
//void  h2multigraph(int n,int r);
//void  hamming_table(int m);
//int   hamming_weight(unsigned _int64 x);
//int   dmin_G(int n, int k, int k0);

//void search_graph(void);
int sign(int x);
void ssort(int n, int A[], int V[]);
void sssort(int n, int A[]);
int nextsorted(int b[], int n, int M );
int nextnonsorted(int b[], int n, int M );

int eq_generation(int n,int r, int g, int *Eqs, int *Nodes);
int min_module(int a[],int eqs, int nodes, int *Mmax);
int check_eq(int a[],int eqs, int nodes);

int divisors(int x, int D[]);
void gen_divtable();
int mini(int a, int b);
int maxi(int a, int b);
