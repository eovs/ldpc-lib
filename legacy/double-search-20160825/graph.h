#define MMAX 400
#define CMAX 80
#define NMAX CMAX*MMAX
#define KMAX NMAX/2     
#define RMAX NMAX*3/4

#define WMAX CMAX/2
#define DMAX 1000
#define NDMAX 64

#define MAX_EQ 100000  // maximum number of equations

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
int HD[CMAX][CMAX]; // monomial matrix
int HD0[CMAX][CMAX]; //  copy



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
int msg[NMAX];

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
int LCONF[4000][12], LCONF0[4000][12];;

// ROUTINES

void  bin_matr_print(int n, int k, int f);
int tanner(int n,int r);
void  H2graph (int n, int r);
//int   graph_girth(int n, int r);
//void  hpoly2tb(int c, int b, int L);

//void search_graph(void);
int sign(int x);
void ssort(int n, int A[], int V[]);
void sssort(int n, int A[]);

int eq_generation(int n,int r, int g, int *Eqs, int *Nodes);
int min_module(int a[],int eqs, int nodes, int *Mmax);
int check_eq(int a[],int eqs, int nodes);
//the same for subcode

int divisors(int x, int D[]);
void gen_divtable();
int mini(int a, int b);
int maxi(int a, int b);

void hd2cv2(int b, int c, int M, int row_weights[], int col_weights[]);
void randnn(float x[], int n);
int bp_decod_lm(double soft[], int n, int r, int cweight[], int rweight[], int maxiter);

int rand_base_matrix(int r, int r0, int n, int msnd, int msnd0, int DD0[], int DD[], int trials[], int score[], int init);
int list_conf(int n, int s);
int list_conf_n(int n, int s);
int generate_code(int r, int n, int g, int M, int N);
int bp_simulation(int b, int b0, int c, int M, int maxiter, int Nerr, int Nexp, float SNR, float PR_ERR[]);
double absd (double x);
float absf(float x);


int count_loops(int n,int r, int score[], int dd[]);
 
int rand_base_matrix_g(int r0, int r, int n, int msnd0, int msnd, int DD0[], int D[], int trials[], int score[], int init);
int random_codeword(int b, int c, int M);
int QCencoder(int b, int c, int M);
int syndrome(int row_weights[], int b, int c, int M);