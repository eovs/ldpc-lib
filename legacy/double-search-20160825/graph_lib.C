#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "graph.h"
#include <conio.h>

/*************************************************/
int    mini(int    a, int    b) {if (a<b) return a; else return b;}
int    maxi(int    a, int    b) {if (a<b) return b; else return a;}
double mind(double a, double b) {if (a<b) return a; else return b;}
double maxd(double a, double b) {if (a<b) return b; else return a;}
float  minf(float  a, float  b) {if (a<b) return a; else return b;}
float  maxf(float  a, float  b) {if (a<b) return b; else return a;}
/*************************************************/
void ssort(int n, int A[], int V[])
{
//function [L,V]=ssort(A)
// very simple sorting
int i,x,flag;
flag=1;
for (i=0; i<n; i++) V[i]=i;
while (flag)
{
    flag=0;
    for (i=0; i<n-1; i++)
        if (A[i+1]<A[i])
		{
            flag=1;
            x=A[i]; A[i]=A[i+1]; A[i+1]=x;
            x=V[i]; V[i]=V[i+1]; V[i+1]=x;
		}
}
} // ssort
/*************************************************/
void sssort(int n, int A[])
{
//function L=ssort(A)
// very simple sorting  INDECES FROM 1!!!!!
int i,x,flag;
flag=1;
while (flag)
{
    flag=0;
//    for (i=0; i<n-1; i++)
	    for (i=1; i<n; i++)
        if (A[i+1]<A[i])
		{
            flag=1;
            x=A[i]; A[i]=A[i+1]; A[i+1]=x;
		}
}
} // ssort
/*************************************************/
int sign(int x)
{
	if (x<0) return -1;
	else if (x==0) return 0;
	return 1;
}

/*************************************************/
float absf(float x)
{
	if (x<0) return -x; else return x;
}
/*************************************************/
double absd (double x)
{
	if (x<0) return -x; else return x;
}


/*
void hamming_table(int m)
{
	int I,J,M;
	WH[0]=0; M=1; 
	for (I=1; I<=m; I++)
	{
		for (J=M; J<2*M; J++) WH[J]=WH[J-M]+1;
		M=M+M;
	}
}
*/
/*************************************************/
/*
int hamming_weight(unsigned _int64 x)
{
	int wt=0;
	while (x!=0)
	{
		wt+=WH[x&MASK]; x>>=16;
	}
	return wt;
}
*/
/*************************************************/
void  bin_matr_print(int n, int k, int f)
{
FILE *f_matr;
int i,j;
errno_t err;

if (f)	err=fopen_s(&f_matr,"matr.dat", "at");
for (i=0; i<k; i++)
{ 
	//word2bin(g[i],n,row);
	for (j=0; j<n; j++) 
	{
		printf("%2d",H[i][j]); 
		if (f) fprintf(f_matr,"%2d",H[i][j]); 
	}
    printf("\n");
	if (f) fprintf(f_matr,"\n"); 
}
printf("\n");
if (f)
{
	fprintf(f_matr,"\n"); 
	fclose(f_matr);
}
}

/*************************************************/

void H2graph (int n, int r)
{
	int i,j,h;
// For bipartite graph 
// transforms incidence marix A to 
// graph description [Graph, Numbers]
// For each node array desription contains:
//   number of achievable nodes, Numbers
//  - list of achievable nodes,   Graph 
// By default(not necessarily) A is parity check matrix H
// i.e. size(A)=r,n, max row weight w
for (i=0; i<r; i++) // over nodes
{
	Numbers[i]=0;
	for (j=0; j<n; j++)
	{
		if (H[i][j]>=1)
		// check the column for more ones
		for (h=0; h<r; h++)
		if ((h!=i) && (H[h][j]>=1))
		{
			Graph[i][Numbers[i]]=h;
			Branches[i][Numbers[i]]=j+1;
			Numbers[i]++;
		}
	}
}
} //inc2graph

/*************************************************/
/*
void hpoly2tb(int c, int b, int L)
{
	// input HD
	// intermediate result is HTL = 2^HD
	// Therefore dgrees above 63 are not allowed
	// output HTB
    int i,j,h,l,wr,wc,wcc, 
		HT[CMAX][CMAX], 
	//_int64	
    	HTL[CMAX][CMAX];

	// in binary form
	for (i=0; i<b; i++)
	for (j=0; j<c; j++) HTL[i][j]=HD[i][j]%L; // MOD L 
	// first row
	wc=0; 
	for (l=0; l<L; l++)
	{
		for (i=0; i<b; i++)
		for (j=0; j<c; j++) 
		{
			if (HTL[i][j]==l) HTB[i][wc+j]=1;
			else        	  HTB[i][wc+j]=0;
		}
		wc=wc+c;
	}
	// other rows
    wr=0; 
	for (l=1; l<L; l++)
	{   // loop over other rows
		wc=0;
		for (h=0; h<L; h++)
		{
			// read block
			for (i=0; i<b; i++)
			for (j=0; j<c; j++) HT[i][j]=HTB[wr+i][wc+j];
			// write it
			if (wc+c>=L*c) wcc=0; 
			else wcc=wc+c; 
			for (i=0; i<b; i++)
				for (j=0; j<c; j++) HTB[wr+b+i][wcc+j]=HT[i][j]; 
			wc=wcc;
		}
		wr=wr+b;
	}
} //hpoly2Tb
*/
/*************************************************/
int tanner(int n,int r)
// constructs Tanner matrix
{
int IC, i, j,s;
//HH=zeros(n+r, sum(sum(H)));
s=0;
for (i=0; i<r; i++)
	for (j=0; j<n; j++) s+=H[i][j];
for (i=0; i<(r+n); i++)
	for (j=0; j<s; j++) HH[i][j]=0;

IC=0; // current column of TH 
for (i=0; i<n; i++)
{
    // along column
	for (j=0; j<r; j++)
	{
		if (H[j][i]==1)
		{
			HH[i][IC]=1;
			HH[n+j][IC]=1;
			IC++;
		}
	}
}
return s;
}// void
   
/*************************************************/
int eq_generation(int n,int r, int g, int *Eqs, int *Nodes)
{
// H is parity-check for base graph
// OUTPUT:
// returns 1 if everything works fine
// otherwise retutns 0;
// EQ[eqs][n]   equations in normal form
// TREE[nodes][2] tree for fast computing
// CYCLES[eqs][2] equaltions in compact form
// Length[eqs]    cycle length
// param=[eqs,nodes] 

int i, j, h, I, J, L; 
//int L,m, Lmin, Lmax; 

// Variables for EQ constructing
int N[CMAX]; // # nodes at tree levels
int eqsp, eqs, nodesp, nodes, source, parent, child,branch; 
int s1,s2,w, flag;
int new_eq[NMAX];       // 

for (i=0; i<r; i++)
for (j=0; j<n; j++) T[i][j]=H[i][j];
H2graph(n,r);

eqsp=0; eqs=0; nodesp=0;

for (source=0; source<r; source++)
{
   // Collection of cycles from source to source
   //Path=zeros(N0,n);    //  Path to node (part of equation)
   CNodes[0]=source;      // PNodes(1)=0;
   N[0]=0;  N[1]=1;
   PBranch[0]=0;
   nodes=0; 
   DEPTH[0]=0;
   for (L=0; L<g/2; L++)
   {
	   for (j=N[L]; j<N[L+1];j++) // over nodes of previous level
	   {
            parent=CNodes[j];               // parent node       
			if (parent>=source)    // TO AVOID IDENTICAL CYCLES
            for (h=0; h<Numbers[parent]; h++)
			{
                child=Graph[parent][h];      // child node
				//if (child<source) break;    // TO AVOID IDENTICAL CYCLES
                branch=Branches[parent][h];  // branch parent->child
                // check if not back
                if  (branch!=PBranch[j])
				{   //  Good. Add to the list 
                    nodes++;
                    CNodes[nodes]=child;
     				for (J=0; J<n; J++) Path[nodes][J]=Path[j][J];
                    Path[nodes][branch-1]=Path[nodes][branch-1]+sign(child-parent);
                    PBranch[nodes]=branch;
					TREE[nodes][0]=j;
					TREE[nodes][1]=sign(child-parent)*branch;
                    DEPTH[nodes]=L+1;
                    // Check for loops
					for(I=0; I<nodes; I++)
					{
                        if (CNodes[I]==CNodes[nodes])
						{
                            // Create new equation
							w=0;
							for (J=0; J<n; J++) 
							{
								new_eq[J]=Path[nodes][J]-Path[I][J]; 
								w+=abs(new_eq[J]);
							}
							if ((w==0)&& (DEPTH[I]+L+1<g)) 
							{
								printf("BALANCED CYCLE"); 
								//for (J=0; J<n; J++) printf(" %d", new_eq[J]); printf("\n");
								return 0;
							}
                            //is it really new equation?
							if ((w>0) && (DEPTH[I]+L+1<g))
							{
								flag=1;
//								Length_new=DEPTH[nodes]+DEPTH[I];
								for (i=0; i<eqs; i++)
								{
									s1=0; 
									while ((EQ [i][s1]==new_eq[s1]) &&(s1<n)) s1++;
									if (s1<n)
									{
										s2=0;
										while ((EQ [i][s2]==-new_eq[s2])&&(s2<n)) s2++;
									}
							
									
									if ((s1==n)||(s2==n))
									{
										flag=0; 
//										if (Length_new<Length[i]) Length[i]=Length_new;
                                        break;
									}
								}
                                if (flag)
								{
                                  for (J=0; J<n; J++)EQ[eqs][J]=new_eq[J];
  //                                Length[eqs]=Length_new;
                                  CYCLES[eqs][0]=nodes;
								  CYCLES[eqs][1]=I;
								  eqs++; //printf("%d,",eqs);
								}
							}
						} 
					}  // for I 
				}//
			}        // for h
	   } // for j
	   N[L+2]=nodes+1;   //     printf("\n Level:%d, eqs:%d",L,eqs);
   }
   nodes++;
   if (eqs>eqsp)
   {
   // new equations found
        for (J=0; J<nodes; J++)
		{
			TREE[J][0]+=nodesp;        
            tree[nodesp+J][0]=TREE[J][0];        
			tree[nodesp+J][1]=TREE[J][1];        
		}
		for (J=0; J<eqs; J++)
		{
			CYCLES[eqsp+J][0]+=nodesp;
			CYCLES[eqsp+J][1]+=nodesp;
		}
		eqsp=eqs;
        nodesp+=nodes;
   }
//   printf(" source, eqs, nodes:");
//    printf(" %d %d %d\n",source, eqs, nodes);
}

for (J=0; J<nodesp; J++)
{
	TREE[J][0]=tree[J][0];
	TREE[J][1]=tree[J][1];
}

// refine tree
// STEP1; Find active nodes of the tree
//ACTIVE=zeros(1, size(TREE,1));
nodes=nodesp;
eqs=eqsp;
for (J=0; J<nodes; J++) ACTIVE[J]=0;
for (J=0; J<eqs; J++)
{
	ACTIVE[CYCLES[J][0]]=1;
	ACTIVE[CYCLES[J][1]]=1;
}


// find
s1=0; 
for (J=0; J<nodes; J++) if (ACTIVE[J]) {f[s1]=J; s1++;}
flag=1;
while (flag)
{
    // mark nodes from tree
	w=0; for (J=0; J<s1; J++) {ACTIVE[TREE[f[J]][0]]=1; 
	//printf("%d : %d ",J, TREE[f[J]][0]); 
	w++;}
    // find again
	s2=0; for (J=0; J<nodes; J++) if (ACTIVE[J]) {f[s2]=J; s2++;}
    if (s1==s2) flag=0;
//	printf(" %d %d\n",s1,s2);
	s1=s2;
	//printf("\n");
}    

// Renumeration
// find
w=0; for (J=0; J<nodes; J++) if (ACTIVE[J]) {f[w]=J; w++;}
for (J=0; J<nodes; J++) ACTIVE[J]=0;
for (i=0; i<w; i++)  ACTIVE[f[i]]=i;
for (J=0; J<eqs; J++)
{
	CYCLES[J][0]=ACTIVE[CYCLES[J][0]];     
	CYCLES[J][1]=ACTIVE[CYCLES[J][1]];     
}
for (J=0; J<nodes; J++)
	TREE[J][0]=ACTIVE[TREE[J][0]];

for (J=0; J<w; J++)
{
	TREE[J][0]=TREE[f[J]][0];
	TREE[J][1]=TREE[f[J]][1];
}

nodes=w;
*Eqs=eqs; 
*Nodes=nodes;

return 1;

} //eq_generation

/*************************************************/

int min_module(int a[],int eqs, int nodes, int *Mmax)
{
// computes girth-length profile for a given input 
// voltage vector a[1...n]
int i, j, h, new_code=0, w, nd, D[DMAX], wmax=0;

for (i=0; i<MAX_EQ; i++) 
	FLAGS[i]=0;

// Computing weights of tree nodes and cycles
j=0; // # already computed tree nodes
if (eqs==0) {*Mmax=0; return 0;}
for (i=0; i<eqs; i++)
{
	w=abs(WT[CYCLES[i][0]]-WT[CYCLES[i][1]]);
	if (w>wmax) wmax=w;
	{
	//	w=abs(w);
		if (w<NDMAX)
		{
			nd=DTABLE[w][0];
			for (h=1; h<=nd; h++) FLAGS[DTABLE[w][h]]++;
		}
		else
		{
			nd=divisors(w,D); 
			for (h=1; h<=nd; h++) FLAGS[D[h]]++;
		}
	}
      
}
h=2;
while (FLAGS[h]!=0) h++; 
i=wmax+1;
while (FLAGS[i]==0) i--; 
*Mmax=i+1;
return h;


return new_code;
}// min_mod

/*************************************************/

int check_eq(int a[],int eqs, int nodes)
{
// voltage vector a[1...n]
int i, j, J, check=1;

// Computing weights of tree nodes and cycles
j=0; // # already computed tree nodes
for (i=0; i<eqs; i++)
{
    for (J=j; J<=CYCLES[i][0]; J++)
	{	
		WT[J]=0;
        if (TREE[J][1]!=0)
           WT[J]=WT[TREE[J][0]]+a[abs(TREE[J][1])]*sign(TREE[J][1]);
	}
    j=CYCLES[i][0]+1;
//    WC[i]=WT[CYCLES[i][0]]-WT[CYCLES[i][1]]; 
    if (WT[CYCLES[i][0]]==WT[CYCLES[i][1]]) return -(i+1);  //gfree<g
}

return check;
}// min_mod
/*************************************************/
int divisors(int x, int D[])
{
 // generate3s list of divisors D and number of divisors nd

	int i, y, nd, maxd;
  
	//for (i=0; i<=NDMAX; i++) D[i]=0;
	nd=0;
	maxd =x;
	i=2;
	while (i<maxd)
	{
		if ((x%i)==0)
		{
			nd++;
			D[nd]=i;
			y=x/i;
			if (y>i)
			{
				nd++;
				D[nd]=y;
				maxd=y;
			}
		}
		i++;
	}
	nd++;
	D[nd]=x;
	D[0]=nd;
	return nd;
}

/*************************************************/
void gen_divtable()
{
	int i, j, D[NDMAX], x;
	
	for (i=0; i<DMAX; i++)
	{
		x=divisors(i,D); 
		for (j=0; j<=D[0]; j++) DTABLE[i][j]=D[j]; 
	} 
}
/*************************************************/
void randnn(float x[], int n) 
{// Box-Muller transform
float u, v, w; 
int i, j=0; 
for (i=0; i<(n+1)/2; i++)
{
	do {
          u = (float) 2.0 * ((float) rand()/RAND_MAX) - (float) 1.0;
          v = (float) 2.0 * ((float) rand()/RAND_MAX) - (float) 1.0;
          w = u * u + v * v;
       } 
	while ( w >= 1.0 );

    w = sqrtf( ((float)-2.0 * logf( w) ) / w );
    x[j] = u * w; j++; 
    x[j] = v * w; j++; 
}
//printf(" %f",s/n);
}
/*************************************************/
int binom(int n, int m)
{
	// calculates number of combinations of length n weight m
	int i, c=1;
	for (i=0; i<m; i++) c=c*(n-i)/(i+1);
	return c;
}
/*************************************************/
int randi(int imax)
  {
     // Generate random number in the half-closed interval [0, imax] 
   return (int) ((double)rand() / (RAND_MAX + 1) * imax); 
  }

/*************************************************/
int randcomb(int n, int w, int pos[])
{
	// generates nonzero positions of random combination of length n and weight w
    int i, j=0,t;
	if (w>=n) {for (i=0; i<n; i++)	pos[i]=i; return w;} 
	for (i=0; i<n-1; i++)
	{   
		t=((n-i)*rand()/RAND_MAX);
		if (t<w)
		{
			pos[j]=i;	j++;  w--; if (w==0) return j;
		}
	}
	if (w==1) pos[j]=n-1; 	
	return j;
}
/*************************************************/
int next_permut(int n, int m, int P[])
{
	// Generates next combination of m elements from n
	// returns 1 if last permutation
	int i, j;
	j=m; 
	while (P[j-1]==n+j-m) 
	{
		j--; 
		if (j==0)  return 1; 
	}
	P[j-1]++;    
    for (i=j; i<=m-1; i++) P[i]=P[i-1]+1;
	return 0;
}
/*************************************************/
int next_permut0(int n, int m, int P[])
{
	// Generates next combination of m elements from n
	// returns 1 if last permutation
	int i, j;
	j=m-1; 
	while (P[j]==n+j-m) 
	{
		j--; 
		if (j<0)  return 1; 
	}
	P[j]++;    
    for (i=j; i<m-1; i++) P[i+1]=P[i]+1;
	return 0;
}
/*************************************************
int rand_base_matrix(int r, int n, int msnd, int DD[], int trials[], int score[], int init)
// generates random base LDPC matrix HB with degree distribution DD
// upper triangular (dual diagonal)  parity-check part
// DD(i) =  number of columns of weight i
// msnd   = max symbol node degree
{
	int i, j, h, I, ii, ij, sum, t, g, nt, rt, 
		N1, N2, N3,      // overall best
		T1, T2, T3,      //best for  for code
		t1, t2, t3, w,   //increments
		best1, best2, best3, //best for column 
		a[CMAX], besta[CMAX], dd[CMAX], flag,
		girth[CMAX], ace[CMAX], num_cycles[CMAX];
	int corcoef0[CMAX][CMAX],corcoef[CMAX][CMAX]; 
	// Check input
	sum=0;
	for (i=1; i<=r; i++) sum+=DD[i];
	if ((sum != n) || (DD[1] > 0) )
	{
		printf("sum(DD) ~= n"); return -1;
	}

//----------------------------------
// Prepare arrays
for (i=0; i<n; i++) girth[i]=ace[i]=num_cycles[i]=100000;
if (init)
{
	for (i=0; i<=r; i++) ncombs[i]=binom(r,i); 
	combs3[0][0]=0; combs3[0][1]=1; combs3[0][2]=2;
	for (i=1; i<ncombs[3]; i++)
	{
		a[0]=combs3[i-1][0]; a[1]=combs3[i-1][1];a[2]=combs3[i-1][2];
		next_permut(r,3,a);
		combs3[i][0]=a[0]; combs3[i][1]=a[1]; combs3[i][2]=a[2];
	}
	// w = 4
	combs4[0][0]=0; combs4[0][1]=1; combs4[0][2]=2; combs4[0][3]=3;	
	for (i=1; i<ncombs[4]; i++)
	{
		a[0]=combs4[i-1][0]; a[1]=combs4[i-1][1];a[2]=combs4[i-1][2];a[3]=combs4[i-1][3];
		next_permut(r,4,a);
		combs4[i][0]=a[0]; combs4[i][1]=a[1]; combs4[i][2]=a[2];  combs4[i][3]=a[3];
	}
	// w = 5
	for (i=0; i<5; i++) combs5[0][i]=i; 	
	for (i=1; i<ncombs[5]; i++)
	{
		for (j=0; j<5; j++) a[j]=combs5[i-1][j];
		next_permut(r,5,a);
		for (j=0; j<5; j++) combs5[i][j]=a[j];
	}
	// w = 6;
	for (i=0; i<6; i++) combs6[0][i]=i; 	
	for (i=1; i<ncombs[6]; i++)
	{
		for (j=0; j<6; j++) a[j]=combs6[i-1][j];
		next_permut(r,6,a);
		for (j=0; j<6; j++) combs6[i][j]=a[j];
	}
	// w = 7;
	for (i=0; i<7; i++) combs7[0][i]=i; 	
	for (i=1; i<ncombs[7]; i++)
	{
		for (j=0; j<7; j++) a[j]=combs7[i-1][j];
		next_permut(r,7,a);
		for (j=0; j<7; j++) combs7[i][j]=a[j];
	}
	// w = 8;
	for (i=0; i<8; i++) combs8[0][i]=i; 	
	for (i=1; i<ncombs[8]; i++)
	{
		for (j=0; j<8; j++) a[j]=combs8[i-1][j];
		next_permut(r,8,a);
		for (j=0; j<8; j++) combs8[i][j]=a[j];
	}
	for (i=0; i<9; i++) combs9[0][i]=i; 	
	for (i=1; i<ncombs[9]; i++)
	{
		for (j=0; j<9; j++) a[j]=combs9[i-1][j];
		next_permut(r,9,a);
		for (j=0; j<9; j++) combs9[i][j]=a[j];
	}
	for (i=0; i<10; i++) combs10[0][i]=i; 	
	for (i=1; i<ncombs[10]; i++)
	{
		for (j=0; j<10; j++) a[j]=combs10[i-1][j];
		next_permut(r,10,a);
		for (j=0; j<10; j++) combs10[i][j]=a[j];
	}
	for (i=0; i<11; i++) combs11[0][i]=i; 	
	for (i=1; i<ncombs[11]; i++)
	{
		for (j=0; j<11; j++) a[j]=combs11[i-1][j];
		next_permut(r,11,a);
		for (j=0; j<11; j++) combs11[i][j]=a[j];
	}
	for (i=0; i<12; i++) combs12[0][i]=i; 	
	for (i=1; i<ncombs[12]; i++)
	{
		for (j=0; j<12; j++) a[j]=combs12[i-1][j];
		next_permut(r,12,a);
		for (j=0; j<12; j++) combs12[i][j]=a[j];
	}
	for (i=0; i<13; i++) combs13[0][i]=i; 	
	for (i=1; i<ncombs[13]; i++)
	{
		for (j=0; j<13; j++) a[j]=combs13[i-1][j];
		next_permut(r,13,a);
		for (j=0; j<13; j++) combs13[i][j]=a[j];
	}
	for (i=0; i<14; i++) combs14[0][i]=i; 	
	for (i=1; i<ncombs[14]; i++)
	{
		for (j=0; j<14; j++) a[j]=combs14[i-1][j];
		next_permut(r,14,a);
		for (j=0; j<14; j++) combs14[i][j]=a[j];
	}
	for (i=0; i<15; i++) combs15[0][i]=i; 	
	for (i=1; i<ncombs[15]; i++)
	{
		for (j=0; j<15; j++) a[j]=combs15[i-1][j];
		next_permut(r,15,a);
		for (j=0; j<15; j++) combs15[i][j]=a[j];
	}
	for (i=0; i<16; i++) combs16[0][i]=i; 	
	for (i=1; i<ncombs[16]; i++)
	{
		for (j=0; j<16; j++) a[j]=combs16[i-1][j];
		next_permut(r,16,a);
		for (j=0; j<16; j++) combs16[i][j]=a[j];
	}

}	// init
//---------------------------------------
for (i=0; i<r; i++) for (j=i; j<r; j++) corcoef0[i][j]=0;
// vector of node symbol weights
sum=0;
for (i=1; i<=msnd; i++) 
{
	for (h=0; h<DD[i]; h++) dd[sum+h]=i;
	sum+=DD[i];
}
// non-random columns of weight 2
for (i=0; i<r; i++) for (j=0; j<n; j++) T[i][j]=0;
for (i=0; i<DD[2]; i++)
{
	j=i%r;
	T[j][j]=1;
	if (i<=r) h=(i+1)%r; 
	else h=(i+r/2)%r; T[h][j]=1;
	T[h][j]=1;
    corcoef0[j][h]=1;
}
// random columns of higher weight
N1=N2=N3=100000;  // Best code performance    
for (t=0; t<trials[0]; t++)
{
	for (i=0; i<r; i++) for (j=i; j<r; j++) corcoef[i][j]=corcoef0[i][j];
	for (i=0; i<r; i++) for (j=DD[2]; j<n; j++) T[i][j]=0;
    T1=T3=T2=0; // current code performance performance 
    for (h=DD[2]; h<n; h++)
	{
        w=dd[h];
		// generate set of random columns
		randcomb(ncombs[w], trials[1], col_set);
        best1=best2=best3=10000;  // best column performance]
	    for (j=0; j<mini(ncombs[w],trials[1]); j++)
		{
           for (i=0; i<r; i++) a[i]=0;
           switch (w)
		   {
               case 3:
				   for (i=0; i<w; i++) a[i]=combs3[col_set[j]][i]; 
               break;
			   case 4:
				   for (i=0; i<w; i++) a[i]=combs4[col_set[j]][i];
               break;
			   case 5:
			   	   for (i=0; i<w; i++) a[i]=combs5[col_set[j]][i];
			   break;
               case 6:
				   for (i=0; i<w; i++) a[i]=combs6[col_set[j]][i];
			   break;
			   case 7:
				   for (i=0; i<w; i++) a[i]=combs7[col_set[j]][i];			   
			   break;
			   case 8:
				   for (i=0; i<w; i++) a[i]=combs8[col_set[j]][i];
			   break;
	 	       case 9:
				   for (i=0; i<w; i++) a[i]=combs9[col_set[j]][i];
			   break;
	           case 10:
				   for (i=0; i<w; i++) a[i]=combs10[col_set[j]][i];
			   break;
			   case 11:
				   for (i=0; i<w; i++) a[i]=combs11[col_set[j]][i];
			   break;
			   case 12:
				   for (i=0; i<w; i++) a[i]=combs12[col_set[j]][i];
			   break;
	 	       case 13:
				   for (i=0; i<w; i++) a[i]=combs13[col_set[j]][i];
			   break;
	           case 14:
				   for (i=0; i<w; i++) a[i]=combs14[col_set[j]][i];
			   break;
			   case 15:
				   for (i=0; i<w; i++) a[i]=combs15[col_set[j]][i];
			   break;
			   case 16:
				   for (i=0; i<w; i++) a[i]=combs16[col_set[j]][i];
			   break;
     	   }
           // column-wise correlation contribution
//		   a[0]=0; a[1]=1; a[2]=7;
           t1=t2=t3=0;
           for (I=0; I<h; I++)
		   {
			   sum=0; for (i=0; i<w; i++) sum+=T[a[i]][I]; 
			   if (I>=h-w+1) t1+=sum; // instead of I>=h-2
			   t2+=binom(sum,2); t3+=binom(sum,3);
			}
           // row-wise correlation contribution
           // only rows with new 1's contribute anything
           for (ii=1; ii<w; ii++)
			   for (ij=0; ij<ii; ij++)
			   {
                  sum=corcoef[a[ij]][a[ii]]+1; // new row-wise correlation
				  if (sum>=3) {  t2+=sum-1; t3+=binom(sum-1,2);}
				
			   }
           // choose the best column
		   if ((t3<best3) || ((t3==best3)&&(t2<best2)) || ((t3==best3)&&(t2==best2)&&(t1<best1)))
           { 
			    for (i=0; i<w; i++) besta[i]=a[i];
				best1=t1;
                best3=t3;
                best2=t2; 
		   }
		}  //for j
		// Check for # cycles:
		// Copy:
		for (ii=0; ii<r; ii++) for (ij=0; ij<h; ij++) H[ii][ij]=T[ii][ij];
		// Add column
  		for (ii=0; ii<r; ii++) H[ii][h]=0;
  		for (ii=0; ii<w; ii++) H[besta[ii]][h]=1;
		rt=h+r+1; 	nt=tanner(h+1,r);
		for (ii=0; ii<rt; ii++) for (ij=0; ij<nt; ij++) H[ii][ij]=HH[ii][ij];
		g=count_loops(nt, rt, score,dd);
		flag=1;
		if (g>girth[h]) flag=0;
		else if (g<girth[h])
		{
			girth[h]=g;
			ace[h]=num_cycles[h]=1000000;
		}
		else //		if (g==girth[h])
			if (score[1]>ace[h]) flag=0;
			else if (score[1]<ace[h])
			{
				ace[h]=score[1];
				num_cycles[h]=1000000;
			}
			else //if (score[1]==ace[h])
				if (score[0]>num_cycles[h]) flag=0;
				else if (score[0]<num_cycles[h]) num_cycles[h]=score[0];
		// add new column
		if (flag) {
        T3+=best3;         T2+=best2; 		T1+=best1;
		for (ii=1; ii<w; ii++)
            for (ij=0; ij<ii; ij++) 
                 corcoef[besta[ij]][besta[ii]]+=1 ;
  		for (i=0; i<w; i++) T[besta[i]][h]=1;
		}
	}//  for cn
 	//printf("T1=%d, T2=%d, T3=%d, t=%d \n", T1, T2, T3, t);
	if ((T3<N3) || ((T3==N3)&&(T2<N2)) || ((T3==N3)&&(T2==N2)&&(T1<N1)))
	{
		
        for (i=0; i<r; i++) for (j=0; j<n; j++) H[i][j]=HM[i][j]=T[i][j];
        
	//	rt=n+r; 	nt=tanner(n,r);
	//	for (i=0; i<rt; i++) for (j=0; j<nt; j++) H[i][j]=HH[i][j];
	//	g=count_loops(nt, rt, score);
		
		N3=T3;  N2=T2; N1=T1;
		printf("g=%d, cycles=%d, ace=%d \n", girth[n-1], num_cycles[n-1],ace[n-1]);
		printf("N1=%d, N2=%d, N3=%d, t=%d \n", N1, N2, N3, t);
    //    disp(N3)
	}
} // % for t
score[0]=N2; score[1]=N3;

//---------------------------------------------
// For low-complexity encoding:
h=DD[2];
while (HM[0][h]=0) h++;
//permute columns h and r-1
for (i=0; i<r; i++)
{
	j=HM[i][h];
	HM[i][h]=HM[i][r-1];
	HM[i][r-1]=j;
}
// Transform HM to mask-matrix: first components of tows/columns=2
for (i=0; i<r; i++)
{
	j=0; while (HM[i][j]==0) j++; HM[i][j]=2;
}
for (i=0; i<n; i++)
{
	j=0; while (HM[j][i]==0) j++; HM[j][i]=2;
}
return 0;
}

*************************************************/
int list_conf(int n, int s)
{
	// generates list of possible configurations of 
	// weights enumerated in 1...s
	// or representing n as a sum of  1...s summands 
	int i, j, h, S, J, P[10];
	S=1;
	for (i=1; i<s; i++) S+=binom(n-1,i);
	for (i=0; i<S; i++) for (j=0; j<=s; j++) LCONF[i][j]=0;
	LCONF[0][0]=n;
	J=0;
	for (i=1; i<s; i++)
	{
		for (j=0; j<i; j++) P[j]=j;
		for (j=0; j<binom(n-1,i); j++)
		{
			J++;
			LCONF[J][0]=P[0]+1;
			for (h=1; h<i; h++) LCONF[J][h]=P[h]-P[h-1];
			LCONF[J][h]=n-P[h-1]-1;
			next_permut(n-1, i, P);
		}
	}
	return J+1;
}
/*************************************************/
int list_conf_n(int n, int s)
{
	// generates list of possible configurations of 
	// weights enumerated in 1...s
	// or representing n as a sum of  1...s summands 
	int  j, h, S, P[20];
	S=binom(n-1,s-1);
	for (j=0; j<s; j++) P[j]=j+1; P[s]=n+1;
	for (j=0; j<binom(n-1,s-1); j++)
	{
		for (h=0; h<s; h++) LCONF[j][h]=P[h+1]-P[h];
		next_permut(n, s, P);
	}
	return S;
}
/*************************************************/
int generate_code(int r, int n, int g, int M, int N)
{
	int eqs, nodes;
	int i, j, I, J, nt, rt,  
	    a[NMAX], am[NMAX], wcol[CMAX], cumwcol[CMAX];
	//int map[NMAX]; // mapping form initial code to base code
	int ok, new_code;//, ans, t; 
	int Mmax, M1;// , flag_continue;
	int max_codes=1; //N=1000000;
	char letter;
	//FILE *flog;

// goto finish;
// column weights (for all, not only used columns)  
cumwcol[0]=0;
for (i=0; i<n; i++)
{
	wcol[i]=0;	for (j=0; j<r; j++)	wcol[i]+=T[j][i];
	cumwcol[i+1]=cumwcol[i]+wcol[i];
}
// Generate mask
I=1;
for (i=0; i<n; i++) for (j=0; j<r; j++)
if (HM[j][i]>0) 
{	
	am[I]=HM[j][i];	I++;   // mask
}	
// First set of equations

// obtain base matrix of proper size for constructing equations
//fopen_s(&flog, "H.txt", "wt");
//I=0;
for (i=0; i<n; i++){ for (j=0; j<r; j++)	
{
		H[j][i]=T[j][i];
		//fprintf(flog,"%d ", H[j][i]); 
}
//fprintf(flog,"\n "); 
}
//fclose(flog);
//	I++;

	
printf("\n n=%d, r=%d, g=%d, M=%d \n",n, r, g, M); 
 
//   Tanner graph 
	rt=n+r; 	nt=tanner(n,r);  	                             //input H, output HH
	for (i=0; i<rt; i++) for (j=0; j<nt; j++) H[i][j]=HH[i][j];  // Tanner graph PCM 
	printf("\n Tanner graph size: nt=%d, rt=%d \n", nt, rt);
   // Generating first set of equations
	printf("Generating first equations  \n");
	ok=0;
	ok=eq_generation(nt,rt, g, &eqs, &nodes); 
	if (ok)  printf("ok: g=%d, eqs=%d, nodes=%d\n",g, eqs, nodes);
	else  return -g;

	for (i=1; i<=nt; i++) a[i]=0; 
	M1=M;

	// Random generating  n columns
//	max_codes=1; N=1000000;
//	err=fopen_s(&f_data, argv[2], "wt" ); fclose(f_data);
	if (1) 
	{
		printf("Creating random %d columns \n", n);
	    I=0; J=0;
		while ((I<N)&& (J<max_codes))
		{
			I++;  if(I%100000==0) printf("%d tested \n",I);
			for (i=1; i<=nt; i++) 
			{   
				switch(am[i])
				{
				case 1: a[i]=rand()%M; 	break;
				case 2: a[i]=0;         break;
				case 3: a[i]=1+rand()%(M-1); 	break;
				}
			}
			new_code=check_eq(a, eqs, nodes);
			
			if (new_code>0)	
	        {
		        Mmax=M;
		        new_code=min_module(a, eqs, nodes, &Mmax);
		        if (new_code<=M)        if (FLAGS[M]==0) J++;
			     
		        
	        }
		
			if (_kbhit())
			{
				letter=_getch();
				if (letter=='x') break;
				if (letter=='i') printf(" n= %d, %d tested %d found \n", n, I, J);
			}
     	}
	}
//}//2

// Write matrix to LOG


printf("\n  %d codes selected, %d tested\n", J, N);

I=0;
for (i=0; i<n; i++)          // column number
for (j=0; j<r; j++)        // index inside column  
{
	HD[j][i]=-1;
	if (HM[j][i]>0)
	{
		I++;
		HD[j][i]=a[I]%M;
	}
}

return J;
}

/*************************************************/

int bp_simulation(int b, int b0, int c, int M, int maxiter, int Nerr, int Nexp, float SNR, float PR_ERR[]) 
{
	int r, n, i,j,J, iter; 
	int cw_length;
	int row_weights[NMAX],    //  
		col_weights[NMAX], nse, 
		NSE, NUE, NDE;
	float bitrate, sigma;
	//double ch_out[NMAX];
	
r=b*M; 
n=c*M; cw_length=n; 
hd2cv2(b, c, M,  row_weights, col_weights);

	// Channel model
bitrate=(float) (c-b) / (float) c;
sigma=(float) sqrt(pow(10,(-SNR/10))/2/bitrate);  
			
// Message
for (i=0; i<n-r; i++) msg[i]=rand()&1;

for (i=0; i<b; i++) 
for (j=0; j<c; j++) HD0[i][j]=HD[i][j];   //keep original matrix
for (i=0; i<b-b0; i++) 
for (j=0; j<c-b0; j++) HD[i][j]=HD0[i+b0][j+b0];   


i=QCencoder(b-b0, c-b0, M);
if (i==1) {printf("Bad encoding 1"); return 0;}
if (i==-1) {printf("bad matrix 1"); return i;}


for (i=0; i<b; i++) 
for (j=0; j<c; j++) HD[i][j]=HD0[i][j];   

for (i=0; i<n-b0*M; i++) msg[i]=codeword[i];

i=QCencoder(b0, c, M);
if (i==1) {printf("Bad encoding 2"); return 0;}
if (i==-1) {printf("bad matrix 2"); return i;}



i=syndrome(row_weights, b, c, M);

//printf(" syndrome =%d codeword=",i);
//for (i=2400; i<2420; i++) printf("%d",codeword[i]);
//printf("\n");

NSE=NUE=NDE=0; //BSC=0;
J=0; //Nexp=2;
while ((NDE<Nerr) && (J<Nexp)) 
{
	if (_kbhit()) if (_getch()=='x') break;
	J++;
//	printf(" J=%d ", J);
	//Noising
	randnn(noise,n);
	//Transmitting
	for (i=0; i<n; i++) 
	{
		received[i]=-(float)2.0*(sigma*noise[i]+
		(float)2.0*(float)codeword[i]-(float)1.0)/(sigma*sigma);
		ch_out[i]=(double)received[i];
	}
	// Decoding
	//printf("ch_out=");
    //for (i=2400; i<2420; i++) printf("%d ",codeword[i]);
    //printf("\n");

	iter=bp_decod_lm(ch_out, hard, n,  r, col_weights, row_weights, maxiter);
	
	//printf("iter=%d \n", iter);
	// check result
	nse=0; 	for (i=0; i<n; i++) 
{		nse+=(hard[i]!=codeword[i]);
	//if (hard[i]!=codeword[i])
    // printf(" i= %d out=%d in=%d \n",i,ch_out[i]<0, codeword[i]);
	}
	//printf("after bp=");
    //for (i=2400; i<2420; i++) printf("%d ",codeword[i]);
    //printf("\n");
	if (nse>0) 
	{
		// counting errors
		NSE+=nse; NDE++; if (iter>=0) NUE++;
		printf("SNR=%5.3f,step=%d,s_ers=%d,f_ers=%d,u_ers=%d,BER=%5.3e,FER=%5.3e\n",SNR, J, NSE, NDE, NUE, (float)NSE/(J)/n,(float)NDE/(J));   
		/*
			fopen_s(&f, "error.txt", "wt");
for (i=0; i<b; i++) {for (j=0; j<c; j++) fprintf(f,"%d ", HD[i][j]); fprintf(f,"\n");}
		for (i=0; i<n; i++)
			fprintf(f,"%f \n", ch_out[i]);
		fclose(f);
		*/
		/*
		if (argc>3)
		{
			err=fopen_s(&fe,argv[3], "at");
			if (err!=0) printf("Problem with opening 'error.txt");
			for (i=1; i<=n; i++) fprintf(fe," %d ", (received[i]<0)!=codeword[i]);
			fprintf(fe,"\n");	
			fclose(fe);
		}
		*/
		if ((NDE>=10)&&((float)NDE/(J+1)>2.5*PR_ERR[1])) break;
	}
}

PR_ERR[0]=(float)NSE/(J)/n; //BER
PR_ERR[1]=(float)NDE/(J);  //FER
return 1;
}


/*************************************************/

int count_loops(int n,int r, int score[], int dd[])
{
// H is binary (0/1) parity-check matrix  for Tanner ase graph

int i, j, h, I, L, Num_cycles=0, ace, ace_min=NMAX;
int col_weights[CMAX*WMAX];
// Variables for EQ constructing
int N[CMAX*WMAX]; // # nodes at tree levels
int nodes, source, parent, child, branch, g=NMAX; 

//map dd into column weights of original matrix
h=0;  //column of original matrix
i=0; // # in the constant-weighh set
for (j=0; j<n; j++) 
{	
	col_weights[j]=dd[h];
	i++;
	if (i==dd[h]) {i=0; h++;}
}
H2graph(n,r);

for (source=0; source<r; source++)
{
   // Collection of cycles from source to source
   CNodes[0]=source;      // PNodes(1)=0;
   N[0]=0;  N[1]=1;
   PBranch[0]=0;          //previous branch
   nodes=0; 
   ACE[0]=0; L=-1;
   while (2*(L+1)<g)
   {
	   L++;
	   for (j=N[L]; j<N[L+1];j++) // over nodes of previous level
	   {
            parent=CNodes[j];               // parent node       
			if (parent>=source)    // TO AVOID IDENTICAL CYCLES
            for (h=0; h<Numbers[parent]; h++)
			{
                child=Graph[parent][h];      // child node
				//if (child<source) break;    // TO AVOID IDENTICAL CYCLES
                branch=Branches[parent][h];  // branch parent->child
                // check if not back
                if  (branch!=PBranch[j])
				{   //  Good. Add to the list 
                    nodes++;
                    CNodes[nodes]=child;
					ACE[nodes]=ACE[j]+col_weights[branch-1];
     				PBranch[nodes]=branch;
                    // Check for loops
					for(I=0; I<nodes; I++)
					{
                        if (CNodes[I]==CNodes[nodes])
						{
						    // Create new cycle
							if (2*(L+1)<g) 
							{
								ace_min=NMAX;  g=2*(L+1);
							}
							if (2*(L+1)==g)
							{
								ace=ACE[I]+ACE[nodes];
								if (ace<ace_min) 
								{
									Num_cycles=1;
									ace_min=ace;
								}
								else if (ace_min==ace) Num_cycles++;
							}
						} 
					}  // for I 
				}//
			}        // for h
	   } // for j
	   N[L+2]=nodes+1;   //     printf("\n Level:%d, eqs:%d",L,eqs);
   } // while or for L (depth
  
} // for source

score[0]=Num_cycles;
score[1]=ace_min-2*g;
return g;
} //count_cycles

/*************************************************/
/*************************************************/

int rand_base_matrix_g(int num_restr, int rw[], int rmin[], int rmax[], int r0, int r, int n, int msnd0, int msnd, 
	                   int DD0[], int DD[], int DDT[], int trials[], int score[], int init)
// generates random base LDPC matrix HB with degree distribution DD
// upper triangular (dual diagonal)  parity-check part
// DD(i) =  number of columns of weight i
// msnd   = max symbol node degree
{
	int i, I, j, jb, jb0, h, ii, ij, sum, t, g, nt, rt, gt, 
		N2, N3,  T2,T3, AT2, AT3,   // overall best
		w, w1, w2, row_weight[CMAX],row_weight0[CMAX], //col_weight[CMAX], 
		max_row_w0, max_row_w,  //increments
		a[CMAX], a1[CMAX], besta[CMAX], dd1[CMAX], dd2[CMAX], dd[CMAX],  
		//flag,row_weights[CMAX], 
		b0[CMAX], b[CMAX],
		girth[CMAX], ace[CMAX], num_cycles[CMAX];
	    // int  DDT[RMAX];
			
	// Check input
    
    for (i=0; i<r; i++) DDT[i]=0;
	sum=0;
	for (i=0; i<=r0; i++) sum+=DD0[i];
	if (sum != n-r0)
	{
		printf("sum(DD0) ~= n-r0"); return -1;
	}
	sum=0;
	for (i=0; i<=r-r0; i++) sum+=DD[i];
	if (sum != n-r)
	{
		printf("sum(DD) ~= n-r"); return -1;
	}
	// vector of symbol node weights
	
	// non-random parts
	for (i=0; i<r0-1; i++) {dd1[i]=2; dd2[i]=0;} 	dd1[r0-1]=3; dd2[r0-1]=0;
	for (i=r0; i<r-1; i++) {dd1[i]=0; dd2[i]=2;}   	dd2[r-1]=3;
	for (i=r0; i<r0+DD0[0]; i++) dd1[i]=0; 
	for (i=r; i<r+DD[0]; i++) dd2[i]=0; 
	
	// random parts
	sum=r0+DD0[0];
	for (i=1; i<=msnd0; i++)                   // msnd = max symbol node degree
	{
		for (h=0; h<DD0[i]; h++) dd1[sum+h]=i;             //    ???
		sum+=DD0[i];
	}
	sum=r+DD[0];
	for (i=1; i<=msnd; i++)                   // msnd = max symbol node degree
	{
		for (h=0; h<DD[i]; h++) dd2[sum+h]=i;             //    ???
		sum+=DD[i];
	}

	for (i=0; i<n; i++) 
	{ 
		dd[i]=dd1[i]+dd2[i]; 
		if (dd[i]<2) return -2;
	}
	for (i=0; i<n; i++) DDT[dd[i]]++;
	for (i=0; i<num_restr; i++)
		if ( (DDT[rw[i]]<rmin[i]) || (DDT[rw[i]]>rmax[i])) return -2;
	

	// non-random columns of weight 2
	for (i=0; i<r; i++) {row_weight0[i]=0; for (j=0; j<n; j++) T[i][j]=0;}
	T[0][0]=T[0][r0-1]=1;
	V[0][0]=0; V[0][1]=r0-1;
	row_weight0[0]=2;
	for (i=1; i<r0; i++) //column number
	{
		T[i][i-1]=T[i][i]=1; 
		V[i][0]=i-1; V[i][1]=i;
		row_weight0[i]=2; 
	}
    T[r0][r0]=T[r0][r-1]=1;
	V[r0][0]=r0; V[0][1]=r-1;
	row_weight0[i]=2;
	for (i=r0+1; i<r; i++) //column number
	{
		T[i][i-1]=T[i][i]=1; 
		V[i][0]=i-1; V[i][1]=i;
		row_weight0[i]=2; 
	}
	max_row_w=2;
   // special weight-3 columns
   i=(r0-1)/2; T[i][r0-1]=1; V[i][2]=i; row_weight0[i]=3; 
   i=(r+r0-1)/2; T[i][r-1]=1; V[i][2]=i; row_weight0[i]=3; 


   // random columns of higher weight
	girth[r0-1]=ace[r0-1]=0; 
	for (i=0; i<=n;i++) num_cycles[i]=1000000;
	
	for (t=0; t<trials[0]; t++)
	{
		max_row_w0=max_row_w=3;
		for (i=0; i<r; i++) row_weight[i]=row_weight0[i]; 
		// loop over codes
		AT2=AT3=0;
		girth[r0+DD0[0]-1]=6;                 // ???
		for (h=mini(r0+DD0[0],r+DD[0]); h<n; h++)
		{
			// loop over columns 

			if (girth[h-1]==0) break;                          // ???
			if (dd[h]>=4) gt=4; else gt=6;
			N2=N3=100000;
			
			w1=dd1[h]; w2=dd2[h]; w=dd[h];
			// generate random column
			//find low-weight rows
			// TRY trial[1] candidates and choose one with min N2,N3 contribution
			// priority to low-weight rows
			jb0=0; for (i=0; i<r0; i++) if (row_weight[i]<max_row_w) { b0[jb0]=i; jb0++;}; 
			jb=0; for (i=r0; i<r; i++)  if (row_weight[i]<max_row_w) { b[jb]=i-r0; jb++;}; 

			for (I=0; I<trials[1]; I++)
			{
				if (w1>0)
				{
				if (jb0>=w1)
				{
					j=0;
					while (j<w1)	j=randcomb(jb0, w1, a);
					for (ii=0; ii<w1; ii++) a[ii]=b0[a[ii]];
				}
				else randcomb(r0, w1, a);
				}
				
				if (h>=r)
				{
					if (w2>0)
					{
						if (jb>=w2)
						{
							j=0;
							while (j<w2)	j=randcomb(jb, w2, a1);
							for (ii=0; ii<w2; ii++) a1[ii]=b[a1[ii]];
						}
						else randcomb(r-r0, w2, a1);
					}
			   		for (j=w1; j<w; j++) a[j]=a1[j-w1]+r0; 
				}
				else
				{
					ii=w1;
					for (j=r0; j<r; j++) if (T[j][h]>0) { a[ii]=j; ii++;}
				}


				//	 compute correlations
				// rows:
				T2=T3=0; // contribution
				for (i=0; i<w-1; i++)
				{ //i-th row
					sum=0;
					for (j=i+1; j<w; j++)
					{ //with j-th row
						for (ij=0; ij < row_weight[a[i]]; ij++)
						sum+=T[a[j]][V[a[i]][ij]];
					}
					T2+=sum;
					T3+=binom(sum,2);
				}
				// columns. T2 is the same
				for (i=0; i<h-1; i++)
				{
					sum=0;
					for (j=0; j<w; j++) sum+=T[a[j]][i];
					T3+=binom(sum,3);
				}

				if ((T3<N3) || ((T3==N3) && (T2<N2)))  
				{
					for (i=0; i<w; i++) besta[i]=a[i]; 
					N3=T3; N2=T2; 
				}
				if ((N2==0)&&(N3==0)) break;
		   }  // for I (column-candidates)
		   // Check for # cycles: If bad then break
		   // Copy:
		   for (ii=0; ii<r; ii++) for (ij=0; ij<h; ij++) H[ii][ij]=T[ii][ij];
		   // Add column
		   for (ii=0; ii<r0; ii++) H[ii][h]=0;
		   if (h>=r) for (ii=r0; ii<r; ii++) H[ii][h]=0;
   		   for (ii=0; ii<w; ii++) H[besta[ii]][h]=1;
		   
		  // for (ii=0; ii<r; ii++) {for (ij=0; ij<n; ij++) 
			//   printf("%d ", H[ii][ij]); printf("\n%");}
		   
		   rt=h+r+1; 	nt=tanner(h+1,r);
		   for (ii=0; ii<rt; ii++) for (ij=0; ij<nt; ij++) H[ii][ij]=HH[ii][ij];
	       // printf("t=%d, h=%d, j=%d \n", t,h,j);
		   g=count_loops(nt, rt, score, dd);
		 //  printf("g=%d, ace=%d, num_cycles=%d \n",g, score[1], score[0]);
		   //flag=0;
		   if (g>=girth[h])
	       if (score[0]<=num_cycles[h])
			 {
				//for (i=0; i<w; i++) besta[i]=a[i];
				girth[h]=g; girth[h+1]=0;
			    ace[h]=score[1]; ace[h+1]=0;  
			    num_cycles[h]=score[0]; //num_cycles[h+1]=1000000;
				// grow matrix
				for (i=0; i<r; i++) T[i][h]=0;
				for (i=0; i<w; i++) 
				{
					T[besta[i]][h]=1;
					V[besta[i]][row_weight[besta[i]]]=h;
					row_weight[besta[i]]++;
					if (row_weight[besta[i]]>max_row_w) max_row_w=row_weight[besta[i]];
				}
				AT2+=T2; AT3+=T3;
				// printf("T=\n"); for (ii=0; ii<r; ii++) {for (ij=0; ij<=h; ij++) 
				// printf("%d ", T[ii][ij]); printf("\n%");}
				 //printf("h=%d, T2=%d, T3=%d \n", h, N2, N3);
				if ((h==n-1)&&(num_cycles[h]<num_cycles[n]))
				{
				 printf("t=%d,h=%d,girth=%d,ace=%d,num=%d,T2=%d, T3=%d, %d, %d \n", 
					 t,h, girth[r-1], ace[h], num_cycles[h], N2, N3, AT2, AT3);
				 if (girth[r-1]>6)
					 ii=0;
				 num_cycles[n]=num_cycles[h];
				 score[2]=N2; score[3]=N3;
				 for (i=0; i<r; i++) for (j=0; j<n; j++) H[i][j]=HM[i][j]=T[i][j];
				}
			}
			// else break;
		}  //for j (all columns of a given w
/*		
		 printf("T=\n"); for (ii=0; ii<r; ii++) {for (ij=0; ij<n; ij++) 
		 printf("%d ", T[ii][ij]); printf("\n%");}
*/
	}//  for code
 	//printf("T1=%d, T2=%d, T3=%d, t=%d \n", T1, T2, T3, t);
 score[0]=ace[n-1]; score[1]=num_cycles[n-1];
 

//---------------------------------------------
// Transform HM to mask-matrix: first components of rows/columns=2
// for (i=0; i<r; i++) for (j=0; j<n; j++) H[i][j]=HM[i][j]=T[i][j];
for (i=0; i<r-1; i++) {HM[i][i]=2; HM[i+1][i]=2;}
HM[r0][r0-1]=0; 
HM[0][r0-1]=2; 
HM[(r0-1)/2][r0-1]=3;
HM[r-1][r-1]=2;
HM[(r+r0-1)/2][r-1]=3;
HM[r0][r-1]=2;

for (i=r; i<n; i++)
{
	j=0; while ( (HM[j][i]==0) && (j<r)) j++; 
	
	if (j==r) 
		printf(" weight=0 ");
	HM[j][i]=2;
}


// modify the first w=3-column


return 0;
}
/*************************************************/


int bp_decod_lm(double soft[], int hard[], int n, int r, int cweight[], int rweight[], int maxiter)
{ 
  // BP decoding (Gallager)
  // y      is channel output 
  // wc(wr) is maximum column(row) weight
  // n,r,      is codelength and redundancy
  // maxiter is maximum number of iterations	
	int i,j, synd, iter=0, bs[NMAX],b;
	double A, s[KMAX], S, x, Eps=2.0e-5,  y[NMAX];
		
	while (iter < maxiter)
	{
		for (i=0; i<n; i++) hard[i]=(soft[i]<0);
		for (i=0; i<r; i++)
		{
			synd=0;
			for (j=0; j<rweight[i]; j++) synd^=(hard[V[i][j]]);
			if (synd) 
				{//printf(" iter, check =%d %d ",iter, i);
			break;}
		}
		if (!synd) return iter;// no errors detected
		if (iter==0) // prepare arrays
		for (i=0; i<n; i++) 
		{
			y[i]= soft[i]=maxd(mind(soft[i],  20.0),-20.0); 
			for (j=0; j<cweight[i]; j++) Z[i][j]=0; 
		}
		// Varable-node activation step
		for (i=0; i<n; i++)    
		{
			for (j=0; j<cweight[i]; j++) 
			{
			
				A=exp(soft[i] - Z[i][j]);
				B[i][j]=A<1;
				x=log(absd((A-1)/(A+1)));
				if (absd(x) < Eps) if (x<=0) x=-Eps; else x=Eps;
				Z[i][j]=x;
			}
		}
		// Check-node activation step
		for (i=0; i<r; i++)
		{  //row-wize products
			s[i]=0.0; bs[i]=0; 
			for (j=0; (j<rweight[i]); j++)
			{
				s[i]+=Z[ V[i][j] ] [ VI[i][j] ] ;
				bs[i]^=B[ V[i][j] ] [ VI[i][j] ];
			}
		}
    	for (i=0; i<n; i++)    
			for (j=0; j<cweight[i]; j++)   // over all ones in a column
			{
				A=exp( s[C[i][j]] - Z[i][j] );
				b=bs[C[i][j]]^B[i][j];
				Z[i][j]=(1-2*b)*log((1+A)/(1-A));
    			Z[i][j]=maxd(mind(Z[i][j], 19.07),-19.07);
			}

		// Output of iteration
		for (i=0; i<n; i++)
		{
			S=y[i]; 
			for (j=0; j<cweight[i]; j++) S+= Z[i][j];
			soft[i]=  S;
		}
		iter++;
	}  //while

	return -iter;
}   // bp
/*************************************************/
void hd2cv2(int b, int c, int M, int row_weights[], int col_weights[])
{
// Transfroms polynomial matrix HD
// to column locations C and row locations V
int i, j, h;
int r, n, c_ind, r_ind;
n=c*M;
r=b*M;

for (i=0; i<n; i++) col_weights[i]=0; 
for (i=0; i<r; i++) row_weights[i]=0;
for (i=0; i<b; i++)
	for (j=0; j<c; j++)
		if (HD[i][j]>=0)
		{
			for (h=0; h<M; h++)
			{
				c_ind=j*M+((HD[i][j]+h)%M);  // column number
				r_ind=i*M+h;                 // row number      
				C[c_ind][col_weights[c_ind]]=r_ind;
				V[r_ind][row_weights[r_ind]]=c_ind;
				VI[r_ind][row_weights[r_ind]]=col_weights[c_ind];
				col_weights[c_ind]++;
				row_weights[r_ind]++;
			}
		}
}

/*************************************************/
/*************************************************/
void cyclic_shift(int x[], int y[], int I, int M)
{
	int i;
	for (i=I; i<I+M; i++)  y[i%M]=x[i-I];
}
/*************************************************/
/*************************************************/
int random_codeword(int b, int c, int M)
{
	int i, j, h, n, k, r, p;
	int synd[RMAX], sumsynd[MMAX], 
		x[MMAX], y[MMAX], cword[NMAX];
//FILE *f;
	r=b*M;
	n=c*M;
	k=n-r;
	// find positive degree in (b-1)-th column
	p=0; 	while ((HD[p][b-1] <= 0)&&(p<b)) p++;
	if (p>=b) 
		return -1;
	
	//generate info part of codeword
	
	for (i=r; i<n; i++) cword[i]=rand()&1;
	/*
	fopen_s(&f, "message.txt", "wt");
	for (i=r; i<n; i++)  fprintf(f,"%d ", cword[i]); 
	fclose(f);
	*/
	// Compute partial syndrome
	for (i=0; i<r; i++) synd[i]=0; 
	for (h=0; h<M; h++) sumsynd[h]=0;
	for (i=0; i<b; i++)
	{
		for (j=b; j<c; j++)
		{
	       // circulate and add
			if (HD[i][j]>=0)
			{
				for (h=0; h<M; h++) y[h]=cword[j*M+h]; //j-th information block 
				cyclic_shift(y,x, M-HD[i][j],M);            //shift cyclic   
				for (h=0; h<M; h++) synd[i*M+h]=synd[i*M+h]^x[h]; //i-th syndrome component
			}
		}
		for (h=0; h<M; h++) sumsynd[h]^=synd[i*M+h];
	}
	// (b-1)st codeblock is created in x
	cyclic_shift(sumsynd, x, HD[p][b-1], M);

	for (h=0; h<M; h++) 
	{ 
		// (b-1)st block
		cword[(b-1)*M+h]=x[h];        
		// 0th block
		cword[h]=synd[h]; 
		if (HD[0][b-1]==0) cword[h]^=x[h];
		if (HD[0][b-1]>0) cword[h]^= sumsynd[h];  
	    // other blocks 
		for (i=1; i<b-1; i++) 
		{
			cword[i*M+h]=synd[i*M+h]^cword[(i-1)*M+h];
			if (HD[i][b-1]==0) cword[i*M+h]^=cword[(b-1)*M+h];
			if (HD[i][b-1]>0)  cword[i*M+h]^=sumsynd[h];
		}
	}
	//  CHECK for being codeword
	for (i=0; i<r; i++) synd[i]=0;
	for (i=0; i<b; i++)
		for (j=0; j<c; j++)
		{
			if (HD[i][j]>=0)
			{
				for (h=0; h<M; h++) y[h]=cword[j*M+h];
				cyclic_shift(y,x,M-HD[i][j],M);
				for (h=0; h<M; h++) synd[i*M+h]^=x[h];
			}
		}
	
	h=0;
	for (i=0; i<r; i++) h+=synd[i];
	for (i=0; i<n; i++) codeword[i]=cword[i];
		if (h!=0)
			h=1;

	return h;
	
}

/*************************************************/
/*************************************************/
int QCencoder(int b, int c, int M)
{
	int i, j, h, n, k, r, p;
	int synd[RMAX], sumsynd[MMAX], 
		x[MMAX], y[MMAX], cword[NMAX];
//FILE *f;
	r=b*M;
	n=c*M;
	k=n-r;
	// find positive degree in (b-1)-th column
	p=0; 	while ((HD[p][b-1] <= 0)&&(p<b)) p++;
	if (p>=b) 
		return -1;
	
	//generate info part of codeword
	
	for (i=r; i<n; i++) cword[i]=msg[i-r];
	/*
	fopen_s(&f, "message.txt", "wt");
	for (i=r; i<n; i++)  fprintf(f,"%d ", cword[i]); 
	fclose(f);
	*/
	// Compute partial syndrome
	for (i=0; i<r; i++) synd[i]=0; 
	for (h=0; h<M; h++) sumsynd[h]=0;
	for (i=0; i<b; i++)
	{
		for (j=b; j<c; j++)
		{
	       // circulate and add
			if (HD[i][j]>=0)
			{
				for (h=0; h<M; h++) y[h]=cword[j*M+h]; //j-th information block 
				cyclic_shift(y,x, M-HD[i][j],M);            //shift cyclic   
				for (h=0; h<M; h++) synd[i*M+h]=synd[i*M+h]^x[h]; //i-th syndrome component
			}
		}
		for (h=0; h<M; h++) sumsynd[h]^=synd[i*M+h];
	}
	// (b-1)st codeblock is created in x
	cyclic_shift(sumsynd, x, HD[p][b-1], M);

	for (h=0; h<M; h++) 
	{ 
		// (b-1)st block
		cword[(b-1)*M+h]=x[h];        
		// 0th block
		cword[h]=synd[h]; 
		if (HD[0][b-1]==0) cword[h]^=x[h];
		if (HD[0][b-1]>0) cword[h]^= sumsynd[h];  
	    // other blocks 
		for (i=1; i<b-1; i++) 
		{
			cword[i*M+h]=synd[i*M+h]^cword[(i-1)*M+h];
			if (HD[i][b-1]==0) cword[i*M+h]^=cword[(b-1)*M+h];
			if (HD[i][b-1]>0)  cword[i*M+h]^=sumsynd[h];
		}
	}
	//  CHECK for being codeword
	for (i=0; i<r; i++) synd[i]=0;
	for (i=0; i<b; i++)
		for (j=0; j<c; j++)
		{
			if (HD[i][j]>=0)
			{
				for (h=0; h<M; h++) y[h]=cword[j*M+h];
				cyclic_shift(y,x,M-HD[i][j],M);
				for (h=0; h<M; h++) synd[i*M+h]^=x[h];
			}
		}
	
	h=0;
	for (i=0; i<r; i++) h+=synd[i];
	for (i=0; i<n; i++) codeword[i]=cword[i];
		if (h!=0)
			h=1;

	return h;
	
}

/*************************************************/
/*************************************************/
 int syndrome(int row_weights[],int b, int c, int M)
 {
	 int i,j,h,r,s, synd[NMAX], x[MMAX], y[MMAX];

	 r=b*M;

	 //  CHECK for being codeword
	for (i=0; i<r; i++) synd[i]=0;
	for (i=0; i<b; i++)
		for (j=0; j<c; j++)
		{
			if (HD[i][j]>=0)
			{
				for (h=0; h<M; h++) y[h]=codeword[j*M+h];
				cyclic_shift(y,x,M-HD[i][j],M);
				for (h=0; h<M; h++) synd[i*M+h]^=x[h];
			}
		}
	
	h=0;
	for (i=0; i<r; i++) if (synd[i]) return i;

	for (i=0; i<r; i++)
		{
			s=0;
			for (j=0; j<row_weights[i]; j++) s^=(codeword[V[i][j]]);
			if (s) return -i;
		}


	return h;
 }