#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "graph.h"
/*************************************************/
int mini(int a, int b)
{
	if (a<b) return a; else return b;
}
/*************************************************/
int maxi(int a, int b)
{
	if (a<b) return b; else return a;
}
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

if (f) f_matr=fopen("matr.dat", "at");
for (i=0; i<k; i++)
{ 
	//word2bin(g[i],n,row);
	for (j=0; j<n; j++) 
	{
		printf("%2d",T[i][j]); 
		if (f) fprintf(f_matr,"%2d",T[i][j]); 
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
		if (H[i][j]==1)
		// check the column for more ones
		for (h=0; h<r; h++)
		if ((h!=i) && (H[h][j]==1))
		{
			Graph[i][Numbers[i]]=h;
			Branches[i][Numbers[i]]=j+1;
			Numbers[i]++;
		}
	}
}
} //inc2graph

/*************************************************/
int dijkstra(int n)
{
/*
n is number of nodes
output: D = length of shortest cycle
Globals:
 Numbers  is the number of branches from i-th vertex
 Graph    is the list of nodes connected with i-th vertex
*/
int i, j, J,
    source, current, 
	num_visited,
	newnode,
	new_metric,
	d,
    D;
int distance[2*NMAX],
	visited[2*NMAX],
	previous[2*NMAX];

D=n;

for (source=0; source<n; source++)
{
	//Search for shortest cycle from source to source
    // Initialization
	for (i=0; i<n; i++)
	{
		distance[i]=n<<1;
		visited[i]=0;
		previous[i]=0;
	}
	distance[source]=0;
    num_visited=0;
    current=source;
	while (num_visited<n)
	{
		visited[current]=1;
        num_visited++;
		for (i=0; i<Numbers[current]; i++)
		{
            newnode=Graph[current][i];
            if (newnode!=previous[current])
			{         // not back
				new_metric=distance[current]+1;
				if (distance[newnode]<(n<<1)) // means that we alredy 
				{                             // explored the node
                   d=distance[newnode]+new_metric; // loop weight
                   if (d<D) 
				   {  // update record
					   D=d; 
					  // recover the path
					  // back from current to source
					   J=new_metric-1;
					   for (j=current; j!=source;j=previous[j])
					   {
						   path[J]=j; J--;
					   }
					   J=new_metric;
					   for (j=newnode; j!=source;j=previous[j])
					   {
						   path[J]=j; J++;
					   }
				   }
               	}
                if ((distance[newnode]>new_metric) && (newnode!=source))
				{
                    // update
                    distance[newnode]=new_metric;
                    previous[newnode]=current;
				}
  			}
		}
		d=D;
		for (i=0; i<n; i++)
		if ((!visited[i])&& (distance[i]<d)) 
		{
			d=distance[i];
			current=i;
		}
	}
}
return D;
} //dijkstra

/*************************************************/
/*
int graph_girth(int n, int r, int b)
{
	int t;
//	int i,j;
	// Computes girth from H
    if (b>2)
	{
		t=h2incidence(n, &r);
		inc2graph (t, n+r);
	    return dijkstra(n+r);
	}
	else
	{
		H2graph(n,r);
		return 2*dijkstra(r);
	}

} //graph girth
*/
/*************************************************/
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
	//for (j=0; j<c; j++) HTL[i][j]=(_int64)((_int64)1)<<(_int64)(HD[i][j]%L);
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
			//(int)HTL[i][j]&1;	HTL[i][j]>>=1;
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
int new_eq[CMAX];       // 

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
								for (J=0; J<n; J++) printf(" %d", new_eq[J]); printf("\n");
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
   printf(" source, eqs, nodes:");
    printf(" %d %d %d\n",source, eqs, nodes);
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

for (i=0; i<MAX_EQ; i++) FLAGS[i]=0;

// Computing weights of tree nodes and cycles
j=0; // # already computed tree nodes
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
int nextsorted(int b[], int n, int M )
{
// generates all M-ary vectors in ascending order
int i,j;
// all indices like in MATLAB are from 1

i=1; n--;
b[n+1]=M;
while ((b[i]==b[i+1]) && (i<=n))	i++;
	
if (i<=n) 
{
	b[i]++;
	for (j=1; j<i; j++) b[j]=0;
	return 0;
}
else
return 1; // finished
}// end next sorted    
/*************************************************/
int nextnonsorted(int b[], int n, int M )
{
// generates all M-ary vectors 
int i,j;
// all indices like in MATLAB are from 1

i=1; //n--;
b[n+1]=M;
while ((b[i]==M) && (i<=n))	i++;
	
if (i<=n) 
{
	b[i]++;
	for (j=1; j<i; j++) b[j]=0;
	return 0;
}
else
return 1; // finished
}// end next non-sorted    
    

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
