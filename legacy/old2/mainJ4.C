#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "graph.h"

//
main(void)
{
	FILE *f_report, *f_data, *f_data1; 	
	FILE *f_data_ext;

	int n, n0,             // length           
		g,                 // girth
	    r, r0, r01;             // check symbols number
	//	c,b,               // conv code parameters
	//	dmin,              // min distance 
	//	k,                 // inf symbols number 
	//	girth;             // girth
	
	int i, ic, h, s, slast, j, I, J, M, M0, M01, Mmax;//, TBL; 
	int  a[CMAX],  am[NMAX], wcol[CMAX], wcol1[CMAX],  pos[CMAX];
		
	int ok, new_code, ans, gtmp; 
	int Hflag, EQflag;
	int exact, flag1, Q, IP,  M1;
	
	// Variables for EQ constructing
	
	int eqs, nodes; 

	// prepare 


srand((unsigned)time(NULL)); //    srand(1);
Hflag=0; EQflag=0;	
gen_divtable();
while (1)
{
printf("\nChoose option:");
printf("\n 1: Read H from file inputh.txt   ");
printf("\n 2: Read Equations from file equations.txt   ");
printf("\n 3: Construct Equations and create file equations.txt");
printf("\n 4: Check voltages a from inputa.txt (degrees by columns) ");	   
printf("\n 5: Random search for a");	   

printf("\n 6: Random search with creating a bank of good codes");	   
printf("\n 7: Random extending by 1 column");	
printf("\n 8: Random extending by 1 row");	
printf("\n 9: Mix 2 codes into 1 ");	
printf("\n 0: Quit program");	   
printf("\n");	   
scanf("%d", &ans);

//*************************************************************//

switch (ans)    //main switch 
{ 
case 1:  // read H
    f_report=fopen("inputh.txt", "rt"); 
	// READ DATA
	fscanf(f_report,"%d", &r);
    fscanf(f_report,"%d", &n);
	fscanf(f_report,"%d", &g);
	fscanf(f_report,"%d", &M);
//	fscanf(f_report,"%d", &TBL); //TailBiting Length

	for (i=0; i<r; i++)
		for (j=0; j<n; j++)	{fscanf(f_report,"%d", &I); HM[i][j]=I; H[i][j]=(I>0);}
   	fclose(f_report);
I=1;
for (i=0; i<n; i++)
for (j=0; j<r; j++)
{
		T[j][i]=H[j][i];
		if (HM[j][i]>0)
		{	
			am[I]=HM[j][i];	I++;   // mask
		}
}
bin_matr_print(n, r, 0);
printf("\n g=%d, M=%d \n", g, M); 
	s=0; 
	for (i=0; i<n; i++)
	{
		wcol[i]=0;
		for (j=0; j<r; j++)	wcol[i]+=T[j][i];
		if (s<wcol[i]) s=wcol[i];
	}
    slast=wcol[n-1];
/*
	ans=1;
	if (s==2)
	{
		printf("Need Tanner graph? 1=YES\n");	   
		scanf("%d", &ans);
	}
	*/
	if (1)//(ans==1)
	{
		//input H, output HH
		r0=r; n0=n;
		r=n+r; //n=3*n;
		n=tanner(n,r0);
		for (i=0; i<r; i++)
			for (j=0; j<n; j++)	{H[i][j]=HH[i][j]; T[i][j]=H[i][j];}
		bin_matr_print(n, r, 0);
		printf("\n n=%d, r=%d \n", n, r);
	}
	
Hflag=1;


break;
//*************************************************************//
case 2:  // Read equations for H
if (Hflag==0){ printf("\n Matrix H is not known yet"); break;}
	
for (i=0; i<r; i++)
for (j=0; j<n; j++) T[i][j]=H[i][j];
bin_matr_print(n, r, 0);
printf("\n g=%d, M=%d \n", g, M); 

//f_report=fopen("equations.txt", "rt"); 
f_report=fopen("tree&cycles.txt", "rt"); 
// READ DATA
fscanf(f_report,"%d", &i);
if (n!=i) {printf("\n Column number does not match"); break;}
fscanf(f_report,"%d", &i);
if (g!=i) {printf("\n Girth does not match"); break;}
fscanf(f_report,"%d", &eqs);
fscanf(f_report,"%d", &nodes);

//read EQuations

//for (j=0; j<eqs; j++) 
//	for (i=0; i<n; i++) {fscanf(f_report,"%d", &I); EQ[j][i]=I;}
//read Length
//for (j=0; j<eqs; j++)   {fscanf(f_report,"%d", &I); Length[j]=I;};
//read Tree
for (j=0; j<nodes; j++) 
	for (i=0; i<2; i++) {fscanf(f_report,"%d", &I); TREE[j][i]=I;};
//read CYCLES
for (j=0; j<eqs; j++) 
	for (i=0; i<2; i++) {fscanf(f_report,"%d", &I); CYCLES[j][i]=I;};
fclose(f_report);


printf("eqs=%d, nodes=%d\n", eqs, nodes);

EQflag=1;
break;
//*************************************************************//
case 3:  // Construct equations for H

if (Hflag==0){ printf("\n Matrix H is not known yet"); break;}

for (i=0; i<r; i++)
for (j=0; j<n; j++) T[i][j]=H[i][j];
bin_matr_print(n, r, 0);
ok=eq_generation(n,r, g, &eqs, &nodes);
if (ok)
printf("ok: eqs=%d, nodes=%d\n", eqs, nodes);
else 
{printf("Problem!"); break;}
// Write to file
//f_report=fopen("equations.txt", "wt"); 
f_report=fopen("tree&cycles.txt", "wt"); 
fprintf(f_report," %d %d \n",n,g);
fprintf(f_report," %d %d \n", eqs, nodes);
//print EQ
/*
for (j=0; j<eqs; j++) 
{
	for (i=0; i<n; i++) fprintf( f_report," %d ", EQ[j][i]);
	fprintf( f_report," \n ");
}
//print Length
for (j=0; j<eqs; j++) 
	fprintf( f_report," %d \n", Length[j]);
*/
//print Tree
for (j=0; j<nodes; j++) 
{
	for (i=0; i<2; i++) fprintf( f_report," %d ", TREE[j][i]);
	fprintf( f_report," \n ");
}
//print CYCLES
for (j=0; j<eqs; j++) 
{
	for (i=0; i<2; i++) fprintf( f_report," %d ", CYCLES[j][i]);
	fprintf( f_report," \n ");
}

fclose(f_report);

EQflag=1;
break;  //case 3
//*************************************************************//
case 4: // Check voltages a from inputa.txt (degrees row by row ");	   
if (Hflag==0){ printf("\n Matrix H is not known yet"); break;}
//if (EQflag==0){ printf("\n Equations are not created yet"); break;}
// Read data
f_report=fopen("inputa.txt", "rt"); 
for (i=1; i<=n; i++){ fscanf(f_report,"%d", &I); a[i]=I;}
fclose(f_report);
//printf("\nType report? 1==YES:"); scanf("%d",&ans);

//I=1; 
gtmp=g;
while (gtmp>2)
{
	ok=eq_generation(n,r, gtmp, &eqs, &nodes);
	if (ok)
	{
		//if (I==1)
		//{
		//	for (i=0; i<eqs; i++)
		//		for (j=1; j<=g/2;j++) FLAGS_ARRAY[i][j]=0;
		//	J=eqs;
		//}
		new_code=check_eq(a, eqs, nodes);
		if (new_code<0) 
			printf(" gfree< %d, wrong eq_number=%d \n", gtmp, -new_code);       
		else 
		{
			new_code=min_module(a, eqs, nodes, &Mmax);
			printf(" g>=%d for  M= %d and all M>Mmax=%d\n", gtmp, new_code, Mmax); 
			break;
		}
		//for (i=0; i<eqs; i++) FLAGS_ARRAY[i][I]=FLAGS[i];
		gtmp-=2; //I++;
	}
}

if (1)// (ans==1)
{
	f_report=fopen("inputa.txt", "at"); 
	fprintf(f_report,"\ng=%d for M= %d and all M>Mmax=%d\n", gtmp, new_code, Mmax);
	I=0;
	for (i=new_code; i<=Mmax; i++)
	{
		if (FLAGS[i]==0) 
		{
			fprintf(f_report,"%d ", i);
			printf("%d ", i);
			I++ ;
		    if (I>20) 
			{
				fprintf(f_report,"\n"); 
				printf("\n"); 
				I=0;
			} 
		//for (j=1; j<I; j++) fprintf(f_report,"%d ", FLAGS_ARRAY[i][j]);
		}
	}
	fprintf(f_report,"...");
	printf("..."); 
	fclose(f_report);
}
	

break; //case 4
//*************************************************************//
case 5: // start random search            
// Search graph with a given base graph 
if (Hflag==0){ printf("\n Matrix H is not known yet\n"); break;}
if (EQflag==0){ printf("\n Equations are not created yet\n"); break;}

for (i=0; i<r; i++)
for (j=0; j<n; j++) T[i][j]=H[i][j];
bin_matr_print(n, r, 0);
printf("\n eqs=%d, nodes=%d \n", eqs, nodes); 
for (i=1; i<=n; i++) a[i]=0;
//M1=M;
M1=10000000;

I=0;
for (;;)
{
	I++;
	if(I%1000000==0) printf("\n %d tested",I);

	for (i=1; i<=n; i++) 
	{                 
		//if (am[i]==1 || s==2)  a[i]=rand()%M;
		if (am[i]==1 )  a[i]=rand()%M;
		else a[i]=0;
	}
	

	new_code=check_eq(a, eqs, nodes);
	if (new_code>0)
	{
		if (new_code<M1)
		{
			new_code=min_module(a, eqs, nodes, &Mmax);
			if (new_code<M1) 
				{
					M1=new_code;
	   //		if (1) //(M1<M) 
		//		{
				//	M=M1;
			printf(" g=%d for M= %d\n", g, M1);
			f_report=fopen("in_out.txt", "at"); 
			fprintf(f_report,"\ng=%d,M=%d, a=",g,M1);	
			for (i=1; i<=n; i++) fprintf(f_report," %d ", a[i]);
			fclose(f_report);
			}
			
		}
	}
} //for I  (random search loop

break; // case 5
//*************************************************************//
case 6:  //Random search with creating a bank of good codes
// Search graph with a given base graph 
if (Hflag==0){ printf("\n Matrix H is not known yet\n"); break;}
if (EQflag==0){ printf("\n Equations are not created yet\n"); break;}

for (i=0; i<r; i++)
for (j=0; j<n; j++) T[i][j]=H[i][j];
bin_matr_print(n, r, 0);
printf("\n eqs=%d, nodes=%d \n", eqs, nodes); 
for (i=1; i<=n; i++) a[i]=0;
M1=M;
exact=0;
printf("\nExact M? 1=YES:"); scanf("%d",&exact);
if (exact) 
{
	printf("\nInput exact M="); scanf("%d",&exact);
}


I=0;
for (;;)
{
	I++;
	if(I%1000000==0) printf("\n %d tested",I);

	for (i=1; i<=n; i++) 
	{                 
		//if (am[i]==1 || s==2)  a[i]=rand()%M;
		if (am[i]==1 )  a[i]=rand()%M;
	//	else a[i]=0;
	}
	
	new_code=check_eq(a, eqs, nodes);
	if (new_code>0)	
	{
		IP++; Mmax=M;
		new_code=min_module(a, eqs, nodes, &Mmax);
		if (new_code<=M)
		{

			if ((FLAGS[exact]==0)&&(exact>0))
			{
				f_data=fopen("data.txt", "at"); 
				fprintf(f_data,"%d ", new_code);
				for (i=1; i<=n; i++) fprintf(f_data," %d ", a[i]);
				fprintf(f_data,"\n");	
				fclose(f_data);
			}
		}
		
		// IF BEST
		if (new_code<M1)	
		{
			f_report=fopen("in_out.txt", "at"); 
			fprintf(f_report,"\ng=%d for M= %d\n", g, new_code);
			for (i=1; i<=n; i++) fprintf(f_report," %d ", a[i]);
			fclose(f_report);
			M1=new_code;
			printf(" g=%d for M= %d\n", g, new_code);
		}
	}
}//2
	
break; // case 6

//*************************************************************//
case 7: // extension of existing H by extra column
if (Hflag==0){ printf("\n Matrix H is not known yet"); break;}
if (EQflag==0){ printf("\n Equations are not created yet"); break;}

f_data=fopen("data.txt", "rt"); 

f_data_ext=fopen("data_ext.txt", "wt"); 
fclose(f_data_ext);


exact=0;
printf("\nExact M? 1=YES:"); scanf("%d",&exact);
if (exact) 
{
	printf("\nInput exact M="); scanf("%d",&exact);
}



I=0; IP=0; flag1=0; J=0;
Q=M; M1=M;
while (1)
{
	// read candidate
	fscanf(f_data,"%d", &M0); 
	if (M0<0)  //refresh
	{
		fclose(f_data);
		f_data=fopen("data.txt", "rt"); 
		fscanf(f_data,"%d", &M0); 
	}
	if (M0<M)
	{
		for (j=1; j<=n-slast; j++) 
		{
			fscanf(f_data,"%d", &J); a[j]=J;
		}
		// extend it
		for (j=0; j<mini(M0,10); j++)
		{
			for (i=n-slast+1; i<=n; i++) 
				if (am[i]==1 || slast==2)  
					a[i]=rand()%M;  
				else a[i]=0; 
			new_code=check_eq(a, eqs, nodes);
			if (new_code>0)
			{
				new_code=min_module(a, eqs, nodes, &Mmax);
				IP++;
				if (new_code<M)	
				{
				if ((FLAGS[exact]==0)&&(exact>0))
					{
						J++;
						f_data_ext=fopen("data_ext.txt", "at"); 
						fprintf(f_data_ext," %d ", new_code);
						for (i=1; i<=n; i++) fprintf(f_data_ext," %d ", a[i]);
						fprintf(f_data_ext,"\n");	
						fclose(f_data_ext);
					}
				}
				if (new_code<M1)	
				{
					printf(" g=%d for M= %d from M=%d \n", g,new_code,M0);  M1=new_code;
				    f_report=fopen("in_out.txt", "at"); 
					fprintf(f_report,"\ng=%d for M= %d\n", g, M1);
					for (i=1; i<=n; i++) fprintf(f_report," %d ", a[i]);
					fclose(f_report);
				}
			} //if new code
		} // iterations for fixed input
	}  // M0<M
}  // while 1
fclose(f_data);
break; // case 7
case 8: // extension of existing H by extra row
if (Hflag==0){ printf("\n Matrix H is not known yet"); break;}
if (EQflag==0){ printf("\n Equations are not created yet"); break;}

f_data=fopen("data.txt", "rt"); 

f_data_ext=fopen("data_ext.txt", "wt"); 
fclose(f_data_ext);

// Column weights without last row:
for (i=0; i<n0; i++)	wcol[i]-= (HM[r0-1][i]>0);
		

I=0; IP=0; flag1=0; J=0;
Q=M; M1=M;
while (1)
{
	// read candidate
	fscanf(f_data,"%d", &M0); 
	if (M0<0)  //refresh
	{
		fclose(f_data);
		f_data=fopen("data.txt", "rt"); 
		fscanf(f_data,"%d", &M0); 
		printf("refresh");
	}

	if (M0<M)
	{
		i=0; ic=0;
		for (j=0; j<n0; j++)
		{   // over columns of base matrix
			for (h=1; h<=wcol[j]; h++)
			{
				i++;
				fscanf(f_data,"%d", &J); a[i]=J;
			}
			if (HM[r0-1][j]>0)
			{
				i++; a[i]=0;
				if (HM[r0-1][j]==1) {pos[ic]=i; ic++;}
			}
		}
 		// extend it
		for (j=0; j<M0; j++)
		{
			for (i=0; i<ic; i++) a[pos[i]]=rand()%M;  
				 
			new_code=check_eq(a, eqs, nodes);
			if (new_code>0)
			{
				new_code=min_module(a, eqs, nodes, &Mmax);
				IP++;
				if (new_code<M)	
				{
					J++;
					f_data_ext=fopen("data_ext.txt", "at"); 
					fprintf(f_data_ext," %d ", new_code);
					for (i=1; i<=n; i++) fprintf(f_data_ext," %d ", a[i]);
					fprintf(f_data_ext,"\n");	
					fclose(f_data_ext);
				}
				if (new_code<M1)	
				{
					printf(" g=%d for M= %d from M=%d \n", g,new_code,M0);  M1=new_code;
				    f_report=fopen("in_out.txt", "at"); 
					fprintf(f_report,"\ng=%d for M= %d\n", g, M1);
					for (i=1; i<=n; i++) fprintf(f_report," %d ", a[i]);
					fclose(f_report);
				}
			} //if new code
		} // iterations for fixed input
	}  // M0<M
}  // while 1
fclose(f_data);
break; // case 8

case 9: // Mix two existing H to one matrix
if (Hflag==0){ printf("\n Matrix H is not known yet"); break;}
if (EQflag==0){ printf("\n Equations are not created yet"); break;}

printf("\nInput number of rows in first sub-matrix"); scanf("%d", &r01);


f_data=fopen("data.txt", "rt"); 
f_data1=fopen("data1.txt", "rt"); 

f_data_ext=fopen("data_ext.txt", "wt"); 
fclose(f_data_ext);

// Column weights 
for (i=0; i<n0; i++)
{
	wcol[i]=0;
	for (j=0; j< r01; j++) if (HM[j][i]>0) wcol[i]++;
	wcol1[i]=0;
	for (j=r0/2; j< (r0-r01); j++) if (HM[j][i]>0) wcol1[i]++;
}
		

I=0; IP=0; flag1=0; J=0;
Q=M; M1=M;
while (1)
{
	for (i=1; i<=n; i++) a[i]=0;
	// read candidate
	fscanf(f_data,"%d", &M0); 
	if (M0<0)  //refresh
	{
		fclose(f_data);
		f_data=fopen("data.txt", "rt"); 
		fscanf(f_data,"%d", &M0); 
		printf("0");
	}
	// Fill upper part of H
	i=0; 
	for (j=0; j<n0; j++)
	{   // over columns of base matrix 
		for (h=1; h<=wcol[j]; h++)
		{
			i++;
			fscanf(f_data,"%d", &J); a[i]=J;
		}
		i+=wcol1[j];
	}
	if (M0<M)
	{
		while (1)
		{
			fscanf(f_data1,"%d", &M01); 
			if (M01<0)  //refresh
			{
				fclose(f_data1);
				f_data1=fopen("data1.txt", "rt"); 
				printf("1");
				break;
			}
			//compose one matrix from two
			i=0; 
			for (j=0; j<n0; j++)
			{   // over columns of base matrix
				i+=wcol[j];
				for (h=1; h<=wcol1[j]; h++)
				{
					i++;
					fscanf(f_data1,"%d", &J); a[i]=J;
				}
			}
			// check combined code	 
			new_code=check_eq(a, eqs, nodes);
			if (new_code>0)
			{
				new_code=min_module(a, eqs, nodes, &Mmax);
				IP++;
				if (new_code<M)	
				{
					J++;
					f_data_ext=fopen("data_ext.txt", "at"); 
					fprintf(f_data_ext," %d ", new_code);
					for (i=1; i<=n; i++) fprintf(f_data_ext," %d ", a[i]);
					fprintf(f_data_ext,"\n");	
					fclose(f_data_ext);
				}
				if (new_code<M1)	
				{
					printf(" g=%d for M= %d from M=%d \n", g,new_code,M0);  M1=new_code;
				    f_report=fopen("in_out.txt", "at"); 
					fprintf(f_report,"\ng=%d for M= %d\n", g, M1);
					for (i=1; i<=n; i++) fprintf(f_report," %d ", a[i]);
					fclose(f_report);
				}
			} //if new code
		} // iterations for given first code
	}  // M0<M
}  // while 1
fclose(f_data);
break; // case 9

case 0: goto finish;  

//*************************************************************//
// case 
} //end while (1)

}  // end main
finish: ;
}



 
	
