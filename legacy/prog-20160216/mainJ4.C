#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "graph.h"
#include "commons_portable.h"

//
int main(void)
{
	FILE *f_report, *f_data;
	FILE *f_data_ext;

	int n,      //  max length of base code
		nb,ns,     //  length of base code
		nt, nts,     //  length of tanner code
		ninp,   //  number of columns to read from input file
		g, gs, gtmp,//  girth
	    r,      //  base code check symbols number
		rt, rts;     //  tanner code check symbols number

	int i,  j, I, J, M, M0, Mmax;
	int  A[NMAX], a[NMAX], as[NMAX], //mask[NMAX],
		 P[CMAX], am[NMAX], wcol[CMAX], cumwcol[CMAX];
	int map[NMAX]; // mapping form initial code to base code
	int ok, new_code, ans, t;
	int Hflag, flag_continue;
	int exact, M1;

	int N, min_codes, max_codes, eqss, nodess;

	// Variables for EQ constructing
	int eqs, nodes;
	char letter;
	char name[]="data00.txt",*mwt="wt", *mat="at", *mrt="rt";
	char nextname[]="data00.txt";


// prepare 
srand((unsigned)time(NULL)); //    srand(1);
Hflag=0; 	
gen_divtable();
dialog:
while (1)
{
printf("\nChoose option:");
printf("\n 1: Read H from file inputh.txt and construct equations ");
printf("\n 2: Check voltages a from inputa.txt (degrees by columns) ");	   
printf("\n 3: Random search with creating a bank of good codes");	   
printf("\n 4: Random extending 'data.txt' by 1 or more columns to 'data_ext.txt' ");
printf("\n 5: Filter the set in 'data.txt' to 'data_f.txt' ");
printf("\n 6: Creating a code according to 'scenario.txt' ");
printf("\n 7: Average girth/spectrum ");
printf("\n 0: Quit program");	   
printf("\n");	   
scanf("%d", &ans);

//*************************************************************//

switch (ans)    //main switch 
{ 
case 1:  // read H
	f_report = fopen_or_die("input_h.txt", "rt", "The file 'input_h.txt' not found");
	// READ DATA
	fscanf_1(f_report, "%d", r);  // number of rows
	skip_until_newline(f_report);
    fscanf_1(f_report,"%d", n);  // total number of columns 
    skip_until_newline(f_report);
	fscanf_1(f_report,"%d", g);  // target girth
	skip_until_newline(f_report);
	fscanf_1(f_report,"%d", M);  // tailbiting length
	skip_until_newline(f_report);
	fscanf_1(f_report,"%d", exact);  exact*=M;// ecact tailbiting length
	skip_until_newline(f_report);
	fscanf_1(f_report,"%d", nb);  // length of base code 
	skip_until_newline(f_report);
	for (i=1; i<=nb; i++) 
	{
		fscanf_1(f_report,"%d", I); P[i]=I; //  column
	}
	P[0]=0;
	skip_until_newline(f_report);
	fscanf_1(f_report,"%d", ninp);  // number of input columns (for analysis and extending)
	skip_until_newline(f_report);
	// scan full HM 
	for (i=0; i<r; i++) for (j=0; j<n; j++)	
	{
		fscanf_1(f_report,"%d", I); 
	    HM[i][j]=I;    // Full matriix with 0,1,2  
		T[i][j]=(I>0); // binary copy
	}
    fclose(f_report);

   	// column weights (for all, not only used columns)  
	cumwcol[0]=0;
	for (i=0; i<n; i++)
	{
		wcol[i]=0;	for (j=0; j<r; j++)	wcol[i]+=T[j][i];
		cumwcol[i+1]=cumwcol[i]+wcol[i];
	}


    // obtain base matrix of proper size for constructing equations
	I=0;
	for (i=1; i<=nb; i++)
	{
		for (j=0; j<r; j++)	H[j][I]=T[j][P[i]-1]; 
		I++;
	};
	// mapping for nonzero positions
	I=0;
	for (i=1; i<=nb; i++)
	{
		for (j=cumwcol[P[i]-1]+1; j<=cumwcol[P[i]]; j++)
		{
			I++; map[I]=j;
		}
	};


	// for (i=nf-1; i<nl; i++) for (j=0; j<r; j++)  H[j][i-nf+1]=T[j][i];
	if (nb < 50) bin_matr_print(nb, r, 0);
	printf("\n n=%d, r=%d, g=%d, M=%d \n",nb, r, g, M); 
	    
	//   Tanner graph 
	rt=nb+r; 	nt=tanner(nb,r);  	                             //input H, output HH
	for (i=0; i<rt; i++) for (j=0; j<nt; j++) H[i][j]=HH[i][j];  // Tanner graph PCM 
	if (nt<50) bin_matr_print(nt, rt, 0);
	printf("\n Tanner graph size: nt=%d, rt=%d \n", nt, rt);
	
	// Generating equations
	printf("Generating equations  \n");
	ok=0;
	while(!ok){ok=eq_generation(nt,rt, g, &eqs, &nodes); g-=2;}; g+=2;
	if (ok)  printf("ok: g=%d, eqs=%d, nodes=%d\n",g, eqs, nodes);
	else  {printf("Problem!"); break;}
	
	Hflag=1; 
break;
//*************************************************************//
case 2: // Check voltages a from inputa.txt (degrees row by row ");	   
if (Hflag==0){ printf("\n Matrix H is not known. Use Option 1 to load "); break;}
// Read data from inputa.txt
f_report=fopen_or_die("input_a.txt", "rt", " File 'input_a.txt' not found"); 
if (P[nb]>ninp){ printf(" Number of input columns ninp too small\n"); break;} 

for (i=1; i<=cumwcol[ninp]; i++){ fscanf_1(f_report,"%d", I); A[i]=I;}
// select subsequence for analysis
for (i=1; i<=nt; i++) a[i]=A[map[i]];

// for (i=cumwcol[nf-1]+1; i<=cumwcol[nl]; i++) a[i-cumwcol[nf-1]]=A[i];

//I=1; 
gtmp=g; Mmax=100000;
ans=1;
while (ans!=0)
{
//	new_code=check_eq(a, eqs, nodes);
	ok=eq_generation(nt,rt, gtmp, &eqs, &nodes);
	if (ok)
	{
		new_code=check_eq(a, eqs, nodes);
		if (new_code<0) 
		{
			printf(" gfree< %d, wrong eq_number=%d \n", gtmp, -new_code);
			if ((ans==1)||(ans==-2)) ans=-2;  // move down
			else ans=0;
		}
		else 
		{
			new_code=min_module(a, eqs, nodes, &Mmax);
			printf(" g>=%d for  M= %d and all M>Mmax=%d\n", gtmp, new_code, Mmax); 
			fprintf(f_report,"\ng=%d for M= %d and all M>Mmax=%d\n", gtmp, new_code, Mmax);
			I=0;
			if (new_code>0){
			for (i=new_code; i<=Mmax; i++)
			{
				if (FLAGS[i]==0) 
				{
					fprintf(f_report,"%d ", i); printf("%d ", i);
					I++ ; if (I>20) {fprintf(f_report,"\n"); printf("\n"); I=0;}	 
				}
			}
			fprintf(f_report,"..."); printf("...\n");}
			if ((ans==1)||(ans==2)) ans=2; //move app
			else ans=0;
			//break;
		}
		gtmp+=ans; 
	}
	else break;
}
fclose(f_report);
Hflag=0;
break; //case 2
//*************************************************************//
case 3: // start random search 
if (Hflag==0){ printf("\n Matrix H is not loaded, use option 1\n"); break;}
for (i=0; i<nb; i++) 	if ((P[i]+1) != (P[i+1]))
{
   printf(" Wrong base matrix. It should consist of nb first columns. Edit configuration. "); 
   goto dialog;
}

printf("New search. Data.txt will be rewritten. Press any key to continue");
ans=get_character();
f_data=fopen_or_die("data.txt", "wt", "Problem with opening 'data.txt'");
fclose(f_data);

// Generate mask
I=1;
for (i=0; i<nb; i++) for (j=0; j<r; j++)
	   if (HM[j][i]>0) 
		{	
			am[I]=HM[j][i];	I++;   // mask
		}	
	
for (i=1; i<=nt; i++) a[i]=0; M1=M;

printf("\n Press 'x' to stop ");

I=0; J=0;
for (;;)
{
	I++;  if(I%1000000==0) printf("\n %d tested",I);
	for (i=1; i<=nt; i++) 
	{                 
		if (am[i]==1 )  a[i]=next_random_int(0, M); 	else a[i]=0;
	}
//	a[5]=16; a[6]=4; a[8]=21; a[9]=18;  a[11]=19; a[12]=20;
	new_code=check_eq(a, eqs, nodes);
	if (new_code>0)	
	{
		Mmax=M;
		new_code=min_module(a, eqs, nodes, &Mmax);
		if (new_code<=M)
		{
			if ((exact==0)||((FLAGS[exact]==0)&&(exact>0)))
			{
				f_data = fopen_or_die("data.txt", "at", "Problem with opening 'data.txt");
				fprintf(f_data,"%d ", new_code);
				for (i=1; i<=nt; i++) fprintf(f_data," %d ", a[i]);
				fprintf(f_data,"\n");	
				fclose(f_data); 
				J++;  if (J%1000==0) printf ("%d codes written to 'data.txt'\n",J);
				if (new_code<M1)	
				{		
					f_report = fopen_or_die("in_out.txt", "at", "Problem with opening 'in_out.txt'");
					fprintf(f_report,"\ng=%d for M= %d\n", g, new_code);
					for (i=1; i<=nt; i++) fprintf(f_report," %d ", a[i]);
					fclose(f_report);
					M1=new_code;
					printf(" g=%d for M= %d\n", g, new_code);
				}
			}
		}
	}
	if (is_keyboard_hit())
	{
		letter=get_character();
		if (letter=='x') 
		{
			f_data_ext = fopen_or_die("data.txt", "at", "Problem with opening 'data.txt'"); 
			fprintf(f_data_ext," %d\n ", -1);
			fclose(f_data_ext);
			break;
		}
	}
}//2
	
break; // case 3
//*************************************************************//
case 4: // extension of existing H by extra column (s)
// Check inputs
if (Hflag==0){printf("\n Matrix H is not loaded yet, use option 1\n"); break;}

for (i=0; i<nb; i++) 	if ((P[i]+1) != (P[i+1]))
{
   printf("Base should consist of nb first columns. Edit configuration. "); 
   goto dialog;
}
	
printf("\nExtending %d colums to  %d columns\n", ninp, nb);
if (ninp>=nb){printf("\nCheck config: nothing to extend"); break;}

// Generate mask for last branch(es)
I=cumwcol[ninp]+1; 
for (i=ninp;i<nb; i++)
for (j=0; j<r; j++)
if (HM[j][i]>0)
{	
	am[I]=HM[j][i];	I++;   // mask
}	
	
for (i=1; i<=nt; i++) a[i]=0;
M1=M;

f_data=fopen_or_die("data.txt", "rt", "Problem with opening 'data.txt");

f_data_ext=fopen_or_die("data_ext.txt", "wt", "Problem with opening 'data_ext.txt");
fclose(f_data_ext);

I=0; 
M1=M;
printf("\n Press 'x' to stop ");

while (1)
{
		// read candidate
	fscanf_1(f_data,"%d", M0); 
	if (M0<0)  //refresh
	{
		fclose(f_data);
		f_data = fopen_or_die("data.txt", "rt", "Problem with opening 'data.txt");
		fscanf_1(f_data,"%d", M0); 
	}

	if (M0<=M)
	{
		for (j=1; j<=cumwcol[ninp]; j++) 
		{
			fscanf_1(f_data,"%d", J); a[j]=J;
		}
		// extend it
		for (j=0; j<mini(M0,10); j++)
		{
			for (i=cumwcol[ninp]+1; i<=nt; i++) 
				if (am[i]==1) a[i]=next_random_int(0, M);  else a[i]=0; 
			new_code=check_eq(a, eqs, nodes);
			if (new_code>0)
			{
				new_code=min_module(a, eqs, nodes, &Mmax);
				if (new_code<=M)	
				{
				    if ((exact==0)||((FLAGS[exact]==0)&&(exact>0)))
					{
	     				f_data_ext=fopen_or_die("data_ext.txt", "at", "Problem with opening 'data_ext.txt'"); 
						fprintf(f_data_ext," %d ", new_code);
						for (i=1; i<=nt; i++) fprintf(f_data_ext," %d ", a[i]);
						//printf("nt=%d",nt);
						fprintf(f_data_ext,"\n");	
						fclose(f_data_ext);
						if (new_code<M1)	
						{
							printf(" g=%d for M= %d from M=%d \n", g,new_code,M0);  M1=new_code;
							f_report=fopen_or_die("in_out.txt", "at", "Problem with opening in_out.txt"); 
							fprintf(f_report,"\ng=%d for M= %d\n", g, M1);
							for (i=1; i<=n; i++) fprintf(f_report," %d ", a[i]);
							fclose(f_report);
						}
					}
				}
		
			} //if new code
		} // iterations for fixed input
	}  // M0<M
	if (is_keyboard_hit())
	{
		letter=get_character();
		if (letter=='x') 
		{
			f_data_ext=fopen_or_die("data_ext.txt", "at", "Problem with opening data_ext.txt"); 
			fprintf(f_data_ext," %d\n ", -1);
			fclose(f_data_ext);
			break;
		}
	}
}  // while 1
fclose(f_data);
break; // case 4
//*************************************************************//
case 5: // filter data.txt  to data_ext.txt
if (Hflag==0){ printf("\n Matrix H is not known yet\n"); break;}

f_data=fopen_or_die("data.txt", "rt", "Problem with opening 'data.txt");

f_data_ext=fopen_or_die("data_ext.txt", "wt", "Problem with opening 'data_ext.txt");
fclose(f_data_ext);

M0=1; I=0; J=0;
while (M0>0)
{
	I++;
	// read candidate
	fscanf_1(f_data,"%d", M0);  if (M0<0) break;
	for (i=1; i<=cumwcol[ninp]; i++){ fscanf_1(f_report,"%d", j); A[i]=j;}
	// select subsequence
	for (i=1; i<=nt; i++) a[i]=A[map[i]];
	//for (i=cumwcol[nf-1]+1; i<=cumwcol[nl]; i++) a[i-cumwcol[nf-1]]=A[i];

	new_code=check_eq(a, eqs, nodes);
	
	if (new_code>0)	
	{
		Mmax=M;
		new_code=min_module(a, eqs, nodes, &Mmax);
		if (new_code<=M)
		{
			if ((exact==0)||((FLAGS[exact]==0)&&(exact>0)))
			{
				J++;
				f_data_ext=fopen_or_die("data_ext.txt", "at", "Problem with opening 'data_ext.txt");
				fprintf(f_data_ext,"%d ", new_code);
				for (i=1; i<=cumwcol[ninp]; i++) fprintf(f_data_ext," %d ", A[i]);
				fprintf(f_data_ext,"\n");	
				fclose(f_data_ext);
			}
		}
	}
	if (is_keyboard_hit())
	{
		letter=get_character();
		if (letter=='x') break; 
	}
}//2
	
f_data_ext=fopen_or_die("data_ext.txt", "at", "Problem with opening 'data_ext.txt'");
fprintf(f_data_ext,"%d ", -1);
fclose(f_data_ext);

fclose(f_data);
printf("%d codes found from %d\n",J,I);
break; // case 5

case 6: //Scenario-based code generation
	f_report=fopen_or_die("scenario.txt", "rt", "The file 'scenario.txt' not found\n");

	// READ DATA
	fscanf_1(f_report,"%d", r);  // number of rows
	skip_until_newline(f_report);
    fscanf(f_report,"%d %d", &n, &ns);  // total number of columns and # of columns in subcode 
    skip_until_newline(f_report);
	fscanf(f_report,"%d %d", &g, &gs);  // target girth ansd subcode girth
	skip_until_newline(f_report);
	fscanf_1(f_report,"%d", M);  // tailbiting length
	skip_until_newline(f_report);
	fscanf_1(f_report,"%d", exact);  exact*=M;// exact tailbiting length
	skip_until_newline(f_report);
	fscanf_1(f_report,"%d", nb);  // starting length of base code 
	if (nb<ns) {printf("\n too short base code"); fclose(f_report); break;}
	if (nb>n) {printf("\n too long base code"); fclose(f_report); break;}
	skip_until_newline(f_report);
	fscanf_1(f_report,"%d", flag_continue);  // Minimum number of codes per iteration 
	skip_until_newline(f_report);
	fscanf_1(f_report,"%d", min_codes);  // Minimum number of codes per iteration 
	skip_until_newline(f_report);
	fscanf_1(f_report,"%d", max_codes);  // Maximum number of codes per iteration 
	skip_until_newline(f_report);
	fscanf_1(f_report,"%d", N);  // Max number of cycles per iteration 
	skip_until_newline(f_report);
	// scan full HM 
	for (i=0; i<r; i++) for (j=0; j<n; j++)	
	{
		fscanf_1(f_report,"%d", I); 
	    HM[i][j]=I;    // Full matriix with 0,1,2  
		T[i][j]=(I>0); // binary copy
	}
    fclose(f_report);

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
	

	// mapping for nonzero positions
	I=0;
	for (i=nb-ns+1; i<=nb; i++)
	{
		for (j=cumwcol[i-1]+1; j<=cumwcol[i]; j++)
		{
			I++; map[I]=j;
		}
	};

if (!flag_continue)
{	
	// First set of equations

    // obtain base matrix of proper size for constructing equations
	I=0;
	for (i=1; i<=nb; i++)
	{
		for (j=0; j<r; j++)	H[j][I]=T[j][i-1]; 
		I++;
	};
	

	printf("\n n=%d, r=%d, g=%d, M=%d \n",nb, r, g, M); 
 
	//   Tanner graph 
	rt=nb+r; 	nt=tanner(nb,r);  	                             //input H, output HH
	for (i=0; i<rt; i++) for (j=0; j<nt; j++) H[i][j]=HH[i][j];  // Tanner graph PCM 
	printf("\n Tanner graph size: nt=%d, rt=%d \n", nt, rt);
   // Generating first set of equations
	printf("Generating first equations  \n");
	ok=0;
	ok=eq_generation(nt,rt, g, &eqs, &nodes); 
	if (ok)  printf("ok: g=%d, eqs=%d, nodes=%d\n",g, eqs, nodes);
	else  {printf("Problem!"); goto finish;}

	// SECOND set of equations

    // obtain base matrix of proper size for constructing equations
	I=0;
	for (i=nb-ns+1; i<=nb; i++)
	{
		for (j=0; j<r; j++)	H[j][I]=HM[j][i-1]>0; 
		I++;
	};
	
	printf("\n n=%d, r=%d, g=%d, M=%d \n",nb, r, gs, M); 
 	//   Tanner graph 
	rts=ns+r; 	nts=tanner(ns,r);  	                             //input H, output HH
	for (i=0; i<rt; i++) for (j=0; j<nts; j++) H[i][j]=HH[i][j];  // Tanner graph PCM 
	printf("\n Tanner graph size: nt=%d, rt=%d \n", nts, rts);
   // Generating first set of equations
	printf("Generating second equations  \n");
	ok=0;
	ok=eq_generations(nts,rts, gs, &eqss, &nodess); 
	if (ok)  printf("ok: g=%d, eqs=%d, nodes=%d\n",gs, eqss, nodess);
	else  {printf("Problem!"); goto finish;}

    for (i=1; i<=nt; i++) as[i]=a[i]=0; M1=M;

	// Random generating first nb columns

	name[0]='d'; name[1]='a'; name[2]='t'; name[3]='a';
	letter=48+(nb)/10;
	name[4]=letter;
	name[5]=48+(nb)%10;
    name[6]='.';name[7]='t';name[8]='x';name[9]='t'; name[10]=NULL;	
	//for (i=0; i<10; i++)
	//printf(" %c ",name[i]);


	f_data=fopen_or_die(name, mwt, "Problem with opening 'data.txt"); fclose(f_data);
	if (1) 
	{
		printf("Creating random %d columns \n", nb);
	    I=0; J=0;
		while ((I<N)&& (J<max_codes))
		{
			I++;  if(I%100000==0) printf("%d tested \n",I);
			for (i=1; i<=nt; i++) 
			{                 
				if (am[i]==1 )  a[i]=next_random_int(0,M); 	else a[i]=0;
			}
			for (i=1; i<=nts; i++) as[i]=a[map[i]];
		    new_code=check_eq(a, eqs, nodes);
			if ((gs>g)&&(new_code>0)) new_code=check_eqs(as, eqss, nodess);
    	    if (new_code>0)	
	        {
		        Mmax=M;
		        new_code=min_module(a, eqs, nodes, &Mmax);
		        if ((gs>g)&&(new_code<=M)) new_code=min_modules(as, eqss, nodess, &Mmax);
		        if (new_code<=M)
		        {
			        if ((exact==0)||((FLAGS[exact]==0)&&(FLAGSS[exact]==0)&&(exact>0)))
			        {
				        f_data=fopen_or_die(name, mat, "Problem with opening 'data.txt");
				        fprintf(f_data,"%d ", new_code);
				        for (i=1; i<=nt; i++) fprintf(f_data," %d ", a[i]);
				        fprintf(f_data,"\n");	
				        fclose(f_data); 
				        J++;  if (J%1000==0) printf ("%d codes written to 'data.txt'\n",J);
				        if (new_code<M1)	
				        {		
					        f_data=fopen_or_die("in_out.txt", "at", "Problem with opening 'in_out.txt'");
					        fprintf(f_report,"\ng=%d for M= %d\n", g, new_code);
					        for (i=1; i<=nt; i++) fprintf(f_report," %d ", a[i]);
					        fclose(f_report);
					        M1=new_code;
					        printf(" g=%d for M= %d\n", g, new_code);
				        }
			        }
		        }
	        }
		
			if (is_keyboard_hit())
			{
				letter=get_character();
				if (letter=='x') break;
				if (letter=='i') printf(" nb= %d, %d tested %d found \n", nb, I, J);
			}
     	}
	}
//}//2

printf("\n  %d codes selected, %d tested\n", J, N);
f_data=fopen_or_die(name, mat, "Problem!");
fprintf(f_data," %d ", -1);
fclose(f_data);
	
// LOOP OVER SUBCODES
//nb++;  
} // if start with random code

	
while (nb<n)
{
	name[1]='d'; 	name[1]='a';name[1]='t';name[1]='a';
	name[4]=48+(nb)/10;
	name[5]=48+(nb)%10;
    name[6]='.';name[7]='t';name[8]='x';name[9]='t'; name[10]=NULL;
	nb++;
	nextname[0]='d';nextname[1]='a';nextname[2]='t';nextname[3]='a';	
	nextname[4]=48+(nb)/10;
	nextname[5]=48+(nb)%10;
	nextname[6]='.';nextname[7]='t';
	nextname[8]='x';nextname[9]='t'; nextname[10]=NULL;
	
    printf("\nExtending %d colums to  %d columns\n", nb-1, nb);
//update equations
   // obtain base matrix of proper size for constructing equations
	I=0;
	for (i=1; i<=nb; i++)
	{
		for (j=0; j<r; j++)	H[j][I]=HM[j][i-1]>0; 
		I++;
	};
	
	printf("\n n=%d, r=%d, g=%d, M=%d \n",nb, r, g, M); 
 
	//   Tanner graph 
	rt=nb+r; 	nt=tanner(nb,r);  	                             //input H, output HH
	for (i=0; i<rt; i++) for (j=0; j<nt; j++) H[i][j]=HH[i][j];  // Tanner graph PCM 
	printf("\n Tanner graph size: nt=%d, rt=%d \n", nt, rt);
   // Generating first set of equations
	printf("Generating first equations  \n");
	ok=0;
	ok=eq_generation(nt,rt, g, &eqs, &nodes); 
	if (ok)  printf("ok: g=%d, eqs=%d, nodes=%d\n",g, eqs, nodes);
	else  {printf("Problem!"); goto finish;}

	// SECOND set of equations
	// obtain base matrix of proper size for constructing equations
	I=0;
	for (i=nb-ns+1; i<=nb; i++)
	{
		for (j=0; j<r; j++)	H[j][I]=HM[j][i-1]>0; 
		I++;
	};
	
	printf("\n n=%d, r=%d, g=%d, M=%d \n",nb, r, gs, M); 
 	//   Tanner graph 
	rts=ns+r; 	nts=tanner(ns,r);  	                             //input H, output HH
	for (i=0; i<rt; i++) for (j=0; j<nts; j++) H[i][j]=HH[i][j];  // Tanner graph PCM 
	printf("\n Tanner graph size: nt=%d, rt=%d \n", nts, rts);
   // Generating first set of equations
	printf("Generating second equations  \n");
	ok=0;
	ok=eq_generations(nts,rts, gs, &eqss, &nodess); 
	if (ok)  printf("ok: g=%d, eqs=%d, nodes=%d\n",gs, eqss, nodess);
	else  {printf("Problem!"); goto finish;}

    for (i=1; i<=nt; i++) as[i]=a[i]=0; 

	// Generate mask for last branch(es)
	I=0;
	for (i=nb-ns+1; i<=nb; i++)
	{
		for (j=cumwcol[i-1]+1; j<=cumwcol[i]; j++)
		{
			I++; map[I]=j;
		}
	};

	f_data=fopen_or_die(name, "rt", "Problem with opening 'data.txt");

	f_data_ext=fopen_or_die(nextname, mwt, "Problem with opening 'data_ext.txt"); 	fclose(f_data_ext);
	//fclose(f_data_ext);

	I=0; J=0; 
	M1=M;
	printf("\n Press 'x' to stop ");
	I=0; J=0;
	while ((I<N)&& (J<max_codes))
	{
		I++;
		// read candidate
		fscanf_1(f_data,"%d", M0); 
		if (M0<0)  //refresh
		{
			fclose(f_data);
			f_data=fopen_or_die(name, mrt, "Problem with opening 'data.txt");
			fscanf_1(f_data,"%d", M0); 
		}

		if (M0<=M)
		{
			for (j=1; j<=cumwcol[nb-1]; j++) 
			{
				fscanf_1(f_data,"%d", ans); a[j]=ans;
			}
			// extend it
			for (j=0; j<mini(M0,5); j++)
			{
				for (i=cumwcol[nb-1]+1; i<=nt; i++) 
					if (am[i]==1) a[i]=next_random_int(0,M);  else a[i]=0; 
				for (i=1; i<=nts; i++) as[i]=a[map[i]];

				new_code=check_eq(a, eqs, nodes);
				if ((gs>g)&&(new_code>0)) new_code=check_eqs(as, eqss, nodess);

				if (new_code>0)
			    {
				new_code=min_module(a, eqs, nodes, &Mmax);
				if ((gs>g)&&(new_code<=M)) new_code=min_modules(as, eqss, nodess, &Mmax);
				if (new_code<=M)	
				{
				    if ((exact==0)||((FLAGS[exact]==0)&&(exact>0)))
					{
						J++;
	    				f_data_ext=fopen_or_die(nextname, mat, "Problem!");
						fprintf(f_data_ext," %d ", new_code);
						for (i=1; i<=nt; i++) fprintf(f_data_ext," %d ", a[i]);
						//printf("nt=%d",nt);
						fprintf(f_data_ext,"\n");	
						fclose(f_data_ext);
						if (new_code<M1)	
						{
							printf(" g=%d for M= %d from M=%d \n", g,new_code,M0);  M1=new_code;
							f_report=fopen_or_die("in_out.txt", mat, "Problem with opening 'in_out.txt'");
							fprintf(f_report,"\ng=%d for M= %d\n", g, M1);
							for (i=1; i<=n; i++) fprintf(f_report," %d ", a[i]);
							fclose(f_report);
						}
					}
				}
		
			} //if new code
		} // iterations for fixed input
	}  // M0<M
	if (is_keyboard_hit())
	{
		letter=get_character();
		if (letter=='x') goto finish;
		if (letter=='i') printf(" nb= %d, %d tested %d found \n", nb, I, J);
    }
	}  // while 1
	f_data_ext=fopen_or_die(nextname, mat, "Problem!"); 
	fprintf(f_data_ext," %d\n ", -1);
	fclose(f_data_ext);
	printf(" %d codes are found from %d \n", J, I);  
	if (J<min_codes) {printf("Too few codes of length %d", nb);goto finish;}

	//nb++;

}
fclose(f_data);
break; // case 6

//*************************************************************//

case 7: // Average girth/spectrum	   
f_report=fopen_or_die("a.cfg", "rt", " File 'a.cfg' not found");


printf("Input # of trials "); 
scanf("%d", &N);
fscanf_1(f_report,"%d", M);  // tailbiting length
skip_until_newline(f_report);
fscanf(f_report,"%d %d",&r, &n);  // number of rows and columns 
skip_until_newline(f_report);
	// scan full HM 
for (i=0; i<r; i++) for (j=0; j<n; j++)	
{
	fscanf_1(f_report,"%d", I); 
    HM[i][j]=I;    // Labeled matriix with -1, 0..M  
}
fclose(f_report);


for (J=0; J<N; J++)
{
	//generate random subcode
	nb=0;
	while (nb<=2) 
	{ 
		nb=0;
		for (I=1; I<=n; I++)
		if (rand()&1) {nb++; P[nb]=I;}
	}

	// obtain base matrix of proper size 
	I=0;
	for (i=1; i<=nb; i++)
	{
		for (j=0; j<r; j++)	
		{
			H[j][I]=HM[j][P[i]-1];
			T[j][I]=H[j][I]>=0;  // Base matrix
		}
		I++;
	};

	//   Tanner graph in poly form
	rt=nb+r; 	nt=tanner(nb,r);  	                             //input H, output HH
	for (i=0; i<rt; i++) for (j=0; j<nt; j++) HD[i][j]=HH[i][j];  // Tanner graph PCM 
	printf("\n Tanner graph size: nt=%d, rt=%d \n", nt, rt);
	 
	// transform to binary form 
	hpoly2tb(nt, rt, M);   // HD --> 
	t= graph_girth(nt, rt);
	
	

}




for (i=1; i<=cumwcol[ninp]; i++){ fscanf_1(f_report,"%d", I); A[i]=I;}
// select subsequence for analysis
for (i=1; i<=nt; i++) a[i]=A[map[i]];

// for (i=cumwcol[nf-1]+1; i<=cumwcol[nl]; i++) a[i-cumwcol[nf-1]]=A[i];

//I=1; 
gtmp=g; Mmax=100000;
ans=1;
while (ans!=0)
{
//	new_code=check_eq(a, eqs, nodes);
	ok=eq_generation(nt,rt, gtmp, &eqs, &nodes);
	if (ok)
	{
		new_code=check_eq(a, eqs, nodes);
		if (new_code<0) 
		{
			printf(" gfree< %d, wrong eq_number=%d \n", gtmp, -new_code);
			if ((ans==1)||(ans==-2)) ans=-2;  // move down
			else ans=0;
		}
		else 
		{
			new_code=min_module(a, eqs, nodes, &Mmax);
			printf(" g>=%d for  M= %d and all M>Mmax=%d\n", gtmp, new_code, Mmax); 
			fprintf(f_report,"\ng=%d for M= %d and all M>Mmax=%d\n", gtmp, new_code, Mmax);
			I=0;
			if (new_code>0){
			for (i=new_code; i<=Mmax; i++)
			{
				if (FLAGS[i]==0) 
				{
					fprintf(f_report,"%d ", i); printf("%d ", i);
					I++ ; if (I>20) {fprintf(f_report,"\n"); printf("\n"); I=0;}	 
				}
			}
			fprintf(f_report,"..."); printf("...\n");}
			if ((ans==1)||(ans==2)) ans=2; //move app
			else ans=0;
			//break;
		}
		gtmp+=ans; 
	}
	else break;
}
fclose(f_report);
Hflag=0;
break; //case 7 

case 0: goto finish;  

//*************************************************************//
// case 
} //end while (1)

}  // end main
finish: ;
}



 
	
