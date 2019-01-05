#include <time.h>
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>
#include "graph.h"

//
main(int argc, char *argv[])
{
	FILE *f_report,// *fe, 
		 *flog; 	
	//FILE *f_data_ext;
	errno_t err;

	int n,       //  max length of base code
		g, gs,//  girth
	    r, r0;      //  base code check symbols numbers
		                 
	
	int i,  j, J,J0, M;//  Mmax; 
	int ok, good_code;//, ans, t; 
//	int flag;//, M1;// , flag_continue;
	int num_restr, rw[20], rmin[20], rmax[20];  
	 
	
	// Variables for EQ constructing
	char letter;

	// SIMULATION 
	int maxiter, num_snrs;

	float SNRS[10], SNR,  PR_ERR[10][2], BEST_FERS[10], best_fer=0.1,
		pr_err[2];
	int    Nerr, Nexp;

	// Constructing Base matrix
		int DD[CMAX],DD0[CMAX],N_codes,
		trials[2], score[4], s, degs[13], d2,d02,
		s0, degs0[13],
		confs0,confs, code,
		msnd0,msnd;                   // max symbol node degree
  int  DDT[RMAX];
/*
err=fopen_s(&flog, "h.txt", "rt");
fscanf_s(flog,"%d %d ", &r,&n);
for (i=0; i<r;i++) for (j=0; j<n; j++) {fscanf_s(flog,"%d", &J); H[i][j]=J;}
fclose(flog);
rt=n+r; 	nt=tanner(n,r);
for (i=0; i<rt; i++) for (j=0; j<nt; j++) H[i][j]=HH[i][j];
g=count_loops(nt, rt, score);
*/
		 

// prepare 
srand((unsigned)time(NULL)); //    
//srand(1);
//Hflag=0; 	
gen_divtable();

//Scenario-based code generation
err=fopen_s(&f_report,argv[1], "rt");
if(err!=0 ) {printf( "The scenario file not found\n" ); };

	// READ DATA
fscanf_s(f_report,"%d %d", &r, &r0);  // number of rows
fscanf_s(f_report,"%c", &letter); while (letter!=0x0a) fscanf_s(f_report,"%c" , &letter,1); 
fscanf_s(f_report,"%d ", &n );  // total number of columns  
fscanf_s(f_report,"%c", &letter); while (letter!=0x0a) fscanf_s(f_report,"%c" , &letter,1); 	
fscanf_s(f_report,"%d ", &g );  // target girth 
fscanf_s(f_report,"%c", &letter); while (letter!=0x0a) fscanf_s(f_report,"%c" , &letter,1);
fscanf_s(f_report,"%d", &M);  // tailbiting length
fscanf_s(f_report,"%c", &letter); while (letter!=0x0a) fscanf_s(f_report,"%c" , &letter,1);
fscanf_s(f_report,"%d", &num_snrs);  // number of SNRS
fscanf_s(f_report,"%c", &letter); while (letter!=0x0a) fscanf_s(f_report,"%c" , &letter,1);
for (i=0; i<num_snrs; i++)
{
fscanf_s(f_report,"%f", &SNR);
SNRS[i]=SNR;
BEST_FERS[i]=best_fer;
}
fscanf_s(f_report,"%f", &SNR);  // SNR
fscanf_s(f_report,"%c", &letter); while (letter!=0x0a) fscanf_s(f_report,"%c" , &letter,1); 
fscanf_s(f_report,"%d %d", &i, &j);//number of trials for base matrix and its column  
trials[0]=i; trials[1]=j;
fscanf_s(f_report,"%c", &letter); while (letter!=0x0a) fscanf_s(f_report,"%c" , &letter,1);
fscanf_s(f_report,"%d", &N_codes);  //  number of tested random codes 
fscanf_s(f_report,"%c", &letter); while (letter!=0x0a) fscanf_s(f_report,"%c" , &letter,1); 
fscanf_s(f_report,"%d", &maxiter);  //max number of BP iterations 
fscanf_s(f_report,"%c", &letter); while (letter!=0x0a) fscanf_s(f_report,"%c" , &letter,1); 
fscanf_s(f_report,"%d %d", &Nexp, &Nerr);//number of exper and number of frame errors in experiments 
fscanf_s(f_report,"%c", &letter); while (letter!=0x0a) fscanf_s(f_report,"%c" , &letter,1); 
fscanf_s(f_report,"%d", &s0);  //number of different degrees
fscanf_s(f_report,"%c", &letter); while (letter!=0x0a) fscanf_s(f_report,"%c" , &letter,1); 
for (i=0; i<s0; i++)
{
	fscanf_s(f_report,"%d", &j);
	degs0[i]=j;
}
fscanf_s(f_report,"%c", &letter); while (letter!=0x0a) fscanf_s(f_report,"%c" , &letter,1); 
fscanf_s(f_report,"%d", &s);  //number of different degrees
fscanf_s(f_report,"%c", &letter); while (letter!=0x0a) fscanf_s(f_report,"%c" , &letter,1); 
for (i=0; i<s; i++)
{
	fscanf_s(f_report,"%d", &j);
	degs[i]=j;
}
fscanf_s(f_report,"%c", &letter); while (letter!=0x0a) fscanf_s(f_report,"%c" , &letter,1); 
fscanf_s(f_report,"%d", &num_restr);  //number of restrictions

for (i=0; i<num_restr; i++)
{
	fscanf_s(f_report,"%c", &letter); while (letter!=0x0a) fscanf_s(f_report,"%c" , &letter,1); 
	fscanf_s(f_report,"%d", &j); rw[i]=j; 
    fscanf_s(f_report,"%d", &j); rmin[i]=j; 	
    fscanf_s(f_report,"%d", &j); rmax[i]=j; 	 
}



fclose(f_report);
err=fopen_s(&flog, argv[2], "wt");
fprintf(flog,"%d %d %d // size of base matrix \n", r,r0, n);
fprintf(flog,"%d     // TB length \n", M);
fprintf(flog,"%f     // SNR \n", SNR);
fprintf(flog,"%d    // maxiter \n", maxiter);
fclose(flog);


//confs=list_conf(n-DD[2], s);
d02=r0-1;                         // # of weight 2 columns in 1st mattrix 
d2=r-r0-1;                      // # of weight 2 columns in 2nd mattrix
confs0=list_conf_n(n-r0, s0);
// Copy LCONF[4000][12];
for (i=0; i<confs0; i++) for (j=0; j<s0; j++) {LCONF0[i][j]=LCONF[i][j]; LCONF[i][j]=0;}  
confs=list_conf_n(n-r, s);

for (J0=confs0-1; J0>0; J0--)
{
	// Transform configurations to column weight distributions
	for (i=0; i<=n; i++) {DD0[i]=0;}
	j=0;
	while (LCONF0[J0][j]>0)
	{
		msnd0=degs0[j];
		DD0[msnd0]=LCONF0[J0][j];
		j++;
	}
	

for (J=confs-1; J>0; J--)
{
	printf("Configurations %d, %d \n", J0, J);
	for (i=0; i<=n; i++) { DD[i]=0;}
// Transform configurations to column weight distributions
	j=0;
	while (LCONF[J][j]>0)
	{
		msnd=degs[j];
		DD[msnd]=LCONF[J][j];
		j++;
	}

	for (i=0; i<=msnd0; i++) printf("%d ", DD0[i]); printf("\n");
	for (i=0; i<=msnd; i++) printf("%d ", DD[i]); printf("\n");
	// Create base matrix
	ok=rand_base_matrix_g(num_restr,rw,rmin,rmax, r0, r, n, msnd0, msnd, DD0, DD, DDT, trials, score, 1);
	if (ok==-1) return;
	if (ok==-2) goto nextconf;
	for (code=0; code<=10; code++) 
	{
		if (score[3]>300) g=6;
		ok=0; gs=g+2;
		while (ok<=0)
		{
			for (i=0; i<r; i++) for (j=0; j<n; j++) T[i][j]=HM[i][j]>0; // binary copy
			gs-=2;
			if (gs<=4) goto nextconf;
			ok=generate_code(r, n, gs, M, N_codes);
			if (ok<=0) 
			{ 
				for (i=0; i<r; i++) for (j=0; j<n; j++) T[i][j]=HM[i][j]>0;
				ok=generate_code(r, n, gs, MMAX-1, N_codes);
				for (i=0; i<r; i++) for (j=0; j<n; j++) HD[i][j]%=M; 
			}
		}
		if (gs==4) break;
		
		good_code=0;

		for (s=0; s<num_snrs; s++)
		{   // LOOP OVER SNRs
		
        
		PR_ERR[s][1]=BEST_FERS[s];
		i=bp_simulation(r, r0, n, M, maxiter, Nerr, Nexp, SNRS[s], pr_err) ;
		PR_ERR[s][0]=pr_err[0]; PR_ERR[s][1]=pr_err[1];
		if (i<0) PR_ERR[s][1]=best_fer*10; // wrong matrix
		if (i==0)                          // wrong coding 
		{
			fopen_s(&flog, "H.txt", "wt");
			for (j=0; j<r; j++) {for (i=0; i<n; i++)  fprintf(flog,"%3d", HD[j][i]); 
			fprintf(flog,"\n "); }	
			fclose(flog);
			J=r*M;
			fopen_s(&flog, "cw.txt", "wt");
			for (j=r; j<n; j++) 
			{for (i=0; i<M; i++) { fprintf(flog,"%2d", codeword[J]); J++;} 
			fprintf(flog,"\n "); }	
			fclose(flog);
			return;
		}

		printf("SBEST_FER=%5.3e , code=%d  J0=%d(%d), J0=%d(%d) \n",BEST_FERS[s], code, J0, confs0, J, confs); 
		for (i=0; i<=msnd0; i++) printf("%d ", DD0[i]); printf("\n");
		for (i=0; i<=msnd; i++) printf("%d ", DD[i]); printf("\n");
	    for (i=0; i<=msnd+msnd0; i++) printf("%d ", DDT[i]); fprintf(flog,"\n");

		// PRINTING
		if ((PR_ERR[s][1] > BEST_FERS[s]*2.0) || ((PR_ERR[s][1] > BEST_FERS[s]) && (code>2))) 
			break;
		if (PR_ERR[s][1] < BEST_FERS[s]*1.25) 
			good_code=1;
		if (PR_ERR[s][1] < BEST_FERS[s]) BEST_FERS[s]=PR_ERR[s][1];
		}   // for s
		if (good_code==1)
		{
			err=fopen_s(&flog, argv[2], "at");
			if (err!=0) printf("Problem with opening 'log.log' ");
			for (i=0; i<r; i++) 
			{             // row number
				for (j=0; j<n; j++)   fprintf(flog,"%3d ", HD[i][j]); fprintf(flog,"\n");
			}
			fprintf(flog," g = %d, J=%d\n", gs, J);

			for (i=0; i<=msnd0; i++) fprintf(flog,"%d ", DD0[i]); fprintf(flog,"\n");
	        for (i=0; i<=msnd; i++) fprintf(flog,"%d ", DD[i]); fprintf(flog,"\n");
	        for (i=0; i<=msnd+msnd0; i++) fprintf(flog,"%d ", DDT[i]); fprintf(flog,"\n");
 
			for (s=0; s<num_snrs; s++)
			fprintf(flog,"SNR=%5.3f, BER=%5.3e,FER=%5.3e\n",SNRS[s], PR_ERR[s][0],PR_ERR[s][1]);   
			fclose(flog);
			//if (PR_ERR[s][1] < BEST_FERS[s]) BEST_FERS[s]=PR_ERR[s][1];

		} // 
		//printf("b,c, SNR= %d , %d, %f \n", b,c, SNR);

	} // over 10 codes
	nextconf: 		;
}


} // loop over configurations
}  // end main





 
	