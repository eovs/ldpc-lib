#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "graph.h"

#include "commons_portable.h"

int main(int argc, char *argv[])
{
    FILE *f_report, *flog;

	int n,       //  max length of base code
		g, gs,//  girth
	    r;      //  base code check symbols number

	int i,  j, J, M;
	int ok, good_code;

	// Variables for EQ constructing
	char letter;

	// SIMULATION 
	int maxiter, num_snrs;

	float SNRS[10], SNR,  PR_ERR[10][2], BEST_FERS[10], best_fer=0.1, pr_err[2];
	int    Nerr, Nexp;

	// Constructing Base matrix
	int DD[CMAX], Jstart, N_codes,
		trials[2], score[4], s, degs[13],
		confs, code,
		msnd;                   // max symbol node degree


    // prepare
    initial_random_seed = 1;

    gen_divtable();

    //Scenario-based code generation
    f_report = fopen_or_die(argv[1], "rt", "The scenario file not found");

	// READ DATA
    fscanf(f_report,"%d", &r);  // number of rows
    skip_until_newline(f_report);
    fscanf(f_report,"%d ", &n );  // total number of columns
    skip_until_newline(f_report);
    fscanf(f_report,"%d ", &g );  // target girth
    skip_until_newline(f_report);
    fscanf(f_report,"%d", &M);  // tailbiting length
    skip_until_newline(f_report);
    fscanf(f_report,"%d", &num_snrs);  // number of SNRS
    skip_until_newline(f_report);
    for (i=0; i<num_snrs; i++)
    {
        fscanf(f_report,"%f", &SNR);
        SNRS[i]=SNR;
        BEST_FERS[i]=best_fer;
    }
    skip_until_newline(f_report);
    fscanf(f_report,"%d %d", &i, &j);//number of trials for base matrix and its column  
    trials[0]=i; trials[1]=j;
    skip_until_newline(f_report);
    fscanf(f_report,"%d", &N_codes);  //  number of tested random codes 
    skip_until_newline(f_report);
    fscanf(f_report,"%d", &maxiter);  //max number of BP iterations 
    skip_until_newline(f_report);
    fscanf(f_report,"%d %d", &Nexp, &Nerr);//number of exper and number of frame errors in experiments 
    skip_until_newline(f_report);
    fscanf(f_report,"%d", &Jstart);  //  start configuration number 
    skip_until_newline(f_report);
    fscanf(f_report,"%d", &j);  DD[2]=j; //  # of columns of weight 2 
    skip_until_newline(f_report);
    fscanf(f_report,"%d", &s);  //number of different degrees
    skip_until_newline(f_report);
    for (i=0; i<s; i++)
    {
	    fscanf(f_report,"%d", &j);
	    degs[i]=j;
    }
    fclose(f_report);

    flog = fopen_or_die(argv[2], "wt", "Log file cannot be opened");
    fprintf(flog,"%d %d  // size of base matrix \n", r,n);
    fprintf(flog,"%d     // TB length \n", M);
    fprintf(flog,"%f     // SNR \n", SNR);
    fprintf(flog,"%d    // maxiter \n", maxiter);
    fclose(flog);

    confs=list_conf_n(n-DD[2], s);

    for (J=confs-1-Jstart; J>=0; J--)
    {
	    DD[0]=DD[1]=0;
	    for (i=3; i<=n; i++) DD[i]=0;
        // Transform configuration to column weight distribution
	    j=0;
	    while (LCONF[J][j]>0)
	    {
		    msnd=degs[j];
		    DD[msnd]=LCONF[J][j];
		    j++;
	    }
	    for (i=0; i<=msnd; i++) printf("%d ", DD[i]); printf("\n");
	    // Create base matrix
	    ok=rand_base_matrix_g(r, n, msnd, DD, trials, score, 1);
	    if (ok<0) return 0;
	    for (code=0; code<=10; code++)
	    {
		    if (score[3]>300) g=6;
		    ok=0; gs=g+2;
		    while (ok<=0)
		    {
			    for (i=0; i<r; i++) for (j=0; j<n; j++) T[i][j]=HM[i][j]>0; // binary copy
			    gs-=2;
			    if (gs<=4) break;
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
	    	    // was a senseless PR_ERR[s][1] = BEST_FERS[s];
        		pr_err[1]=BEST_FERS[s];
		        i=bp_simulation(r, n, M, maxiter, Nerr, Nexp, SNRS[s], pr_err) ;
		        PR_ERR[s][0]=pr_err[0]; PR_ERR[s][1]=pr_err[1];
		        if (i<0) PR_ERR[s][1]=best_fer*10; // wrong matrix
		        if (i==0)                          // wrong coding
		        {
			        flog = fopen_or_die("H.txt", "wt", "H.txt cannot be opened");
			        for (j=0; j<r; j++) {for (i=0; i<n; i++)  fprintf(flog,"%3d", HD[j][i]); 
			        fprintf(flog,"\n "); }	
			        fclose(flog);
			        J=r*M;
			        flog = fopen_or_die("cw.txt", "wt", "cw.txt cannot be opened");
			        for (j=r; j<n; j++) {
			            for (i=0; i<M; i++) {
			                fprintf(flog,"%2d", codeword[J]);
			                J++;
			            }
    			        fprintf(flog,"\n ");
			        }
			        fclose(flog);
			        return 0;
		        }
		        printf("SBEST_FER=%5.3e, code=%d, pr_err[0] = %f, pr_err[1] = %f\n",
		                BEST_FERS[s], code, PR_ERR[s][0], PR_ERR[s][1]);
		        for (i=0; i<=msnd; i++) printf("%d ", DD[i]); printf("\n");

    	    	// PRINTING
	    	    if ((PR_ERR[s][1] > BEST_FERS[s]*2.0) || ((PR_ERR[s][1] > BEST_FERS[s]) && (code>2)))
			        break;
        		if (PR_ERR[s][1] < BEST_FERS[s]*1.25)
		        	good_code=1;
        		if (PR_ERR[s][1] < BEST_FERS[s]) BEST_FERS[s]=PR_ERR[s][1];
		    }   // for s
		    if (good_code==1)
		    {
			    flog = fopen_or_die(argv[2], "at", "log.log cannot be opened");
			    for (i=0; i<r; i++)
			    {             // row number
				    for (j=0; j<n; j++)   fprintf(flog,"%3d ", HD[i][j]); fprintf(flog,"\n");
			    }
			    fprintf(flog," g = %d, J=%d\n", gs, J);
			    for (i=0; i<=msnd; i++) fprintf(flog,"%d ", DD[i]); fprintf(flog,"\n");
			    for (s=0; s<num_snrs; s++) {
			        fprintf(flog,"SNR=%5.3f, BER=%5.3e,FER=%5.3e\n",SNRS[s], PR_ERR[s][0],PR_ERR[s][1]);
			    }
			    fclose(flog);
			    printf("Good code written!\n");
    		} //
        }
    } // loop over configurations
}  // end main
