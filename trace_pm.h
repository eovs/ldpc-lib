#ifndef _TRACE_PM_H_

#define _TRACE_PM_H_

#include "3d_array.h"

#define NATURAL_POLYNOM

typedef struct  
{
	ARRAY array1;
	ARRAY array2;
}TANNER_MON_RES;

typedef struct  
{
	ARRAY array1;
	ARRAY array2;
}HP2A_MON_RES;
#ifdef NATURAL_POLYNOM
typedef struct  
{
	ARRAY ret1;
	ARRAY ret2;
	ARRAY ret3;
	ARRAY ret4;
	ARRAY pol;
	ARRAY ase;
}MUL_MAT_MAT_MON_RES;
#else
typedef struct  
{
	ARRAY ret1;
	ARRAY ret2;
	ARRAY ret3;
	ARRAY ret4;
	int	  ret5;
}MUL_MAT_MAT_MON_RES;
#endif

typedef struct 
{
	int degree;
	int *pos;
	int *coef;
	int *ace;
} POLYNOM; 

void trace_bound_pol_mon_pm( ARRAY HD, int M, int gmax, int gtarget, int *SN, int *SA );
TANNER_MON_RES tanner_mon( ARRAY HD );
HP2A_MON_RES hp2a_mon_pm( ARRAY HD, ARRAY W, int M );

#ifdef NATURAL_POLYNOM
MUL_MAT_MAT_MON_RES mul_mat_mat_mon( int bb, ARRAY polB, ARRAY aceB, ARRAY A, ARRAY AC, int M);
#else
MUL_MAT_MAT_MON_RES mul_mat_mat_mon( ARRAY WB, ARRAY DB, ARRAY CB, ARRAY BS, ARRAY A, ARRAY AC, int M);
#endif

void unpack_pol( int wx, int dx[], int cx[], int xs[], int M, int pol[], int ase[] );
int pack_pol( int pol[], int ase[], int M, int dx[], int cx[], int xs[] );

#endif //_TRACE_PM_H_