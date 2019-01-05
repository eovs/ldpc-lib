#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#ifndef SKIP_MEX
#include <mex.h>
#endif

#include "decoders.h"

char const * const DEC_FULL_NAME[] = 
{
    "Belief Propagation",
    "Sum-Product",
    "Advanced Sum-Product",
    "Min-Sum",
    "Integer Min-Sum",
    "Integer Advanced Sum-Product",
    "FHT Sum-Product",
    "TDMP Advanced Sum-Product"
};

typedef struct  
{
	double val;
	int pos;
}ELEMENT;

//#define MAP_GRAPH_USE_INT
//#define MAP_MULT 10//12


#define MAKE_LIST
#define NMAIN 16

#define   COLUMN_BY_COLUMN
//#define NORM_MAX
#define NORM_SHIFT
#define SCALABLE


#define INP_FPP  16
#define ONE_INP (1 << INP_FPP)

#define PROB_FPP 16//16 //28  // less than 30 !!!
#define ONE_PROB (1 << PROB_FPP)

#define MAP_FPP  16//14//24
#define ONE_MAP (1 << MAP_FPP)

#define HAD_FPP 16//14
#define ONE_HAD (1 << HAD_FPP)

#define TMP_FPP 20
#define ONE_TMP (1 << TMP_FPP)

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define QMAX  64//256
#define RWMAX 64

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#define div_power2( x, n )    ((x) >> (n))
#define div_power2r( x, n )   (((x) + (1 << (n-1))) >> (n))

#define SOFT_FPP 12
#define P_FPP    32//32


#define MAX_UI16 0x7fff
#define ONE_SOFT (1 << SOFT_FPP)


#define MAX_SOFT (ONE_SOFT - 1) 
#define IASP_DEC_MAX_VAL (ONE_SOFT - 1)


#define INPUT_LIMIT 20.0

#define SP_DEC_MIN_VAL 0.000001
#define SP_DEC_MAX_VAL (1.0 - SP_DEC_MIN_VAL)


#define maxi( a, b )  (a) < (b) ? (b) : (a)
#define mini( a, b )  (b) < (a) ? (b) : (a)


static double mind(double a, double b) {if (a<b) return a; else return b;}
static double maxd(double a, double b) {if (a<b) return b; else return a;}
static double absd (double x)
{
	if (x<0) return -x; else return x;
}

UFLT16 fixed2uflt16( int xfpp , unsigned short x );
ui32 uflt2fixed( UFLT16 x, int fpp );

int pol_bank( int c[], int *pol, int *M );
int gf2( int m, int p, short *gf_log, short *gf_alog );
static void p2table( short **t, short **ti, short *gf_log, short *gf_alog, int list[], int q_bits, int n_used );

void hadamar( double x[], double y[], int q, int c );
void normalize( double *x, int n, int step, int fpp );
void normalize_abs( double *x, int n, int step, int fpp );
ui32 get_bitsize( ui32 d );
void sort( ELEMENT x[], int n, int k );


#ifndef SKIP_MEX
void unpackMatrix(double m[], int height, int width, double *result[] )
{
    int i,j;
	for( i = 0; i < height; ++i)
	{
		for( j = 0; j < width; ++j)
		{
			result[i][j] = m[j * height + i];
		}
	}
}

void unpackRow(double m[], int height, int width, double result[] )
{
    int j;
	for( j = 0; j < width; ++j)
	{
		result[j] = m[j];
	}
}


void unpackMatrix_double2int(double m[], int height, int width, int *result[] )
{
    int i,j;
	for( i = 0; i < height; ++i)
	{
		for( j = 0; j < width; ++j)
		{
			result[i][j] = (int)m[j * height + i];
		}
	}
}
void unpackMatrix_double2short(double m[], int height, int width, short *result[] )
{
    int i,j;
	for( i = 0; i < height; ++i)
	{
		for( j = 0; j < width; ++j)
		{
			result[i][j] = (short)m[j * height + i];
		}
	}
}

#endif  //SKIP_MEX


double** Alloc2d_double( int b, int c )
{
	double **p;
	int i;

	p = (double**)calloc( b, sizeof(double*) );
	assert(p);
	p[0] = (double*)calloc( b*c, sizeof(double) );
	for( i = 1; i < b; i++ )
	{
		p[i] = p[i-1] + c;
		assert( p[i] );
	}
	return p;
}

short** Alloc2d_short( int b, int c )
{
	short **p;
	int i;

	p = (short**)calloc( b, sizeof(short*) );
	assert(p);
	p[0] = (short*)calloc( b*c, sizeof(short) );
	for( i = 1; i < b; i++ )
	{
		p[i] = p[i - 1] + c;
		assert( p[i] );
	}
	return p;
}

int** Alloc2d_int( int b, int c )
{
	int **p;
	int i;

	p = (int**)calloc( b, sizeof(int*) );
	assert(p);
	p[0] = (int*)calloc( b*c, sizeof(int) );
	for( i = 1; i < b; i++ )
	{
		p[i] = p[i-1] + c;
		assert( p[i] );
	}
	return p;
}

UFLT16** Alloc2d_FLT16( int b, int c )
{
	UFLT16 **p;
	int i;

	p = (UFLT16**)calloc( b, sizeof(UFLT16*) );
	assert(p);
	p[0] = (UFLT16*)calloc( b*c, sizeof(UFLT16) );
	for( i = 1; i < b; i++ )
	{
		p[i] = p[i-1] + c;
		assert( p[i] );
	}
	return p;
}

void free2d_int( int **p )
{
	free( p[0] );
	free( p );
}

void free2d_double( double **p )
{
	free( p[0] );
	free( p );
}

void free2d_short( short **p )
{
	free( p[0] );
	free( p );
}

void free2d_FLT16( UFLT16 **p )
{
	free( p[0] );
	free( p );
}



void rotate_sign( short x[], short y[], int shift, int M )
{
	while( shift < 0 )
		shift += M;

	while( shift >= M )
		shift -= M;

	memcpy( y, x + shift, (M-shift)*sizeof(y[0]));
	memcpy( y+(M-shift), x, shift*sizeof(y[0]));
}

void rotate_data( double x[], double y[], int shift, int M )
{
	while( shift < 0 )
		shift += M;

	while( shift >= M )
		shift -= M;

	memcpy( y, x + shift, (M-shift)*sizeof(y[0]));
	memcpy( y+(M-shift), x, shift*sizeof(y[0]));
}

void rotate( void *x, void *y, int shift, int size, int M )
{
	int shift1;
	int shift2;
/*
	int shift1 = size * shift;
	int shift2 = size * (M-shift);
*/
	while( shift < 0 )
		shift += M;

	while( shift >= M )
		shift -= M;

	shift1 = size * shift;
	shift2 = size * (M-shift);

	memcpy( (char*)y, (char*)x + shift1, shift2 );
	memcpy( (char*)y + shift2, (char*)x, shift1 );
}

DEC_STATE* decod_open( int codec_id, int q_bits, int mh, int nh, int M )
{
    DEC_STATE* st;
    
	int N = nh * M;
	int R = mh * M;
	int BBsize = (mh > nh - mh) ? mh : nh - mh;


    st = (DEC_STATE*)calloc( 1, sizeof(DEC_STATE) );
    if( !st ) 
        return NULL;

	st->q_bits  = q_bits;
	st->nh      = nh;
	st->rh      = mh;
	st->m       = M;
    st->n       = N;
    st->codec_id = codec_id;

	switch( codec_id )
	{
	case BP_DEC: 
	case SP_DEC:
	case ASP_DEC: 
	case MS_DEC:
	case IMS_DEC:
	case IASP_DEC:
	case TASP_DEC: 
		st->bin_codec = 1;
		break;
	case FHT_DEC:
		st->bin_codec = 0;
		break;
	}

	st->q = st->bin_codec ? 1 : 1 << q_bits;

	if( st->bin_codec )
	{
		st->y = (double*)calloc(N, sizeof(st->y[0]));
		st->decword = (double*)calloc(N, sizeof(st->decword[0]) );
		if(!st->y || !st->decword )
			return NULL;

		st->hd = Alloc2d_short( mh, nh );
		if( st->hd==NULL)
			return NULL;
	}
	else
	{
		st->qy = Alloc2d_double( st->q, N );
		st->qdecword = Alloc2d_double( st->q, N );
		if(!st->qy || !st->qdecword )
			return NULL;

		st->hb = Alloc2d_short( mh, nh );
		if( st->hb==NULL)
			return NULL;

		st->hc = Alloc2d_short( mh, nh );
		if( st->hc==NULL)
			return NULL;

		st->qhard = (short*)calloc(N, sizeof(st->qhard[0]));
		if( st->qhard==NULL)
			return NULL;
	}

	st->syndr = (short*)calloc(R, sizeof(st->syndr[0]) );
	if( !st->syndr )
		return NULL;
    

	switch( codec_id )
	{
	case BP_DEC:
		st->bp_ZZ = Alloc2d_double( mh, N );
		if( st->bp_ZZ==NULL)
			return NULL;
		st->bp_BB = Alloc2d_short(BBsize, N );
		if( st->bp_BB==NULL)
			return NULL;
		st->bp_data = (double*)calloc( M, sizeof(st->bp_data[0]) );
		if( !st->bp_data )
		   return NULL;
		st->bp_sign = (short*)calloc( M, sizeof( st->bp_sign[0]) );
		if( !st->bp_sign )
			return NULL;
		st->bp_bs = (short*)calloc( R, sizeof(st->bp_bs[0]));
		if( !st->bp_bs )
			return NULL;
		st->bp_yd = (double*)calloc(N, sizeof(st->bp_yd[0]));
		if( !st->bp_yd )
			return NULL;
		st->bp_s = (double*)calloc( R, sizeof(st->bp_s[0]));
		if( !st->bp_s )
			return NULL;
		break;

	case SP_DEC:
		st->sp_ZZ = Alloc2d_double( mh, N );
		if( st->sp_ZZ==NULL)
			return NULL;
		st->sp_ZZ0 = Alloc2d_double( mh, M );
		if( st->sp_ZZ0==NULL)
			return NULL;
		st->sp_data = (double*)calloc( M, sizeof(st->sp_data[0]) );
		if( !st->sp_data )
			return NULL;
		st->sp_yd = (double*)calloc(N, sizeof(st->sp_yd[0]));
		if( !st->sp_yd )
			return NULL;
		st->sp_s = (double*)calloc( R, sizeof(st->sp_s[0]));
		if( !st->sp_s )
			return NULL;
		st->sp_AA = (double*)calloc( N, sizeof(st->sp_AA[0]));
		if( !st->sp_AA )
			return NULL;
		break;

	case ASP_DEC:
		st->asp_data0 = (double*)calloc(M, sizeof( st->asp_data0[0] ) );
		st->asp_data1 = (double*)calloc(M, sizeof( st->asp_data1[0] ) );
		if(!st->asp_data0 || !st->asp_data1 )
			return NULL;

		st->asp_p0 = (double*)calloc(M, sizeof( st->asp_p0[0] ) );
		st->asp_p1 = (double*)calloc(M, sizeof( st->asp_p1[0] ) );
		if(!st->asp_p0 || !st->asp_p1 )
			return NULL;

		st->asp_soft_out = (double*)calloc(N, sizeof(st->asp_soft_out[0]) );
		if( !st->asp_soft_out )
			return NULL;

		st->asp_posh = (int*)calloc(mh, sizeof(st->asp_posh[0]) );
		if( !st->asp_posh )
			return NULL;

		st->asp_rw = (int*)calloc(nh, sizeof(st->asp_rw[0]) );
		if( !st->asp_rw )
			return NULL;

		st->asp_hc_ri = Alloc2d_int( nh, 2 );
		if( st->asp_hc_ri==NULL)
			return NULL;

		st->asp_state = Alloc2d_double( nh, mh*M );
		if( st->asp_state==NULL)
			return NULL;

		st->asp_all_cw_2 = 0;

		break;

	case IASP_DEC:
#ifdef IASP_FIXED_POINT		
		st->iasp_y = (ui16*)calloc(N, sizeof(st->iasp_y[0]) );
		if( !st->iasp_y )
			return NULL;

		st->iasp_data0 = (ui16*)calloc(M, sizeof( st->iasp_data0[0] ) );
		st->iasp_data1 = (ui16*)calloc(M, sizeof( st->iasp_data1[0] ) );
		if(!st->iasp_data0 || !st->iasp_data1 )
			return NULL;

		st->iasp_state = (ui16**)Alloc2d_short( nh, mh*M );
		if( st->iasp_state==NULL)
			return NULL;

		st->iasp_p0 = (ui32*)calloc(M, sizeof( st->iasp_p0[0] ) );
		st->iasp_p1 = (ui32*)calloc(M, sizeof( st->iasp_p1[0] ) );
		if(!st->iasp_p0 || !st->iasp_p1 )
			return NULL;

#else
		st->iasp_y = (UFLT16*)calloc(N, sizeof(st->iasp_y[0]) );
		if( !st->iasp_y )
			return NULL;

		st->iasp_data0 = (UFLT16*)calloc(M, sizeof( st->iasp_data0[0] ) );
		st->iasp_data1 = (UFLT16*)calloc(M, sizeof( st->iasp_data1[0] ) );
		if(!st->iasp_data0 || !st->iasp_data1 )
			return NULL;

		st->iasp_state = (UFLT16**)Alloc2d_FLT16( nh, mh*M );
		if( st->iasp_state==NULL)
			return NULL;

		st->iasp_data = (ui16*)calloc(M, sizeof( st->iasp_data0[0] ) );
		if(!st->iasp_data)
			return NULL;

		st->iasp_p0 = (UFLT16*)calloc(M, sizeof( st->iasp_p0[0] ) );
		st->iasp_p1 = (UFLT16*)calloc(M, sizeof( st->iasp_p1[0] ) );
		if(!st->iasp_p0 || !st->iasp_p1 )
			return NULL;

#endif

		st->iasp_soft_out = (ui16*)calloc(N, sizeof(st->iasp_soft_out[0]) );
		if( !st->iasp_soft_out )
			return NULL;

		st->iasp_posh = (int*)calloc(mh, sizeof(st->iasp_posh[0]) );
		if( !st->iasp_posh )
			return NULL;

		st->iasp_rw = (int*)calloc(nh, sizeof(st->iasp_rw[0]) );
		if( !st->iasp_rw )
			return NULL;

		st->iasp_hc_ri = Alloc2d_int( nh, 2 );
		if( st->iasp_hc_ri==NULL)
			return NULL;


		st->iasp_all_cw_2 = 0;

		break;

	case MS_DEC:
		st->ms_soft = (MS_DATA*)calloc(N, sizeof(st->ms_soft[0]) );
		if( !st->ms_soft )
			return NULL;

		st->ms_BnNS = (short*)calloc(mh*N, sizeof( st->ms_BnNS[0] ) );
		if( !st->ms_BnNS )
			return NULL;

		st->ms_dcs = (MS_DEC_STATE*)calloc(R, sizeof( st->ms_dcs[0] ) );
		if( !st->ms_dcs )
			return NULL;

		st->ms_tmps = (MS_DEC_STATE*)calloc(M, sizeof( st->ms_tmps[0] ) );
		if( !st->ms_tmps )
			return NULL;

		st->ms_buffer = (MS_DATA*)calloc(M, sizeof(st->ms_buffer[0]) );
		if( !st->ms_buffer )
			return NULL;

		st->ms_rbuffer = (MS_DATA*)calloc(M, sizeof(st->ms_rbuffer[0]) );
		if( !st->ms_rbuffer )
			return NULL;

		st->ms_rsoft = (MS_DATA*)calloc(M, sizeof(st->ms_rsoft[0]) );
		if( !st->ms_rsoft )
			return NULL;
		break;

	case IMS_DEC:
		st->ims_y = (IMS_DATA*)calloc(N, sizeof(st->ims_y[0]) );
		if( !st->ims_y )
			return NULL;

		st->ims_soft = (IMS_DATA*)calloc(N, sizeof(st->ims_soft[0]) );
		if( !st->ims_soft )
			return NULL;

		st->ims_BnNS = (short*)calloc(mh*N, sizeof( st->ims_BnNS[0] ) );
		if( !st->ims_BnNS )
			return NULL;

		st->ims_dcs = (IMS_DEC_STATE*)calloc(R, sizeof( st->ims_dcs[0] ) );
		if( !st->ims_dcs )
			return NULL;

		st->ims_tmps = (IMS_DEC_STATE*)calloc(M, sizeof( st->ims_tmps[0] ) );
		if( !st->ims_tmps )
			return NULL;

		st->ims_buffer = (IMS_DATA*)calloc(M, sizeof(st->ims_buffer[0]) );
		if( !st->ims_buffer )
			return NULL;

		st->ims_rbuffer = (IMS_DATA*)calloc(M, sizeof(st->ims_rbuffer[0]) );
		if( !st->ims_rbuffer )
			return NULL;

		st->ims_rsoft = (IMS_DATA*)calloc(M, sizeof(st->ims_rsoft[0]) );
		if( !st->ims_rsoft )
			return NULL;
		break;

	case FHT_DEC:
		st->fht_pos = (short*)calloc(nh, sizeof(st->fht_pos[0]) );
		if( !st->fht_pos )
			return NULL;

		st->fht_rw = (int*)calloc(nh, sizeof(st->fht_rw[0]) );
		if( !st->fht_rw )
			return NULL;

		st->fht_list = (int*)calloc(st->q, sizeof(st->fht_list[0]) );
		if( !st->fht_list )
			return NULL;

		st->fht_ilist = (int*)calloc(st->q, sizeof(st->fht_ilist[0]) );
		if( !st->fht_ilist )
			return NULL;

		st->fht_buf = (short*)calloc(M, sizeof(st->fht_buf[0]) );
		if( !st->fht_buf )
			return NULL;

		st->fht_HAD = (double*)calloc(st->q*st->q, sizeof(st->fht_HAD[0]) );
		if( !st->fht_HAD )
			return NULL;

		st->fht_smask = (short*)calloc(N, sizeof(st->fht_smask[0]) );
		if( !st->fht_smask )
			return NULL;

		st->fht_hb_ci = Alloc2d_short(nh*M, 2 );
		if( st->fht_hb_ci==NULL)
			return NULL;

		st->fht_hb_cj = Alloc2d_short(st->nh, 2 );
		if( st->fht_hb_cj==NULL)
			return 0;

		st->fht_soft_out = Alloc2d_double(st->q, N );
		if( st->fht_soft_out==NULL)
			return NULL;

		st->fht_buf0 = Alloc2d_double(st->q, M );
		if( st->fht_buf0==NULL)
			return NULL;

		st->fht_buf1 = Alloc2d_double(st->q, M );
		if( st->fht_buf1==NULL)
			return NULL;
		break;

	case TASP_DEC:
		st->tasp_data0 = (double*)calloc(M, sizeof( st->tasp_data0[0] ) );
		if(!st->tasp_data0 )
			return NULL;

		st->tasp_tmp = (double*)calloc(nh, sizeof( st->tasp_tmp[0] ) );
		if(!st->tasp_tmp )
			return NULL;

//		st->tasp_p0 = (double*)calloc(M, sizeof( st->tasp_p0[0] ) );
//		st->tasp_p1 = (double*)calloc(M, sizeof( st->tasp_p1[0] ) );
//		if(!st->tasp_p0 || !st->tasp_p1 )
//			return NULL;

		st->tasp_soft_out = (double*)calloc(N, sizeof(st->tasp_soft_out[0]) );
		if( !st->tasp_soft_out )
			return NULL;

//		st->tasp_posh = (int*)calloc(mh, sizeof(st->tasp_posh[0]) );
//		if( !st->tasp_posh )
//			return NULL;

//		st->tasp_rw = (int*)calloc(nh, sizeof(st->tasp_rw[0]) );
//		if( !st->tasp_rw )
//			return NULL;

//		st->tasp_hc_ri = Alloc2d_int( nh, 2 );
//		if( st->tasp_hc_ri==NULL)
//			return NULL;

//		st->tasp_state = Alloc2d_double( nh, mh*M );
		st->tasp_state = Alloc2d_double( mh*M, nh*M );
		if( st->tasp_state==NULL)
			return NULL;

		//st->tasp_all_cw_2 = 0;
		break;

	default: return NULL;
	}

    return st;
    
}



int find_row_weight( short *hb[], int rh, int nh, int rw[] )
{
	int i, j;
	int max_w;

	max_w = 0;
	for( i = 0; i < rh; i++ )
	{
		int cnt = 0;
		for( j = 0; j < nh; j++ )
			cnt += hb[i][j] != -1;
		rw[i] = cnt;

		if( max_w < cnt )
			max_w = cnt;
	}
	return max_w;
}


int find_column_weight( short *hb[], int rh, int nh, short *hb_ri[] )
{
	int i, j;

	int all_cw2 = 1;
	for( i = 0; i < nh; i++ )
	{
		int cnt = 0;
		for( j = 0; j < rh; j++ )
			cnt += hb[j][i] != -1;

		if( cnt != 2 )
			all_cw2 = 0;
	}

	for( j = 0; j < nh; j++ )
	{
		int cnt = 0;

		for( i = 0; i < rh; i++ )
			if( hb[i][j] != -1 )
				hb_ri[j][cnt++] = i;
	}

	return all_cw2;
}


int find_list_of_symbols( int q, short *hc[], int rh, int nh, int list[], int ilist[] )
{
	int i, j;

	for( i = 0; i < q; i++ )
		list[i] = 0;

	for( i = 0; i < rh; i++ )
	{
		for( j = 0; j < nh; j++ )
		{
			int val = hc[i][j];

			if( val != -1 )
			{
				list[val] = 1;
			}
		}
	}

	j = 0;
	for( i = 0; i < q; i++ )
	{
		if( list[i] )
		{
			ilist[j] = i;
			list[i]  = j;	
			j++;
		}
	}

	return j;
}

int prepare_mul_tab
(
	 short *hc[], 
	 int rh, 
	 int nh, 
	 int list[], 
	 short *hc_rl[], 
	 int q_bits, 
	 short gf2log[], 
	 short gf2alog[], 
	 short *T[], 
	 short *TI[],
	 int ilist[],
	 int ilist_size
	 )
{
	int i, j;
	int q = 1 << q_bits;
	int t[2];
	int pol;
	int unused_M;
	
	for( i = 0; i < rh; i++ )
	{
		int cnt = 0;
		for( j = 0; j < nh; j++ )
		{
			int val = hc[i][j];

			if( val != -1 )
			{
				hc_rl[i][cnt++] = list[val];
			}
		}
		hc_rl[i][cnt] = -1;
	}


	t[0] = q_bits;
	t[1] = 1;

	if( pol_bank( t, &pol, &unused_M ) == 0 )
		return 0;

	gf2( q_bits, pol, gf2log, gf2alog );

	p2table( T, TI, gf2log, gf2alog, ilist, q_bits, ilist_size );

	return 1;

}

void prepare_cj( short **hb, int rh, int nh, short **cj )
{
	int i, j, k, cnt, cnum;

	for( i = 0; i < nh; i++ )
	{
		cnum = 0;
		for( j = 0; j < rh; j++ )
		{
			if( hb[j][i] != -1 )
			{
				cnt = 0;
				for( k = 0; k < i; k++ )
				{
					if( hb[j][k] != -1 )
						cnt++;
				}
				cj[i][cnum++] = cnt;
			}
		}
	}
}

void prepare_HAD( int c, double *HAD )
{
	double buffer[1024];
	int i, j;
	int q = 1 << c;

	for( i = 0; i < q; i++ )
	{
		for( j = 0; j < q; j++ )
			buffer[j] = 0;
		buffer[i] = 1;

		hadamar( buffer, &HAD[i * q], q, c );
	}
}

void prepare_row_indexes( short **hc, int rh, int nh, short **hc_ri )
{
	int i, j;

	for( i = 0; i < rh; i++ )
	{
		int cnt = 0;
		for( j = 0; j < nh; j++ )
		{
			if( hc[i][j] != -1 )
				hc_ri[i][cnt++] = j;
		}
	}
}


int decod_init( void *state )
{
	int i, j;
	DEC_STATE* st = (DEC_STATE*)state;
	short *gf2log;
	short *gf2alog;

	if( !st ) 
		return 1;


	switch( st->codec_id )
	{
	case BP_DEC:
		break;

	case SP_DEC:
		break;

	case ASP_DEC:
		st->asp_all_cw_2 = 1;

		for( i = 0; i < st->nh; i++ )
		{
			int cnt = 0;

			for( j = 0; j < st->rh; j++ )
				cnt +=  st->hd[j][i] != -1;

			if( cnt != 2 )
			{
				st->asp_all_cw_2 = 0;
				break;
			}
		}


		if( st->asp_all_cw_2 )
		{
			for( i = 0; i < st->nh; i++ )
			{
				int cnt = 0;

				for( j = 0; j < st->rh; j++ )
				{
					if( st->hd[j][i] != -1 )
					{
						st->asp_hc_ri[i][cnt] = j;
						cnt++;
					}
				}
			}
		}
		break;

	case IASP_DEC:
		st->iasp_all_cw_2 = 1;

		for( i = 0; i < st->nh; i++ )
		{
			int cnt = 0;

			for( j = 0; j < st->rh; j++ )
				cnt +=  st->hd[j][i] != -1;

			if( cnt != 2 )
			{
				st->iasp_all_cw_2 = 0;
				break;
			}
		}


		if( st->iasp_all_cw_2 )
		{
			for( i = 0; i < st->nh; i++ )
			{
				int cnt = 0;

				for( j = 0; j < st->rh; j++ )
				{
					if( st->hd[j][i] != -1 )
					{
						st->iasp_hc_ri[i][cnt] = j;
						cnt++;
					}
				}
			}
		}
		break;

	case MS_DEC:
		break;

	case IMS_DEC:
		break;

	case FHT_DEC:
		st->max_rw = find_row_weight( st->hb, st->rh, st->nh, st->fht_rw );

		st->fht_hb_ri = Alloc2d_short( st->rh, st->max_rw );
		if( st->fht_hb_ri==NULL)
			return 0;

		prepare_row_indexes( st->hb, st->rh, st->nh, st->fht_hb_ri );

		st->fht_all_cw_2 = find_column_weight( st->hb, st->rh, st->nh, st->fht_hb_ci );


		if( st->fht_all_cw_2 == 0 )
			return 0;

		st->fht_soft_in = Alloc2d_double(st->q, st->rh * st->m * st->max_rw );
		if( st->fht_soft_in==NULL)
			return 0;

		st->fht_soft_outs = Alloc2d_double(st->q, st->rh * st->m * st->max_rw );
		if( st->fht_soft_outs==NULL)
			return 0;

		st->fht_hc_rl = Alloc2d_short(st->rh, st->max_rw+1 );
		if( st->fht_hc_rl==NULL)
			return 0;

		prepare_cj( st->hb, st->rh, st->nh, st->fht_hb_cj );

		st->fht_ilist_size = find_list_of_symbols(st->q, st->hc, st->rh, st->nh, st->fht_list, st->fht_ilist );
		prepare_HAD( st->q_bits, st->fht_HAD );

#ifdef ORIG_TABLES
		st->fht_t = Alloc2d_short( st->q, st->fht_ilist_size );
		st->fht_ti = Alloc2d_short( st->q, st->fht_ilist_size );
#else
		st->fht_t = Alloc2d_short( st->fht_ilist_size, st->q );
		st->fht_ti = Alloc2d_short( st->fht_ilist_size, st->q );
#endif
		if( st->fht_t==NULL || st->fht_ti==NULL)
			return 0;

		gf2log  = (short*)calloc( st->q, sizeof(short) );
		gf2alog = (short*)calloc( st->q, sizeof(short) );

		prepare_mul_tab( st->hc, st->rh, st->nh, st->fht_list, st->fht_hc_rl, st->q_bits, gf2log, gf2alog, st->fht_t, st->fht_ti, st->fht_ilist, st->fht_ilist_size );

		free( gf2log );
		free( gf2alog );

		break;

	case TASP_DEC:
/*
		st->tasp_all_cw_2 = 1;

		for( i = 0; i < st->nh; i++ )
		{
			int cnt = 0;

			for( j = 0; j < st->rh; j++ )
				cnt +=  st->hd[j][i] != -1;

			if( cnt != 2 )
			{
				st->tasp_all_cw_2 = 0;
				break;
			}
		}


		if( st->asp_all_cw_2 )
		{
			for( i = 0; i < st->nh; i++ )
			{
				int cnt = 0;

				for( j = 0; j < st->rh; j++ )
				{
					if( st->hd[j][i] != -1 )
					{
						st->asp_hc_ri[i][cnt] = j;
						cnt++;
					}
				}
			}
		}
*/
		break;



	default:;
	}

	return 1;
}


void decod_close( DEC_STATE* st )
{
	int nh = st->nh;
	int mh = st->rh;
	int N = st->n;
	int M = st->m;
	int q_bits = st->q_bits;
	int q = st->q;
	int codec_id = st->codec_id;
	int R = mh * M;

	if( st->bin_codec )
	{
		if(st->hd)		{ free2d_short( st->hd );          st->hd       = NULL; }
		if(st->y)		{ free(st->y);                     st->y        = NULL; }
		if(st->decword) { free(st->decword);               st->decword  = NULL; }
	}
	else
	{
		if(st->hb)		 { free2d_short( st->hb );         st->hb        = NULL; }
		if(st->hc)		 { free2d_short( st->hc );         st->hc        = NULL; }
		if(st->qy)		 { free2d_double(st->qy );		   st->qy        = NULL; }
		if(st->qdecword) { free2d_double(st->qdecword);    st->qdecword  = NULL; }
	}

	if(st->syndr)	{ free( st->syndr ); st->syndr = NULL; }

	switch( codec_id )
	{
	case BP_DEC:
		if( st->bp_ZZ )    { free2d_double( st->bp_ZZ );    st->bp_ZZ   = NULL; }
		if( st->bp_BB )    { free2d_short( st->bp_BB  );    st->bp_BB   = NULL; }
		if( st->bp_data )  { free( st->bp_data );           st->bp_data = NULL; }
		if( st->bp_sign )  { free( st->bp_sign );           st->bp_sign = NULL; }
		if( st->bp_bs )    { free( st->bp_bs );             st->bp_bs   = NULL; }
		if( st->bp_yd )    { free( st->bp_yd );             st->bp_yd   = NULL; }
		if( st->bp_s )     { free( st->bp_s );              st->bp_s    = NULL; }
		break;

	case SP_DEC:
		if( st->sp_ZZ )   { free2d_double( st->sp_ZZ );     st->sp_ZZ   = NULL; }
		if( st->sp_ZZ0 )  { free2d_double( st->sp_ZZ0 );    st->sp_ZZ0  = NULL; }
		if( st->sp_data ) { free( st->sp_data );            st->sp_data = NULL; }
		if( st->sp_yd )   { free( st->sp_yd );              st->sp_yd   = NULL; }
		if( st->sp_s )    { free( st->sp_s );               st->sp_s    = NULL; }
		if( st->sp_AA )   { free( st->sp_AA );              st->sp_AA   = NULL; }
		break;

	case ASP_DEC:
		if( st->asp_data0 )     { free(st->asp_data0);            st->asp_data0    = NULL; }
		if( st->asp_data1 )     { free(st->asp_data1);            st->asp_data1    = NULL; }
		if( st->asp_p0 )        { free(st->asp_p0);               st->asp_p0       = NULL; }
		if( st->asp_p1 )        { free(st->asp_p1);               st->asp_p1       = NULL; }
		if( st->asp_soft_out )  { free(st->asp_soft_out);         st->asp_soft_out = NULL; }
		if( st->asp_posh )      { free(st->asp_posh);             st->asp_posh     = NULL; }
		if( st->asp_rw )        { free(st->asp_rw);               st->asp_rw       = NULL; }
		if( st->asp_state )     { free2d_double( st->asp_state ); st->asp_state    = NULL; }
		if( st->asp_hc_ri )     { free2d_int( st->asp_hc_ri );    st->asp_hc_ri    = NULL; }
		break;

	case IASP_DEC:
		if( st->iasp_y )         { free(st->iasp_y);                         st->iasp_y        = NULL; }
#ifndef IASP_FIXED_POINT
		if( st->iasp_data )      { free(st->iasp_data);                      st->iasp_data     = NULL; }
#endif
		if( st->iasp_data0 )     { free(st->iasp_data0);                     st->iasp_data0    = NULL; }
		if( st->iasp_data1 )     { free(st->iasp_data1);                     st->iasp_data1    = NULL; }
		if( st->iasp_p0 )        { free(st->iasp_p0);                        st->iasp_p0       = NULL; }
		if( st->iasp_p1 )        { free(st->iasp_p1);                        st->iasp_p1       = NULL; }
		if( st->iasp_soft_out )  { free(st->iasp_soft_out);                  st->iasp_soft_out = NULL; }
		if( st->iasp_posh )      { free(st->iasp_posh);                      st->iasp_posh     = NULL; }
		if( st->iasp_rw )        { free(st->iasp_rw);                        st->iasp_rw       = NULL; }
#ifdef IASP_FIXED_POINT		
		if( st->iasp_state )     { free2d_short( (short**)st->iasp_state );  st->iasp_state    = NULL; }
#else
		if( st->iasp_state )     { free2d_FLT16( (UFLT16**)st->iasp_state);  st->iasp_state    = NULL; }
#endif
		if( st->iasp_hc_ri )      { free2d_int( st->iasp_hc_ri );            st->iasp_hc_ri    = NULL; }
		break;

	case MS_DEC:
		if( st->ms_soft )    { free(st->ms_soft);    st->ms_soft    = NULL; }
		if( st->ms_BnNS )    { free(st->ms_BnNS);    st->ms_BnNS    = NULL; }
		if( st->ms_dcs )     { free(st->ms_dcs);     st->ms_dcs     = NULL; }
		if( st->ms_tmps )    { free(st->ms_tmps);    st->ms_tmps    = NULL; }
		if( st->ms_buffer )  { free(st->ms_buffer);  st->ms_buffer  = NULL; }
		if( st->ms_rbuffer ) { free(st->ms_rbuffer); st->ms_rbuffer = NULL; }
		if( st->ms_rsoft )   { free(st->ms_rsoft);   st->ms_rsoft   = NULL; }
		break;

	case IMS_DEC:
		if( st->ims_y )       { free(st->ims_y);       st->ims_y       = NULL; }
		if( st->ims_soft )    { free(st->ims_soft);    st->ims_soft    = NULL; }
		if( st->ims_BnNS )    { free(st->ims_BnNS);    st->ims_BnNS    = NULL; }
		if( st->ims_dcs )     { free(st->ims_dcs);     st->ims_dcs     = NULL; }
		if( st->ims_tmps )    { free(st->ims_tmps);    st->ims_tmps    = NULL; }
		if( st->ims_buffer )  { free(st->ims_buffer);  st->ims_buffer  = NULL; }
		if( st->ims_rbuffer ) { free(st->ims_rbuffer); st->ims_rbuffer = NULL; }
		if( st->ims_rsoft )   { free(st->ims_rsoft);   st->ims_rsoft   = NULL; }
		break;

	case FHT_DEC:
		if( st->qhard )		   { free(st->qhard);                   st->qhard        = NULL; }
		if( st->fht_buf )	   { free(st->fht_buf);                 st->fht_buf      = NULL; }
		if( st->fht_HAD )	   { free(st->fht_HAD);                 st->fht_HAD      = NULL; }
		if( st->fht_smask )	   { free(st->fht_smask);               st->fht_smask    = NULL; }
		if( st->fht_pos )      { free(st->fht_pos);                 st->fht_pos      = NULL; }
		if( st->fht_rw )       { free(st->fht_rw);                  st->fht_rw       = NULL; }
		if( st->fht_list )     { free(st->fht_list);                st->fht_list     = NULL; }
		if( st->fht_ilist )    { free(st->fht_ilist);               st->fht_ilist    = NULL; }
		if( st->fht_hb_ci )    { free2d_short( st->fht_hb_ci );     st->fht_hb_ci    = NULL; }
		if( st->fht_hb_ri )    { free2d_short( st->fht_hb_ri );     st->fht_hb_ri    = NULL; }
		if( st->fht_hc_rl )    { free2d_short( st->fht_hc_rl );     st->fht_hc_rl    = NULL; }
		if( st->fht_hb_cj )    { free2d_short( st->fht_hb_cj );     st->fht_hb_cj    = NULL; }
		if( st->fht_t )        { free2d_short( st->fht_t );         st->fht_t        = NULL; }
		if( st->fht_ti )       { free2d_short( st->fht_ti );        st->fht_ti       = NULL; }
		if( st->fht_soft_in )  { free2d_double( st->fht_soft_in );  st->fht_soft_in  = NULL; }
		if( st->fht_soft_outs ){ free2d_double( st->fht_soft_outs );st->fht_soft_outs= NULL; }
		if( st->fht_soft_out ) { free2d_double( st->fht_soft_out ); st->fht_soft_out = NULL; }
		if( st->fht_buf0 )     { free2d_double( st->fht_buf0 );     st->fht_buf0   = NULL; }
		if( st->fht_buf1 )     { free2d_double( st->fht_buf1 );     st->fht_buf1   = NULL; }
		break;

	case TASP_DEC:
		if( st->tasp_data0 )     { free(st->tasp_data0);            st->tasp_data0    = NULL; }
		if( st->tasp_tmp )     { free(st->tasp_tmp);		        st->tasp_tmp      = NULL; }
//		if( st->tasp_p0 )        { free(st->tasp_p0);               st->tasp_p0       = NULL; }
//		if( st->tasp_p1 )        { free(st->tasp_p1);               st->tasp_p1       = NULL; }
		if( st->tasp_soft_out )  { free(st->tasp_soft_out);         st->tasp_soft_out = NULL; }
//		if( st->tasp_posh )      { free(st->tasp_posh);             st->tasp_posh     = NULL; }
//		if( st->tasp_rw )        { free(st->tasp_rw);               st->tasp_rw       = NULL; }
		if( st->tasp_state )     { free2d_double( st->tasp_state ); st->tasp_state    = NULL; }
//		if( st->tasp_hc_ri )     { free2d_int( st->tasp_hc_ri );    st->tasp_hc_ri    = NULL; }
		break;

	default:;
	}

    free( st );
}

int bp_decod_qc_lm( DEC_STATE* st, double soft[], double decword[], int maxiter, int decision )
{ 
	// BP decoding (Gallager)
	// soft    - channel output 
	// maxiter - maximum number of iterations	
	int nh    = st->nh;
	int rh    = st->rh;
	int m     = st->m;
	int Ncode = nh * m;
	int Rcode = rh * m;
    int Kcode = Ncode - Rcode;
	int i, j, k, synd, iter=0;
    short b;
	double A;
    double x, Eps=2.0e-5;
	double *data = st->bp_data;
	short  *sign = st->bp_sign;
	short  *bs   = st->bp_bs;
	double *yd   = st->bp_yd;
	double *s    = st->bp_s;
	double **ZZ  = st->bp_ZZ;
	short  **BB  = st->bp_BB;
	short  **hd  = st->hd;

    for( i = 0; i < rh; i++ )
		for( j = 0; j < Ncode; j++ )
			ZZ[i][j] = 0.0;
    
	
	for( i = 0; i < Ncode; i++ ) 
		yd[i] = soft[i] = maxd(mind(soft[i],  INPUT_LIMIT),-INPUT_LIMIT); 


	// check input codeword
	for( j = 0; j < rh; j++ )
	{
		for( i = 0; i < nh; i++ ) 
		{
			int pos_r = j * m;
			int pos_n = i * m;
//			int circ = hd[j*stride + i];
			int circ = hd[j][i];

			if( circ != - 1 )
			{
				rotate_data( &soft[pos_n], data, circ, m );

				for( k = 0; k < m; k++ )
					st->syndr[pos_r + k] ^= data[k] < 0;
			}
		}
	}

	synd = 0;
	for( i = 0; i < Rcode; i++ )
		synd |= st->syndr[i];

	if (!synd) 
	{
		if( decision )
		{
			for( i = 0; i < Ncode; i++ )
				decword[i] = soft[i];
		}
		else
		{
			for( i = 0; i < Ncode; i++ )
				decword[i] = (soft[i] < 0 ) ? 1 : 0;
		}

		return 0;// no errors detected
	}

	while (iter < maxiter)
	{
		//for( i = 0; i < Rcode; i++ ) st->syndr[i] = 0;
		//for( i = 0; i < Rcode; i++ ) st->bs[i] = 0;
		//for( i = 0; i < Rcode; i++ ) st->s[i] = 0.0;
		memset( st->syndr, 0, Rcode*sizeof(st->syndr[0]));
		memset( bs, 0, Rcode*sizeof(bs[0]));
		memset(  s, 0, Rcode*sizeof(s[0]) );

		for( i = 0; i < nh; i++ )    
		{
			for( j = 0; j < rh; j++ ) 
			{
				int pos_r = j * m;
				int pos_n = i * m;
//				int circ = hd[j*stride + i];
				int circ = hd[j][i];
		
				if( circ != -1 )
				{
					// Variable-node activation step
					for( k = 0; k < m; k++ )
					{
						A = exp(soft[pos_n + k] - ZZ[j][pos_n + k]);
						x = log(absd((A-1)/(A+1)));
						
						BB[j][pos_n + k] = A < 1; 
#ifdef BP_USE_EPS
						ZZ[j][pos_n + k] = absd(x) < Eps ? (x<=0 ? -Eps : Eps) : x;
#else
						ZZ[j][pos_n + k] = x;
#endif
					}

					// Check-node activation step
					rotate_data( &ZZ[j][ pos_n], data, circ, m );
					rotate_sign( &BB[j][ pos_n], sign, circ, m );


					for( k = 0; k < m; k++ )
					{
						s [pos_r + k] += data[k];
						bs[pos_r + k] ^= sign[k];
					}
				}
			}
		}


		//for( i = 0; i < Ncode; i++ )
		//	soft[i] = st->yd[i]; 
		memcpy( soft, yd, Ncode*sizeof(soft[0]));

		for( i = 0; i < nh; i++ )    
		{
			for( j = 0; j < rh; j++ )
			{
				int pos_n = i * m;
				int pos_r = j * m;
				int circ = hd[j][i];

				if( circ != -1 )
				{
					rotate_data( &s[pos_r],  data, m - circ, m );
					rotate_sign( &bs[pos_r], sign, m - circ, m );

					for( k = 0; k < m; k++ )
					{
						A = exp( data[k] - ZZ[j][pos_n + k] );

						b = sign[k] ^ BB[j][pos_n + k];

						A = (1-2*b)*log((1+A)/(1-A));

						ZZ[j][pos_n + k] = maxd(mind(A, 19.07),-19.07);
					}

					for( k = 0; k < m; k++ )
						soft[pos_n + k] += ZZ[j][pos_n + k];
				}
			}
		}


		for( i = 0; i < nh; i++ ) 
		{
			for( j = 0; j < rh; j++ )
			{
				int pos_r = j * m;
				int pos_n = i * m;
//				int circ = hd[j*stride + i];
				int circ = hd[j][i];

				if( circ != - 1 )
				{
					rotate_data( &soft[pos_n], data, circ, m );

					for( k = 0; k < m; k++ )
						st->syndr[pos_r + k] ^= data[k] < 0;
				}
			}
		}

		synd = 0;
		for( i = 0; i < Rcode; i++ )
			synd |= st->syndr[i];

		if (!synd) 
		{   
			if( decision )
			{
				for( i = 0; i < Ncode; i++ )
					decword[i] = soft[i];
			}
			else
			{
				for( i = 0; i < Ncode; i++ )
					decword[i] = (soft[i] < 0 ) ? 1 : 0;
			}
			iter++;
			return iter; // errors corrected
		}

		iter++;
	}  //while

	if( decision )
	{
		for( i = 0; i < Ncode; i++ )
			decword[i] = soft[i];
	}
	else
	{
		for( i = 0; i < Ncode; i++ )
			decword[i] = (soft[i] < 0 ) ? 1 : 0;
	}

        
	return -iter;
}   // bp


int sum_prod_decod_qc_lm( DEC_STATE* st, double soft[], double decword[], int maxiter, int decision )    
{
	int nh = st->nh;
	int rh = st->rh;
	int m  = st->m;

	short **hd   = st->hd;
	double *yd   = st->sp_yd;
	double *s    = st->sp_s;
	double **ZZ  = st->sp_ZZ;
	double **ZZ0 = st->sp_ZZ0;
	short *syndr = st->syndr;
	double *data = st->sp_data;
	double *AA   = st->sp_AA;

	int i, j, k, synd, iter=0;

	int Ncode  = nh * m;
	int Rcode  = rh * m;

#if 0
	for( i = 0; i < Ncode; i++ )
		st->yd[i] = soft[i] = pow( 2.0, soft[i] );
#else
	for( i = 0; i < Ncode; i++ )
	{
		double y = maxd( mind( soft[i],  INPUT_LIMIT ), -INPUT_LIMIT );
		yd[i] = soft[i] = exp( y ); 
	}
#endif

	for( i = 0; i < Rcode; i++ )
		s[i] = 0;

	for( i = 0; i < rh; i++ )
		for( j = 0; j < Ncode; j++ )
			ZZ[i][j] = 1.0;


	// check input codeword

	for( i = 0; i < Rcode; i++ )
		syndr[i] = 0;

	for( j = 0; j < rh; j++ )
	{
		for( i = 0; i < nh; i++ ) 
		{
			int pos_r = j * m;
			int pos_n = i * m;
			int circ = hd[j][i];

			if( circ != - 1 )
			{
				rotate( &soft[pos_n], data, circ, sizeof(data[0]), m );

				for( k = 0; k < m; k++ )
					syndr[pos_r + k] ^= data[k] < 0.5;
			}
		}
	}

	synd = 0;
	for( i = 0; i < Rcode; i++ )
		synd |= syndr[i];

	if (!synd) 
	{
		if( decision )
		{
			for( i = 0; i < Ncode; i++ ) 
				decword[i] = soft[i];
		}
		else
		{
			for( i = 0; i < Ncode; i++ ) 
				decword[i] = soft[i] < 0.5;
		}
		return 0;// no errors detected
	}



	while( iter < maxiter )
	{

		for( i = 0; i < Rcode; i++ ) syndr[i] = 0;
		for( j = 0; j < Rcode; j++ ) s[j]     = 1.0;
		for( i = 0; i < Ncode; i++ ) soft[i]  = yd[i]; 

		for( i = 0; i < nh; i += 1 )    
		{
			int pos_n = i * m;
#if 1
			for( j = 0; j < rh; j++ )
			{
				int pos_r = j * m;
				int circ = hd[j][i];

				if( circ != -1 )
				{
					int jj;
				
					// Variable-node activation step
					memcpy( AA, &yd[pos_n], m*sizeof(AA[0]) );

					for( jj = 0; jj < rh; jj++ )
					{
						int circ = hd[jj][i];

						if( jj == j )
							continue;

						if( circ != -1 )
						{
							for( k = 0; k < m; k++ )
								AA[k] *= ZZ[jj][pos_n + k];		// product of all elements in column except of j-th
						}
					}

					for( k = 0; k < m; k++ )
						ZZ0[j][k]  = (AA[k]-1)/(AA[k]+1);

					// Check-node activation step
					rotate( &ZZ0[j][0], data, circ, sizeof(data[0]), m );

					for( k = 0; k < m; k++ )
						s[pos_r + k] *= data[k];
				}
			}

			for( j = 0; j < rh; j++ )
			{
				int circ = hd[j][i];

				if( circ != -1 )
					memcpy( &ZZ[j][pos_n], &ZZ0[j][0], m*sizeof(ZZ[0][0]) );
			}

#else
			for( k = 0; k < m; k++ )
				AA[k] = yd[pos_n + k];

			// Variable-node activation step
			for( j = 0; j < rh; j++ )
			{
				int circ = hd[j][i];

				if( circ != -1 )
				{
					for( k = 0; k < m; k++ )
						AA[k] *= ZZ[j][pos_n + k];				// product of all elements in column except of j-th
				}
			}

			for( j = 0; j < rh; j++ )
			{
				int pos_r = j * m;
				int circ = hd[j][i];

				if( circ != -1 )
				{
					// Variable-node activation step
					for( k = 0; k < m; k++ )
					{
						double A = AA[k] / ZZ[j][pos_n + k];	// product of all elements in column except of j-th
						ZZ[j][pos_n + k]  = (A-1)/(A+1);
					}

					// Check-node activation step
					rotate( &ZZ[j][pos_n], data, circ, sizeof(data[0]), m );

					for( k = 0; k < m; k++ )
						s[pos_r + k] *= data[k];
				}
			}
#endif
		}


		for( i = 0; i < nh; i++ )    
		{
			for( j = 0; j < rh; j++ )
			{
				int circ = hd[j][i];
				int pos_r = j * m;
				int pos_n = i * m;

				if( circ != -1 )
				{
					rotate( &s[pos_r], data, m-circ, sizeof(data[0]), m );
				
					for( k = 0; k < m; k++ )
					{
						double A = data[k] / ZZ[j][pos_n + k];

						A = (1+A)/(1-A);
						A = maxd( mind( A,  1.9e+8 ), -5.2e-9 );
//						A = maxd( mind( A,  1.9e+8 ), 5.2e-9 );

						ZZ[j][pos_n + k] = A;
						soft[pos_n + k] *= A;
					}
				}
			}

			for( j = 0; j < rh; j++ )
			{
				int pos_r = j * m;
				int pos_n = i * m;
				int circ = hd[j][i];

				if( circ != - 1 )
				{
					rotate( &soft[pos_n], data, circ, sizeof(data[0]), m );

					for( k = 0; k < m; k++ )
						syndr[pos_r + k] ^= data[k] < 0.5;
				}
			}

		}


		synd = 0;
		for( i = 0; i < Rcode; i++ )
			synd |= syndr[i];

		if (!synd) 
		{   
			if( decision )
			{
				for( i = 0; i < Ncode; i++ ) 
					decword[i] = soft[i];
			}
			else
			{
				for( i = 0; i < Ncode; i++ ) 
					decword[i] = soft[i] < 0.5;
			}

			iter++;
			//mexPrintf("-------> iter = %d\n", iter );
			return iter; // errors corrected

		}

		iter++;
	}  //while

	if( decision )
	{
		for( i = 0; i < Ncode; i++ ) 
			decword[i] = soft[i];
	}
	else
	{
		for( i = 0; i < Ncode; i++ ) 
			decword[i] = soft[i] < 0.5;
	}

	return -iter;
} 



#define ROW_WEIGHT_MAX 1024

void map_bin( double soft[], int rw, int step)
{
	int i;
	double SF[ROW_WEIGHT_MAX];
	double SB[ROW_WEIGHT_MAX];
	double P[ROW_WEIGHT_MAX];

//	for( i = 0; i < rw; i++ )	soft_out[i] = 0;
/*
	for( i = 0; i < rw; i++ )	SF[i] = 0;
	for( i = 0; i < rw; i++ )	SB[i] = 1;
*/
	SB[0]    = 1.0;
	SF[rw-1] = 0.0;

	for( i = 0; i < rw; i++ )	P[i] = 1 - 2 * (double)soft[i*step];

	// Alpha
	SF[0] = P[0]; 
	for( i = 1; i < rw-1; i++ )
		SF[i] = P[i] * SF[i-1];

	// Beta
	SB[rw-1] = P[rw-1]; 
	for( i = rw-2; i > 0; i-- )
		SB[i] = P[i] * SB[i+1]; 

	// Sigma
	soft[0*step] = (1-SB[1])/2;

	for( i = 1; i < rw-1; i++ ) //loop over symbols in the current check
	{
		double Z = SF[i-1] * SB[i+1];
	
		soft[i*step] = (1-Z)/2;
	}
	soft[(rw-1)*step] = (1 - SF[rw-2])/2;
}






void imap_bin( ui16 soft[], int rw, int step)
{
	int i;
	i16 SF[ROW_WEIGHT_MAX];
	i16 SB[ROW_WEIGHT_MAX];
	i16  P[ROW_WEIGHT_MAX];

//	SB[0]    = ONE_SOFT;
//	SF[rw-1] = 0;

	for( i = 0; i < rw; i++ )
		P[i] = ONE_SOFT - 2 * soft[i*step];

	// Alpha
	SF[0] = P[0]; 
	for( i = 1; i < rw-1; i++ )
		SF[i] = div_power2r((int)P[i] * SF[i-1], SOFT_FPP);

	// Beta
	SB[rw-1] = P[rw-1]; 
	for( i = rw-2; i > 0; i-- )
		SB[i] = div_power2r((int)P[i] * SB[i+1], SOFT_FPP); 

	// Sigma
	soft[0*step] = div_power2r(ONE_SOFT - SB[1], 1);
	soft[0*step] = maxi( soft[0*step], 1 );

	for( i = 1; i < rw-1; i++ ) //loop over symbols in the current check
	{
		int Z = div_power2r((int)SF[i-1] * SB[i+1], SOFT_FPP);

		soft[i*step] = div_power2r(ONE_SOFT - Z, 1);
		soft[i*step] = maxi( soft[i*step], 1 );
	}
	soft[(rw-1)*step] = div_power2r(ONE_SOFT - SF[rw-2], 1);
	soft[(rw-1)*step] = maxi( soft[(rw-1)*step], 1 );
}


int check_syndrome( short *syndr, int r, short **hd, int rh, int nh, int m, double *soft, double *buf, double thr )
{
	int synd;
	int i, j, k;

	for( i = 0; i < r; i++ )
		syndr[i] = 0;

	for( j = 0; j < rh; j++ )
	{
		for( i = 0; i < nh; i++ ) 
		{
			int circ = hd[j][i];

			if( circ != - 1 )
			{
				int pos_r = j * m;
				int pos_n = i * m;

				rotate( &soft[pos_n], buf, circ, sizeof(soft[0]), m );

				for( k = 0; k < m; k++ )
					syndr[pos_r + k] ^= buf[k] > thr;
			}
		}
	}

	synd = 0;
	for( i = 0; i < r; i++ )
		synd |= syndr[i];

	return synd;
}

void make_output( int mode, int n, double *outword, double *soft, double thr )
{
	int i;

	if( mode )
	{
		for( i = 0; i < n; i++ )
			outword[i] = soft[i];
	}
	else
	{
		for( i = 0; i < n; i++ )
			outword[i] = soft[i] > thr;
	}
}

int sum_prod_gf2_decod_qc_lm( DEC_STATE* st, double soft[], double decword[], int maxsteps, int decision )
{ 
	int i, j, k;
	int steps;
	short **hd = st->hd;
	double **state = st->asp_state;
	double *soft_out = st->asp_soft_out;
	short *syndr    = st->syndr;
	double *data0 = st->asp_data0;
	double *data1 = st->asp_data1;
	double *P1 = st->asp_p1;
	double *P0 = st->asp_p0;
	int *posh = st->asp_posh;
	int *rw = st->asp_rw;
	int **hci = st->asp_hc_ri;


	int synd;

	int rh = st->rh;
	int nh = st->nh;
	int m  = st->m;
	int r = rh * m;
	int n = nh * m;
	int all_cw_2 = st->asp_all_cw_2;


	for( i = 0; i < n; i++ )
	{
		double x = soft[i] * 0.5;
		double y = maxd( mind( x,  INPUT_LIMIT ), -INPUT_LIMIT );
		double e0 = exp(  y );
		double e1 = exp( -y );
		soft[i] = e1 / (e0 + e1 );
	}

	for( j = 0; j < rh; j++ )
	{
		int cnt = 0;

		for( i = 0; i < nh; i++ )
		{
			int pos_n = i * m;
			int pos_r = j * m;
			int circ = hd[j][i];

			if( circ != -1 )
			{
				rotate( &soft[pos_n], &state[cnt][pos_r], circ, sizeof(soft[0]), m );
				cnt += 1;
			}
		}
		rw[j] = cnt;
	}


	// just to compute syndrome before iterations
	for( i = 0; i < n; i++ )
		soft_out[i] = soft[i];

	// check input codeword
	synd = check_syndrome( syndr, r, hd, rh, nh, m, soft_out, data0, 0.5 );

	if( synd == 0 ) 
	{
		make_output( decision, n, decword, soft_out, 0.5 );
		return 0; 
	}

	steps = 0; // number of iterations
	while( steps < maxsteps )
	{
		//	START ITERATIONS

		// check nodes processing
		for( i = 0; i < rh; i++ ) 
		{
			// decode constituent code of each row
			for( k = 0; k < m; k++ )
				map_bin( &state[0][i*m+k],rw[i], r );
		}


		//symbol nodes
		if( all_cw_2 )
		{
			//overall products
			for( j = 0; j < rh; j++ )
				posh[j] = 0;

			for( i = 0; i < nh; i++ )
			{
				int pos_n  = i * m;
				int j0     = hci[i][0];
				int j1     = hci[i][1];
				int pos_r0 = j0 * m;
				int pos_r1 = j1 * m;
				int circ0  = hd[j0][i];
				int circ1  = hd[j1][i];
				int ph0    = posh[j0];
				int ph1    = posh[j1];

				rotate( &state[ph0][pos_r0], data0, m - circ0, sizeof(state[0][0]), m );
				rotate( &state[ph1][pos_r1], data1, m - circ1, sizeof(state[0][0]), m );

				for( k = 0; k < m; k++ )
				{
					double p1 = soft[pos_n + k]; 
					double q10 = p1;   
					double q11 = p1;
					double p0  = 1.0 - p1;
					double q00 = 1.0 - p1;    
					double q01 = 1.0 - p1;    

					q10 = q10 *      data1[k];
					q00 = q00 * (1 - data1[k]);
					q11 = q11 *      data0[k]; 
					q01 = q01 * (1 - data0[k]);
					p1  = q10 *      data0[k];
					p0  = q00 * (1 - data0[k]);

					// Normalizations
					soft_out[pos_n + k] = p1 / (p0 + p1);

					data0[k] = q10 / (q10 + q00);
					data1[k] = q11 / (q11 + q01);
				}

				rotate( data0, &state[ph0][pos_r0], circ0, sizeof(state[0][0]), m );
				rotate( data1, &state[ph1][pos_r1], circ1, sizeof(state[0][0]), m );

				posh[j0] += 1;
				posh[j1] += 1;
			}
		}
		else   
		{
			//overall products
			for( j = 0; j < rh; j++ )
				posh[j] = 0;

			for( i = 0; i < nh; i++ )
			{
				int pos_n = i * m;

				for( k = 0; k < m; k++ )
				{
					P1[k] = soft[pos_n + k];
					P0[k] = 1 - soft[pos_n + k];
				}

				for( j = 0; j < rh; j++ )
				{
					int circ = hd[j][i];

					if( circ != -1 ) 
					{
						int pos_r = j * m;
						int ph   = posh[j];

						rotate( &state[ph][pos_r], data0, m - circ, sizeof(state[0][0]), m );

						for( k = 0; k < m; k++ )
						{
							P1[k] *= data0[k];
							P0[k] *= 1 - data0[k];
						}

						posh[j] += 1;
					}
				}

				for( k = 0; k < m; k++ )
					soft_out[pos_n + k] = P1[k] / (P0[k] + P1[k]);
			}

			//local data updating
			for( j = 0; j < rh; j++ )
				posh[j] = 0;

			for( i = 0; i < nh; i++ )
			{
				for( j = 0; j < rh; j++ )
				{
					int circ = hd[j][i];

					if( circ != -1 ) 
					{
						int ph    = posh[j];
						int pos_n = i * m;
						int pos_r = j * m;

						rotate( &state[ph][pos_r], data0, m - circ, sizeof(state[0][0]), m );

						for( k = 0; k < m; k++ )
						{
							double so  = soft_out[pos_n + k];
							double sos = data0[k];
							double p1 = so / sos;
							double p0 = (1 - so) / ( 1 - sos );
							double d = p1 / (p1 + p0);

							data0[k] = maxd( mind( d,  SP_DEC_MAX_VAL ), SP_DEC_MIN_VAL );
						}

						rotate( data0, &state[ph][pos_r], circ, sizeof(state[0][0]), m );

						posh[j] += 1;
					}
				}
			}
		}   

		//check syndrome
		synd = check_syndrome( syndr, r, hd, rh, nh, m, soft_out, data0, 0.5 );

		if( synd == 0 ) 
		{
			make_output( decision, n, decword, soft_out, 0.5 );

			return steps+1; 
		}

		steps = steps+1;
	}

	make_output( decision, n, decword, soft_out, 0.5 );

	return -steps;  // errors detected but not corrected
}


int tdmp_sum_prod_gf2_decod_qc_lm( DEC_STATE* st, double soft[], double decword[], int maxsteps, int decision )
{ 
	int i, j, k;
	int steps;
	short **hd = st->hd;
	double **Z = st->tasp_state;
	double *soft_out = st->tasp_soft_out;
	short *syndr    = st->syndr;
	double *data0 = st->tasp_data0;
//	double *data1 = st->tasp_data1;
//	double *P1 = st->tasp_p1;
//	double *P0 = st->tasp_p0;
//	int *posh = st->tasp_posh;
	double T = 0.0001;//1e-5;
	double TT = 0;
	double *y = st->tasp_tmp;

	int synd;

	int rh = st->rh;
	int nh = st->nh;
	int m  = st->m;
	int r = rh * m;
	int n = nh * m;


	for( i = 0; i < n; i++ )
	{
		double x = soft[i] * 0.5;
		double y = maxd( mind( x,  INPUT_LIMIT ), -INPUT_LIMIT );
		double e0 = exp(  y );
		double e1 = exp( -y );
		soft[i] = e1 / (e0 + e1 );
	}
	

	for( j = 0; j < rh; j++ )
	{
		int pos_r  = j * m;
		int cnt;

		cnt = 0;
		for( i = 0; i < nh; i++ )
		{
			int pos_n = i * m;
			int circ = hd[j][i];
			if( circ == 0 )
				circ = m;

			if( circ != -1 )
			{
				for( k = 0; k < m; k++ )
					Z[pos_r + k][cnt] = 0.5;
				cnt++;
			}
		}
	}


	// just to compute syndrome before iterations
	for( i = 0; i < n; i++ )
		soft_out[i] = soft[i];

	synd = check_syndrome( syndr, r, hd, rh, nh, m, soft_out, data0, 0.5 );
	if( synd == 0 )
	{
		for( i = 0; i < n; i++ )
			decword[i] = soft_out[i] > 0.5;

		return 0;
	}

	steps = 0; // number of iterations
	while( steps < maxsteps )
	{
		//	START ITERATIONS
		//check syndrome

		for( i = 0; i < rh; i++ )	//loop over checks
		{
			int cnt;

			for( k = 0; k < m; k++ )
			{
				double *a = Z[i * m + k];

				cnt = 0;
				for( j = 0; j < nh; j++ )
				{
					int circ = hd[i][j];

					if( circ != -1 )
					{
						int idx = j * m + ((k + circ) % m);
						double x = soft_out[idx];          // gamma

						y[cnt] = x * (1.0 - a[cnt]) / (a[cnt] + x- 2.0 * a[cnt] * x);   // rho=gamma-lambda

						cnt++;
					}
				}

				for( j = 0; j < cnt; j++ )
				{
					if( y[j] < TT )	y[j] = TT;
					if( y[j] > 1 - TT ) y[j] = 1 - TT;
					a[j] = y[j];
				}

				map_bin( a, cnt, 1 );

				for( j = 0; j < cnt; j++ )
				{
					if( a[j] < T )  a[j] = T;
					if( a[j] > 1.0 - T) a[j] = 1.0 - T;
				}
				
				cnt = 0;
				for( j = 0; j < nh; j++ )
				{
					int circ = hd[i][j];

					if( circ != -1 )
					{
						int idx = j * m + ((k + circ) % m);

						soft_out[idx] = y[cnt] * a[cnt] / (1.0 - y[cnt] - a[cnt] + 2 * y[cnt] * a[cnt]);  // gamma=rho+lambda

						cnt++;
					}
				}
			}

			synd = check_syndrome( syndr, r, hd, rh, nh, m, soft_out, data0, 0.5 );
		}


		steps = steps+1;

		if( synd == 0 ) 
			break;
	}

	for( i = 0; i < n; i++ )
		decword[i] = soft_out[i] > 0.5;

	if( synd == 1 )
		steps = -steps;  // errors detected but not corrected
	
	return steps;
}

#ifndef IASP_FIXED_POINT


double dp0[1000];
double dp1[1000];



//#define KILL_MNT 6

UFLT16 fixed2uflt16( int xfpp , unsigned short x )
{
	UFLT16 y;
	i16 m = (i16)x << (16 - xfpp);
	i16 p = 0;

	while( m >= 0 )
	{
		m <<= 1;
		p  -= 1;
	}
#ifdef UFLT_MNT_16
	y.m = (UFLT_MNT)m;

#ifdef KILL_MNT
	y.m >>= KILL_MNT;
	y.m <<= KILL_MNT;
#endif

#else
	y.m = (UFLT_MNT)((ui16)m >> 8);
#endif
	y.p = (FLT_POW)p;
	
	return y;
}


FLT16 fixed2flt16( int xfpp , short x )
{
	FLT16 y;
	i16 sign = (ui16)x >> 15;
	i16 m = x << (16 - 1 - xfpp);
	i16 p = 0;

	if( m == 0 )
	{
		y.m = 0;
		y.p = 0;
		return y;
	}

	if( sign )
		m = -m;

	while( !(m & 0x4000)  )
	{
		m <<= 1;
		p  -= 1;
	}

	if( sign )
		m = -m;
#ifdef UFLT_MNT_16
	y.m = (FLT_MNT)m;

#ifdef KILL_MNT
	y.m >>= KILL_MNT;
	y.m <<= KILL_MNT;
#endif

#else
	y.m = (FLT_MNT)((ui16)m >> 8);
#endif
	y.p = (FLT_POW)p;

	return y;
}


UFLT16 mul_uflt16( UFLT16 x, UFLT16 y )
{
	UFLT16 z;
#ifdef UFLT_MNT_16
	ui32 m;
	m = (ui32)x.m * (ui32)y.m;
#else
	ui16 m;
	m = (ui16)x.m * (ui16)y.m;
#endif
  
#ifdef UFLT_MNT_16
	if( (i32)m > 0 )
#else
	if( (i16)m > 0 )
#endif
	{
		z.m = (UFLT_MNT)(m >> (sizeof(UFLT_MNT)*8-1)); 
		z.p = (FLT_POW)(x.p + y.p - 1);
	}
	else
	{
		z.m = (UFLT_MNT)(m >>= (sizeof(UFLT_MNT)*8)); 
		z.p = (FLT_POW)(x.p + y.p);
	}

#ifdef KILL_MNT
	z.m >>= KILL_MNT;
	z.m <<= KILL_MNT;
#endif

	return z;
}

FLT16 mul_flt16( FLT16 x, FLT16 y )
{
	FLT16 z; 
#ifdef FLT_MNT_16
	i32 m = (i32)x.m * (i32)y.m;
	i32 absm = m < 0 ? -m : m;
#else
	i16 m = (i16)x.m * (i16)y.m;
	i16 absm = m < 0 ? -m : m;
#endif

#ifdef FLT_MNT_16
	if( !(absm & 0x20000000) )
#else
	if( !(absm & 0x2000) )
#endif
	{
		z.m = (FLT_MNT)(m >> (sizeof(FLT_MNT)*8-2)); 
		z.p = (FLT_POW)(x.p + y.p - 1);
	}
	else
	{
		z.m = (FLT_MNT)(m >>= (sizeof(FLT_MNT)*8-1)); 
		z.p = (FLT_POW)(x.p + y.p);
	}

#ifdef KILL_MNT
	z.m >>= KILL_MNT;
	z.m <<= KILL_MNT;
#endif

	return z;
}


UFLT16 add_uflt16( UFLT16 x, UFLT16 y )
{
	UFLT16 z;
	ui32 m;
	FLT_POW p;

	short d = x.p - y.p;
	if( d > 0 )
	{
		d = mini(d, (short)(sizeof(UFLT_MNT)*8));
		m = (ui32)x.m + (ui32)(y.m >> d);
		p = x.p;
	}
	else
	{
		d = mini((-d), (short)(sizeof(UFLT_MNT)*8));
		m = (ui32)y.m + (ui32)(x.m >> d);
		p = y.p;
	}

	if( m & (1 << (sizeof(UFLT_MNT)*8)) )
	{
		m >>= 1;
		p += 1;
	}

	z.m = (UFLT_MNT)m;
	z.p = (FLT_POW)p;

#ifdef KILL_MNT
	z.m >>= KILL_MNT;
	z.m <<= KILL_MNT;
#endif


	return z;
}



UFLT16 div_uflt16( UFLT16 x, UFLT16 y )
{
	UFLT16 z;
	int p = x.p - y.p + 1;
	unsigned int m = (unsigned int)x.m << (sizeof(UFLT_MNT)*8 - 1);
	unsigned int m1 = (unsigned int)m / (unsigned int)y.m;
	
#ifdef UFLT_MNT_16
	if( (i16)m1 > 0 )
#else
	if( (i8)m1 > 0 )
#endif
	{
		m1 <<= 1;
		p -= 1;
	}
	z.m = (UFLT_MNT)m1;
	z.p = (FLT_POW)p;
	
#ifdef KILL_MNT
	z.m >>= KILL_MNT;
	z.m <<= KILL_MNT;
#endif

	return z;
}

ui32 uflt2fixed( UFLT16 x, int fpp )
{
#ifdef UFLT_MNT_16
	ui32 m = x.m;
#else
	ui32 m = (ui32)x.m << 8;
#endif
	int shift = mini( 16, abs( x.p ) ); 

	m = x.p > 0 ? m << shift : m >> shift;  
	m >>= 16 - fpp;

	return m;
}

ui32 flt2fixed( FLT16 x, int fpp )
{
#ifdef FLT_MNT_16
	i32 m = x.m;
#else
	i32 m = (i32)x.m << 8;
#endif
	int shift = mini( 16, abs( x.p ) ); 

	m = x.p > 0 ? m << shift : m >> shift;  
	m >>= 16 - 1 - fpp;

	return m;
}
double uflt2double( UFLT16 x )
{
#ifdef UFLT_MNT_16
	ui32 m = x.m;
#else
	ui32 m = (ui32)x.m << 8;
#endif
	double r = (double)m * pow(2.0, x.p-16);

	return r;
}

float uflt2float( UFLT16 x )
{
#ifdef UFLT_MNT_16
	ui32 m = x.m;
#else
	ui32 m = (ui32)x.m << 8;
#endif
	float r = (float)m * (float)pow(2.0, x.p-16);

	return r;
}

float flt2float( FLT16 x )
{
#ifdef FLT_MNT_16
	i32 m = x.m;
#else
	i32 m = (i32)x.m << 8;
#endif
	float r = (float)m * (float)pow(2.0, x.p-(16-1));

	return r;
}

UFLT16 kill_mnt( UFLT16 x, int n )
{
	UFLT16 z = x;
    z.m = div_power2(z.m, n);
	z.m <<= n;
	return z;
}


//#define USE_FLOAT
void imap_bin_flt16( UFLT16 soft[], int rw, int step)
{
	int i;
	ui16 tmp;
	const int ONE_15 = 1 << 15;
#ifdef USE_FLOAT
	float fSF[ROW_WEIGHT_MAX];
	float fSB[ROW_WEIGHT_MAX];
	float  fP[ROW_WEIGHT_MAX];
	short  res[ROW_WEIGHT_MAX];
#endif
	//i16 acc;
	FLT16 acc;
	FLT16 SF[ROW_WEIGHT_MAX];
	FLT16 SB[ROW_WEIGHT_MAX];
	FLT16 P[ROW_WEIGHT_MAX];


	for( i = 0; i < rw; i++ )
	{
		i16 p = ONE_15 - 2 * uflt2fixed( soft[i*step], 15 );
//		ui32 k = uflt2fixed( soft[i*step], 15 );
//		double f = uflt2double(soft[i*step]); 
//		k *= 2;
//		k = ONE_15 - k;
		P[i] = fixed2flt16( 15, p );
	}

#ifdef USE_FLOAT
	for( i = 0; i < rw; i++ )
		fP[i] = 1.0f - 2.0f * uflt2float( soft[i*step] );
	
	fSF[0] = fP[0]; 
	for( i = 1; i < rw-1; i++ )
		fSF[i] = fP[i] * fSF[i-1];

	fSB[rw-1] = fP[rw-1]; 
	for( i = rw-2; i > 0; i-- )
		fSB[i] = fP[i] * fSB[i+1];

#endif

	// Alpha
	acc = SF[0] = P[0]; 
	for( i = 1; i < rw-1; i++ )
	{
//		float f;
//		float p = flt2float( P[i] );
//		float a = flt2float( acc );
//		float r = p * a;
		SF[i] = acc = mul_flt16( P[i], acc );
//		f = flt2float( SF[i] );
//		r=r;
	}

	// Beta
	acc = SB[rw-1] = P[rw-1]; 
	for( i = rw-2; i > 0; i-- )
	{
//		float f;
//		float p = flt2float( ifP[i] );
//		float a = flt2float( ifacc );
//		float r = p * a;
		SB[i] = acc = mul_flt16( P[i], acc );
//		f = flt2float( ifSB[i] );
//		r=r;
	}

	// Sigma
	tmp = div_power2r(ONE_15 - flt2fixed(SB[1], 15), 1);
	soft[0*step] = fixed2uflt16( 15, maxi( tmp, 1 ) );


#ifdef USE_FLOAT
	res[0] = ( ((1.0f)-fSB[1])/2.0f ) * ONE_15;
	soft[0*step] = fixed2uflt16( 15, maxi( res[0], 1 ) );
#endif
	for( i = 1; i < rw-1; i++ ) //loop over symbols in the current check
	{
		FLT16 Z = mul_flt16(SF[i-1], SB[i+1]);

		tmp = div_power2r(ONE_15 - flt2fixed(Z, 15), 1);
		soft[i*step] = fixed2uflt16( 15, maxi( tmp, 1 ) );
#ifdef USE_FLOAT
		{
			float Z = fSF[i-1] * fSB[i+1] ;// / (1 << 15);

			res[i] = (((1.0f)-Z)/2.0f) * ONE_15;
			soft[i*step] = fixed2uflt16( 15, maxi( res[i], 1 ) );
		}
#endif
	}
	tmp = div_power2r(ONE_15 - flt2fixed(SF[rw-2], 15), 1);
	soft[(rw-1)*step] = fixed2uflt16( 15, maxi( tmp, 1 ) );

#ifdef USE_FLOAT
	res[(rw-1)] = (((1.0f) - fSF[rw-2])/2.0f) * ONE_15;
	soft[(rw-1)*step] = fixed2uflt16( 15, maxi( res[(rw-1)], 1 ) );
#endif
}

int isum_prod_gf2_decod_qc_lm( DEC_STATE* st, double soft[], double decword[], int maxsteps, int decision )
{ 
	int i, j, k;
	int steps;
	short **hd = st->hd;
	short *syndr    = st->syndr;
	UFLT16 **state = st->iasp_state;
	UFLT16 *y = st->iasp_y;
	ui16 *soft_out = st->iasp_soft_out;
	ui16 *data = st->iasp_data;
	UFLT16 *data0 = st->iasp_data0;
	UFLT16 *data1 = st->iasp_data1;
	UFLT16 *P1 = st->iasp_p1;
	UFLT16 *P0 = st->iasp_p0;
	int *posh = st->iasp_posh;
	int *rw = st->iasp_rw;
	int **hci = st->iasp_hci;


	int synd;

	int rh = st->rh;
	int nh = st->nh;
	int m  = st->m;
	int r = rh * m;
	int n = nh * m;
	int all_cw_2 = st->iasp_all_cw_2;

	for( i = 0; i < n; i++ )
	{
		double x = soft[i];
		double y = maxd( mind( x,  INPUT_LIMIT ), -INPUT_LIMIT );
		soft[i] = 1.0 / (1.0 + exp( y ));
	}

	for( i = 0; i < n; i++ )
	{
		ui32 x = (int)(soft[i] * ONE_SOFT + 0.5);
		x = mini( x, MAX_SOFT );
		x = maxi( x, 1 );
		y[i] = fixed2uflt16( SOFT_FPP, x );
		soft_out[i] = (ui16)x << (16 - SOFT_FPP);
	}


	for( j = 0; j < rh; j++ )
	{
		int cnt = 0;

		for( i = 0; i < nh; i++ )
		{
			int pos_n = i * m;
			int pos_r = j * m;
			int circ = hd[j][i];

			if( circ != -1 )
			{
				rotate( &y[pos_n], &state[cnt][pos_r], circ, sizeof(y[0]), m );
				cnt += 1;
			}
		}
		rw[j] = cnt;
	}


	steps = 0; // number of iterations
	while( steps < maxsteps )
	{
		//	START ITERATIONS
		//check syndrome
		for( i = 0; i < r; i++ )
			syndr[i] = 0;

		for( j = 0; j < rh; j++ )
		{
			for( i = 0; i < nh; i++ ) 
			{
				int circ = hd[j][i];

				if( circ != - 1 )
				{
					int pos_r = j * m;
					int pos_n = i * m;

					rotate( &soft_out[pos_n], data, circ, sizeof(soft_out[0]), m );

					for( k = 0; k < m; k++ )
						syndr[pos_r + k] ^= data[k] >> 15;
				}
			}
		}

		synd = 0;
		for( i = 0; i < r; i++ )
			synd |= syndr[i];

		if( synd == 0 ) 
		{
			if( decision )
			{
				for( i = 0; i < n; i++ )
					decword[i] = (double)soft_out[i] / (1 << 16);
			}
			else
			{
				for( i = 0; i < n; i++ )
					decword[i] = soft_out[i] >> 15;
			}

			return steps; 
		}

		// check nodes processing
		for( i = 0; i < rh; i++ ) 
		{
			// decode constituent code of each row
			for( k = 0; k < m; k++ )
				imap_bin_flt16( &state[0][i*m+k], rw[i], r );
		}


		//symbol nodes
		if( all_cw_2 )
		{
			//overall products
			for( j = 0; j < rh; j++ )
				posh[j] = 0;
/*
			for( i = 0; i < nh; i++ )
			{
				int pos_n  = i * m;
				int j0     = hci[i][0];
				int j1     = hci[i][1];
				int pos_r0 = j0 * m;
				int pos_r1 = j1 * m;
				int circ0  = hd[j0][i];
				int circ1  = hd[j1][i];
				int ph0    = posh[j0];
				int ph1    = posh[j1];

				rotate( &state[ph0][pos_r0], data0, m - circ0, sizeof(state[0][0]), m );
				rotate( &state[ph1][pos_r1], data1, m - circ1, sizeof(state[0][0]), m );

				for( k = 0; k < m; k++ )
				{
					ui16 ip1 = y[pos_n + k]; 
					ui16 ip0 = (ui16)((ui32)(1 << 16) - ip1); 

					ui16 d1 = data1[k] << (16 - SOFT_FPP);
					ui16 d0 = data0[k] << (16 - SOFT_FPP);
					ui16 t1 = ((ONE_SOFT) << (16-SOFT_FPP)) - d1;
					ui16 t0 = ((ONE_SOFT) << (16-SOFT_FPP)) - d0;

					ui16 q10 = (ui16)div_power2r( (ui32)ip1 * d1, 16 );
					ui16 q11 = (ui16)div_power2r( (ui32)ip1 * d0, 16 ); 
					ui16 q00 = (ui16)div_power2r( (ui32)ip0 * t1, 16 );
					ui16 q01 = (ui16)div_power2r( (ui32)ip0 * t0, 16 );

					ui16 p1  = (ui16)div_power2r((ui32)q10 * d0, 16 );
					ui16 p0  = (ui16)div_power2r((ui32)q00 * t0, 16 );

					p0 = p1 + p0;
					p0 = maxi( p0, 1 );
					// Normalizations
					soft_out[pos_n + k] = ((ui32)p1 << 16) / p0;
					soft_out[pos_n + k] = maxi( soft_out[pos_n + k], 1 << (16-SOFT_FPP) );

					q00 = q00 + q10;
					q01 = q01 + q11; 
					q00 = maxi( q00, 1 );
					q01 = maxi( q01, 1 );
					data0[k] =  ((ui32)q10 << SOFT_FPP) / q00;
					data1[k] =  ((ui32)q11 << SOFT_FPP) / q01;

					data0[k] = maxi( data0[k], 1 );
					data1[k] = maxi( data1[k], 1 );
				}

				rotate( data0, &state[ph0][pos_r0], circ0, sizeof(state[0][0]), m );
				rotate( data1, &state[ph1][pos_r1], circ1, sizeof(state[0][0]), m );

				posh[j0] += 1;
				posh[j1] += 1;
			}
*/
		}
		else   
		{
			//overall products
			for( j = 0; j < rh; j++ )
				posh[j] = 0;

			for( i = 0; i < nh; i++ )
			{
				int pos_n = i * m;

				for( k = 0; k < m; k++ )
				{
#if 0
					{
					double a1 = uflt2double( y[pos_n + k] );
					double a0 = 1.0 - a1;

					if( a0 == 0 )
						printf("a0=0\n");
					dp1[k] = a1;
					dp0[k] = a0;
					}
#else
					{
					P1[k]  = y[pos_n + k];
					ui16 a1 = uflt2fixed( P1[k], 16 ); 
//					ui16 a0 = (ui16)((ui32)(1 << 16) - a1);
					ui16 a0 = ~a1 + 1;
					P0[k]  = fixed2uflt16( 16, a0 );
					}
#endif
				}

				for( j = 0; j < rh; j++ )
				{
					int circ = hd[j][i];

					if( circ != -1 ) 
					{
						int pos_r = j * m;
						int ph   = posh[j];

						rotate( &state[ph][pos_r], data0, m - circ, sizeof(state[0][0]), m );

//						for( k = 0; k < m; k++ )	data0[k] = kill_mnt( data0[k], 6 );

						for( k = 0; k < m; k++ )
						{
#if 0
							{
							double d1 = uflt2double( data0[k] );
							double d0 = 1.0 - d1;
							if( d0 == 0 )
								printf("d0=0\n");
							dp1[k] *= d1;
							dp0[k] *= d0;
							}
#else
							{
							UFLT16 fd1 = data0[k];
							ui16 d1   = uflt2fixed( fd1, 16 );
							ui16 d0   = ~d1 + 1;
							UFLT16 fd0 = fixed2uflt16( 16, d0 );
							

							P1[k] = mul_uflt16( P1[k], fd1 ); 
							P0[k] = mul_uflt16( P0[k], fd0 ); 
							}
#endif
						}

						posh[j] += 1;
					}
				}

				for( k = 0; k < m; k++ )
				{
#if 0
					ui32 s = (int)((dp1[k] / (dp0[k] + dp1[k])) * ONE_SOFT);
#else
					UFLT16 t = add_uflt16( P0[k], P1[k] );
					UFLT16 q = div_uflt16( P1[k], t );

					ui32 s = uflt2fixed( q, SOFT_FPP );
#endif
					s = mini( s, MAX_SOFT );
					s = maxi( s, 1 );

					soft_out[pos_n + k] = s << (16 - SOFT_FPP);
				}
			}

			//local data updating
			for( j = 0; j < rh; j++ )
				posh[j] = 0;

			for( i = 0; i < nh; i++ )
			{
				for( j = 0; j < rh; j++ )
				{
					int circ = hd[j][i];

					if( circ != -1 ) 
					{
						int ph    = posh[j];
						int pos_n = i * m;
						int pos_r = j * m;

						rotate( &state[ph][pos_r], data0, m - circ, sizeof(state[0][0]), m );

//						for( k = 0; k < m; k++ )	data0[k] = kill_mnt( data0[k], 6 );

						for( k = 0; k < m; k++ )
						{
#if 1
							double dso  = (double)soft_out[pos_n + k] / (1 << 16);
							double dsos = uflt2double( data0[k] );
							double p1 = dso / dsos;
							double p0 = (1.0 - dso)/(1.0 - dsos );
							double fd = p1 / (p1 + p0);
							ui32 d = (ui32)(fd * (1 << SOFT_FPP));
							if( dso == 0 || dsos == 0 )
								printf("ds=0\n");

							d = mini( maxi( d, 1 ), IASP_DEC_MAX_VAL );
							data0[k] = fixed2uflt16( SOFT_FPP, d );
#else
							ui16 so  = soft_out[pos_n + k];

							UFLT16 fso  = fixed2uflt16( 16, so );
							UFLT16 fsos = data0[k];
							UFLT16 fp1  = div_uflt16( fso, fsos );

							ui16 iso  = uflt2fixed( fso, 16 );
							ui16 isos = uflt2fixed( fsos, 16 );
							ui16 so1  = ~iso + 1;				//(1<<16) - iso
							ui16 sos1 = ~isos + 1;				//(1<<16) - isos

							UFLT16 fso1   = fixed2uflt16( 16,so1 );
							UFLT16 fsos1  = fixed2uflt16( 16,sos1 );
							UFLT16 fp0    = div_uflt16( fso1, fsos1 );

							UFLT16 t = add_uflt16( fp0, fp1 );
							UFLT16 q = div_uflt16( fp1, t );
							ui32 d = uflt2fixed( q, SOFT_FPP );

							d = mini( maxi( d, 1 ), IASP_DEC_MAX_VAL );
							data0[k] = fixed2uflt16( SOFT_FPP, d );
#endif
						}

//						for( k = 0; k < m; k++ )	data0[k] = kill_mnt( data0[k], 6 );

						rotate( data0, &state[ph][pos_r], circ, sizeof(state[0][0]), m );

						posh[j] += 1;
					}
				}
			}
		}   
		steps = steps+1;
	}

	if( decision )
	{
		for( i = 0; i < n; i++ )
			decword[i] = (double)soft_out[i] / (1 << 16);
	}
	else
	{
		for( i = 0; i < n; i++ )
			decword[i] = soft_out[i] >> 15;
	}


	return -steps;  // errors detected but not corrected
}
#else //IASP_FIXED_POINT

int icheck_syndrome( short *syndr, int r, short **hd, int rh, int nh, int m, ui16 *soft, ui16 *buf )
{
	int synd;
	int i, j, k;

	for( i = 0; i < r; i++ )
		syndr[i] = 0;

	for( j = 0; j < rh; j++ )
	{
		for( i = 0; i < nh; i++ ) 
		{
			int circ = hd[j][i];

			if( circ != - 1 )
			{
				int pos_r = j * m;
				int pos_n = i * m;

				rotate( &soft[pos_n], buf, circ, sizeof(soft[0]), m );

				for( k = 0; k < m; k++ )
					syndr[pos_r + k] ^= buf[k] >> 15;
			}
		}
	}

	synd = 0;
	for( i = 0; i < r; i++ )
		synd |= syndr[i];

	return synd;
}

void imake_output( int mode, int n, double *outword, ui16 *soft )
{
	int i;

	if( mode )
	{
		for( i = 0; i < n; i++ )
			outword[i] = (double)soft[i] / (1 << 16);
	}
	else
	{
		for( i = 0; i < n; i++ )
			outword[i] = soft[i] >> 15;
	}
}

int isum_prod_gf2_decod_qc_lm( DEC_STATE* st, double soft[], double decword[], int maxsteps, int decision )
{ 
	int i, j, k;
	int steps;
	short **hd = st->hd;
	short *syndr    = st->syndr;
	ui16 **state = st->iasp_state;
	ui16 *y = st->iasp_y;
	ui16 *soft_out = st->iasp_soft_out;
	ui16 *data0 = st->iasp_data0;
	ui16 *data1 = st->iasp_data1;
	ui32 *P1 = st->iasp_p1;
	ui32 *P0 = st->iasp_p0;
	int *posh = st->iasp_posh;
	int *rw = st->iasp_rw;
	int **hci = st->iasp_hc_ri;


	int synd;

	int rh = st->rh;
	int nh = st->nh;
	int m  = st->m;
	int r = rh * m;
	int n = nh * m;
	int all_cw_2 = st->iasp_all_cw_2;

	for( i = 0; i < n; i++ )
	{
		double x = soft[i];
		double y = maxd( mind( x,  INPUT_LIMIT ), -INPUT_LIMIT );
		soft[i] = 1.0 / (1.0 + exp( y ));
	}

	for( i = 0; i < n; i++ )
	{
		int x = (int)(soft[i] * ONE_SOFT + 0.5);
		x = mini( x, MAX_SOFT );
		y[i] = maxi( x, 1 );
	}


	for( i = 0; i < n; i++ )	soft_out[i] = y[i];

//	if( !all_cw_2  )
		for( i = 0; i < n; i++ )	y[i] = y[i] << (16 - SOFT_FPP);

	for( j = 0; j < rh; j++ )
	{
		int cnt = 0;

		for( i = 0; i < nh; i++ )
		{
			int pos_n = i * m;
			int pos_r = j * m;
			int circ = hd[j][i];

			if( circ != -1 )
			{
				rotate( &soft_out[pos_n], &state[cnt][pos_r], circ, sizeof(soft_out[0]), m );
				cnt += 1;
			}
		}
		rw[j] = cnt;
	}

//	if( !all_cw_2  )
		for( i = 0; i < n; i++ )	soft_out[i] = soft_out[i] << (16 - SOFT_FPP);

	synd = icheck_syndrome( syndr, r, hd, rh, nh, m, soft_out, data0 );

	if( synd == 0 ) 
	{
		imake_output( decision, n, decword, soft_out );

		return 0; 
	}

	steps = 0; // number of iterations
	while( steps < maxsteps )
	{
		//	START ITERATIONS

		// check nodes processing
		for( i = 0; i < rh; i++ ) 
		{
			// decode constituent code of each row
			for( k = 0; k < m; k++ )
				imap_bin( &state[0][i*m+k],rw[i], r );
		}


		//symbol nodes
		if( all_cw_2 )
		{
			//overall products
			for( j = 0; j < rh; j++ )
				posh[j] = 0;

			for( i = 0; i < nh; i++ )
			{
				int pos_n  = i * m;
				int j0     = hci[i][0];
				int j1     = hci[i][1];
				int pos_r0 = j0 * m;
				int pos_r1 = j1 * m;
				int circ0  = hd[j0][i];
				int circ1  = hd[j1][i];
				int ph0    = posh[j0];
				int ph1    = posh[j1];

				rotate( &state[ph0][pos_r0], data0, m - circ0, sizeof(state[0][0]), m );
				rotate( &state[ph1][pos_r1], data1, m - circ1, sizeof(state[0][0]), m );

				for( k = 0; k < m; k++ )
				{
					ui16 ip1 = y[pos_n + k]; 
					ui16 ip0 = ((ONE_SOFT) << (16-SOFT_FPP)) - ip1; 

					ui16 d1 = data1[k] << (16 - SOFT_FPP);
					ui16 d0 = data0[k] << (16 - SOFT_FPP);
					ui16 t1 = ((ONE_SOFT) << (16-SOFT_FPP)) - d1;
					ui16 t0 = ((ONE_SOFT) << (16-SOFT_FPP)) - d0;

					ui16 q10 = (ui16)div_power2r( (ui32)ip1 * d1, 16 );
					ui16 q11 = (ui16)div_power2r( (ui32)ip1 * d0, 16 ); 
					ui16 q00 = (ui16)div_power2r( (ui32)ip0 * t1, 16 );
					ui16 q01 = (ui16)div_power2r( (ui32)ip0 * t0, 16 );
					
					ui16 p1  = (ui16)div_power2r((ui32)q10 * d0, 16 );
					ui16 p0  = (ui16)div_power2r((ui32)q00 * t0, 16 );

					p0 = p1 + p0;
					p0 = maxi( p0, 1 );
					// Normalizations
					soft_out[pos_n + k] = ((ui32)p1 << 16) / p0;
					soft_out[pos_n + k] = maxi( soft_out[pos_n + k], 1 << (16-SOFT_FPP) );

					q00 = q00 + q10;
					q01 = q01 + q11; 
					q00 = maxi( q00, 1 );
					q01 = maxi( q01, 1 );
					data0[k] =  ((ui32)q10 << SOFT_FPP) / q00;
					data1[k] =  ((ui32)q11 << SOFT_FPP) / q01;
					
					data0[k] = maxi( data0[k], 1 );
					data1[k] = maxi( data1[k], 1 );
				}

				rotate( data0, &state[ph0][pos_r0], circ0, sizeof(state[0][0]), m );
				rotate( data1, &state[ph1][pos_r1], circ1, sizeof(state[0][0]), m );

				posh[j0] += 1;
				posh[j1] += 1;
			}
		}
		else   
		{
			//overall products
			for( j = 0; j < rh; j++ )
				posh[j] = 0;

			for( i = 0; i < nh; i++ )
			{
				int pos_n = i * m;

				for( k = 0; k < m; k++ )
				{
					P1[k] = (ui32)y[pos_n + k] << (P_FPP - 16);
					P0[k] = (ui32)((ONE_SOFT << (16-SOFT_FPP)) - y[pos_n + k]) << (P_FPP - 16);
				}

				for( j = 0; j < rh; j++ )
				{
					int circ = hd[j][i];

					if( circ != -1 ) 
					{
						int pos_r = j * m;
						int ph   = posh[j];

						rotate( &state[ph][pos_r], data0, m - circ, sizeof(state[0][0]), m );

//						for( k = 0; k < m; k++ ) data0[k] >>= (16-SOFT_FPP);

						for( k = 0; k < m; k++ )
						{
							ui16 d1 = data0[k] << (16 - SOFT_FPP);
							ui16 d0 = (MAX_SOFT - data0[k]) << (16 - SOFT_FPP);		//MAX_SOFT not ONE_SOFT!!!
#if 1
							ui64 pp1 = (ui64)P1[k] * d1;
							ui64 pp0 = (ui64)P0[k] * d0;

							P1[k] = (ui32)div_power2( pp1, 16 );
							P0[k] = (ui32)div_power2( pp0, 16 );
#else
							// it is bad
							ui64 pp1 = (P1[k] >> 10) * (d1 >> 6);
							ui64 pp0 = (P0[k] >> 10) * (d0 >> 6);

							P1[k] = (ui32)pp1;
							P0[k] = (ui32)pp0;
#endif
						}

						posh[j] += 1;
					}
				}

				for( k = 0; k < m; k++ )
				{
//					soft_out[pos_n + k] = (int)(((double)iP1[k] / (iP0[k] + iP1[k])) * ONE_SOFT);
					int s;
					ui32 x = P1[k] >> 1;
					ui32 y = (P0[k] >> 1) + x;
					int flg = y > (ONE_SOFT << 4);	// why 4? 

					if( flg )
						y = div_power2( y, SOFT_FPP );
					else
						x = x << SOFT_FPP;

					y = maxi( y, 1 );
					s = x / y;

					s = mini( s, MAX_SOFT );
					soft_out[pos_n + k] = maxi( s, 1 );
					
					soft_out[pos_n + k] <<= (16 - SOFT_FPP);
				}
			}

			//local data updating
			for( j = 0; j < rh; j++ )
				posh[j] = 0;

			for( i = 0; i < nh; i++ )
			{
				for( j = 0; j < rh; j++ )
				{
					int circ = hd[j][i];

					if( circ != -1 ) 
					{
						int ph    = posh[j];
						int pos_n = i * m;
						int pos_r = j * m;

						rotate( &state[ph][pos_r], data0, m - circ, sizeof(state[0][0]), m );
						
//						for( k = 0; k < m; k++ ) data0[k] >>= (16-SOFT_FPP);

						for( k = 0; k < m; k++ )
						{
#if 01
							int so  = soft_out[pos_n + k] << (SOFT_FPP - (16-SOFT_FPP));
							int sos = maxi(data0[k], 1);
							int p1  = so / sos;
							int t   = maxi( ONE_SOFT - sos, 1 );
							int p0  = (ONE_SOFT*ONE_SOFT - so) / t;
							int y   = div_power2r(p1 + p0, SOFT_FPP/2);
							int y1  = maxi( y, 1 );
							int d   = (p1 << (SOFT_FPP - SOFT_FPP/2)) / y1;
#else
							int so  = soft_out[pos_n + k]
							int sos = data0[k];
							double p1 = (double)so / (double)sos;
							double p0 = (double)(ONE_SOFT - so) / (double)( ONE_SOFT - sos );
							int d = (int)(p1 / (p1 + p0) * ONE_SOFT + 0.5);
#endif
							data0[k] = mini( maxi( d, 1 ), IASP_DEC_MAX_VAL );
						}

//						for( k = 0; k < m; k++ ) data0[k] <<= (16-SOFT_FPP);

						rotate( data0, &state[ph][pos_r], circ, sizeof(state[0][0]), m );

						posh[j] += 1;
					}
				}
			}
		}   
		
		//check syndrome
		synd = icheck_syndrome( syndr, r, hd, rh, nh, m, soft_out, data0 );

		if( synd == 0 ) 
		{
			imake_output( decision, n, decword, soft_out );

			return steps+1; 
		}

		steps = steps+1;
	}

	imake_output( decision, n, decword, soft_out );

	return -steps;  // errors detected but not corrected
}
#endif //IASP_FIXED_POINT



/*
int min_sum_decod_qc_lm( DEC_STATE* st, double y[], short hard[], int maxsteps, double alpha )
{
	int   i, j, k, n;
	int iter;
	int parity;
	short **matr       = st->hd;
	short *synd        = st->syndr; 
	MS_DATA *soft      = st->ms_soft;
	short *BnNS        = st->ms_BnNS;  
	MS_DATA *buffer    = st->ms_buffer;
	MS_DATA *rbuffer   = st->ms_rbuffer;
	MS_DATA *rsoft     = st->ms_rsoft;
	MS_DEC_STATE *dcs  = st->ms_dcs;		
	MS_DEC_STATE *tmps = st->ms_tmps;		


	int m  = st->m;
	int rh = st->rh;
	int nh = st->nh;
	
	int c_ldpc = m;
	int r_ldpc = m * rh;
	int n_ldpc = m * nh;


	for( i = 0; i < r_ldpc; i++ )
	{
		dcs[i].min1 = 0;
		dcs[i].min2 = 0;
		dcs[i].pos  = 0;
		dcs[i].sign = 0;
	}

	memset( BnNS, 0, rh*n_ldpc*sizeof(BnNS[0]) );


	for( iter = 0; iter < maxsteps; iter++ )
	{

		int memOffset;


		// INIT_STAGE
		memset( synd, 0, r_ldpc*sizeof(synd[0]) );
		memset( soft, 0, n_ldpc*sizeof(soft[0]) );
		memOffset = 0;

		// STATE 1: compute sum
		for( j = 0; j < rh; j++ ) 
		{
			int stateOffset = j*c_ldpc;

			for( k = 0; k < nh; k++ )
			{
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					for( n = 0; n < c_ldpc; n++ )
					{
						MS_DATA tmp = dcs[stateOffset+n].pos == k ? dcs[stateOffset+n].min2 : dcs[stateOffset+n].min1;

						buffer[n] = (BnNS[memOffset+n] ^ dcs[stateOffset+n].sign) ? -tmp : tmp;
					}

					rotate( buffer, rbuffer, c_ldpc - circ, sizeof(buffer[0]), c_ldpc ); 

					for( n = 0; n < c_ldpc; n++ )
					{
						MS_DATA t = soft[k * c_ldpc + n] + rbuffer[n];
						soft[k * c_ldpc + n] = t;//limit_val( t );
					}

				}

				memOffset += c_ldpc;
			}

		} 

		// STATE 2: 
		for( k = 0; k < n_ldpc; k++ )
		{
			soft[k] = y[k] + soft[k] * alpha; //alpha-normalization
			hard[k]   = soft[k] < 0;
		}



		// STATE 3: update statistic
		memOffset = 0;
		for( j = 0; j < rh; j++ )
		{
			int stateOffset = j*c_ldpc;

			for( k = 0; k < c_ldpc; k++ )
			{
				tmps[k].min1 = INT_MAX;
				tmps[k].min2 = INT_MAX;
				tmps[k].pos  = 0;
				tmps[k].sign = 0;
			}

			for( k = 0; k < nh; k++ )
			{
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					rotate( &soft[k*c_ldpc], rsoft, circ, sizeof(soft[0]), c_ldpc ); 

					for( n = 0; n < c_ldpc; n++ )
						synd[stateOffset+n] ^= rsoft[n] < 0;

					for( n = 0; n < c_ldpc; n++ )
						buffer[n] = dcs[stateOffset+n].pos == k ? dcs[stateOffset+n].min2 : dcs[stateOffset+n].min1;

					for( n = 0; n < c_ldpc; n++ )
					{
						short sign = BnNS[memOffset+n] ^ dcs[stateOffset+n].sign;
						MS_DATA val  = buffer[n] * alpha;
						MS_DATA t    = (BnNS[memOffset+n] ^ dcs[stateOffset+n].sign) ? -val : val;

						t  = rsoft[n] - t;

						sign = t < 0;
						BnNS[memOffset+n] = sign;
						tmps[n].sign ^= sign;

//						val = abs( t );	// incorrect for double t
						val = t < 0.0 ? -t : t;
//						val = (val > MAX_VAL) ? MAX_VAL : val;


						if( val < tmps[n].min1 )      
						{
							tmps[n].pos = k;
							tmps[n].min2 = tmps[n].min1; 
							tmps[n].min1 = val;
						}
						else
						{
							if( val < tmps[n].min2 )
								tmps[n].min2 = val;
						}

					}  
				}


				memOffset += c_ldpc;
			}

			for( k = 0; k < c_ldpc; k++ )
				dcs[stateOffset+k] = tmps[k];

		}

		for( parity = 0, i = 0; i < r_ldpc; i++ )
			parity |= synd[i];


		if( parity == 0 )
			break; 

	} // stage

	return parity ? -iter : iter+1; 
}
*/

#if 1
#define BPS	15 // bits per soft
#define BPV 15 // bits per inner variables
#define MAX_VAL ((1L << BPV) - 1)
#else
#define BPS	3 // bits per soft
#define BPV 5 // bits per inner variables
#define MAX_VAL ((1 << BPV) - 1)
#endif

#define limit_val( x, max_val )	(x) > max_val ? max_val : ((x) < -max_val ? -max_val : (x))




#ifndef MS_MUL_CORRECTION
int min_sum_decod_qc_lm( DEC_STATE* st, double y[], double decword[], int maxsteps, int decision, double alpha )
{
	int   i, j, k, n;
	int iter;
	int parity;
	short **matr       = st->hd;
	short *synd        = st->syndr; 
	MS_DATA *soft      = st->ms_soft;
	short *BnNS        = st->ms_BnNS;  
	MS_DATA *buffer    = st->ms_buffer;
	MS_DATA *rbuffer   = st->ms_rbuffer;
	MS_DATA *rsoft     = st->ms_rsoft;
	MS_DEC_STATE *dcs  = st->ms_dcs;		
	MS_DEC_STATE *tmps = st->ms_tmps;		


	int m  = st->m;
	int rh = st->rh;
	int nh = st->nh;
	
	int c_ldpc = m;
	int r_ldpc = m * rh;
	int n_ldpc = m * nh;


	for( i = 0; i < r_ldpc; i++ )
	{
		dcs[i].min1 = 0;
		dcs[i].min2 = 0;
		dcs[i].pos  = 0;
		dcs[i].sign = 0;
	}
/*	
	for( i = 0; i < n_ldpc; i++ )
	{
		double val = abs(y[i]) * 4;
		int sign = y[i] < 0;
		val = (int)(val + 0.5);
		y[i] = sign ? -val : val;
	}
*/

    {
        double coef;
        double en = 0.0;
  
        for( i = 0; i < n_ldpc; i++ )
            en += y[i] * y[i];
        
        coef = sqrt( n_ldpc / en );
        for( i = 0; i < n_ldpc; i++ )
            y[i] = y[i] * coef;
    }   
    
    memset( BnNS, 0, rh*n_ldpc*sizeof(BnNS[0]) );




	// check input codeword
	for( j = 0; j < rh; j++ )
	{
		for( i = 0; i < nh; i++ ) 
		{
			int pos_r = j * m;
			int pos_n = i * m;
			int circ = matr[j][i];

			if( circ != - 1 )
			{
				rotate( &y[pos_n], rsoft, circ, sizeof(y[0]), m );

				for( k = 0; k < m; k++ )
					synd[pos_r + k] ^= rsoft[k] < 0;
			}
		}
	}

	parity = 0;
	for( i = 0; i < r_ldpc; i++ )
		parity |= synd[i];


	for( iter = 0; iter < maxsteps; iter++ )
	{

		int memOffset;


		// INIT_STAGE
		memset( synd, 0, r_ldpc*sizeof(synd[0]) );
		memset( soft, 0, n_ldpc*sizeof(soft[0]) );
		memOffset = 0;

		// STATE 1: compute sum
		for( j = 0; j < rh; j++ ) 
		{
			int stateOffset = j*c_ldpc;

			for( k = 0; k < nh; k++ )
			{
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					for( n = 0; n < c_ldpc; n++ )
					{
						MS_DATA tmp = dcs[stateOffset+n].pos == k ? dcs[stateOffset+n].min2 : dcs[stateOffset+n].min1;

                        tmp -= alpha;
                        if( tmp < 0 )
                            tmp = 0;
                        
						buffer[n] = (BnNS[memOffset+n] ^ dcs[stateOffset+n].sign) ? -tmp : tmp;
					}

					rotate( buffer, rbuffer, c_ldpc - circ, sizeof(buffer[0]), c_ldpc ); 

					for( n = 0; n < c_ldpc; n++ )
					{
						MS_DATA t = soft[k * c_ldpc + n] + rbuffer[n];
						soft[k * c_ldpc + n] = t;//limit_val( t, MAX_VAL );
					}

				}

				memOffset += c_ldpc;
			}

		} 

		// STATE 2: 
		if( decision )
		{
			for( k = 0; k < n_ldpc; k++ )
			{
				soft[k] = y[k] + soft[k];
				decword[k] = soft[k];
			}
		}
		else
		{
			for( k = 0; k < n_ldpc; k++ )
			{
				soft[k] = y[k] + soft[k];
				decword[k] = soft[k] < 0;
			}
		}


		// STATE 3: update statistic
		memOffset = 0;
		for( j = 0; j < rh; j++ )
		{
			int stateOffset = j*c_ldpc;

			for( k = 0; k < c_ldpc; k++ )
			{
				tmps[k].min1 = MAX_VAL;
				tmps[k].min2 = MAX_VAL;
				tmps[k].pos  = 0;
				tmps[k].sign = 0;
			}

			for( k = 0; k < nh; k++ )
			{
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					rotate( &soft[k*c_ldpc], rsoft, circ, sizeof(soft[0]), c_ldpc ); 

					for( n = 0; n < c_ldpc; n++ )
						synd[stateOffset+n] ^= rsoft[n] < 0;

					for( n = 0; n < c_ldpc; n++ )
                    {
                        MS_DATA tmp = dcs[stateOffset+n].pos == k ? dcs[stateOffset+n].min2 : dcs[stateOffset+n].min1;
                        
                        tmp -= alpha;
                        if( tmp < 0 )
                            tmp = 0;
                        
						buffer[n] = tmp;
                    }

					for( n = 0; n < c_ldpc; n++ )
					{
						short sign = BnNS[memOffset+n] ^ dcs[stateOffset+n].sign;
						MS_DATA val  = buffer[n];// * alpha;
						MS_DATA t    = (BnNS[memOffset+n] ^ dcs[stateOffset+n].sign) ? -val : val;

						t  = rsoft[n] - t;

						sign = t < 0;
						BnNS[memOffset+n] = sign;
						tmps[n].sign ^= sign;

//						val = abs( t );	// incorrect for double t
						val = t < 0.0 ? -t : t;
						val = (val > MAX_VAL) ? MAX_VAL : val;


						if( val < tmps[n].min1 )      
						{
							tmps[n].pos = k;
							tmps[n].min2 = tmps[n].min1; 
							tmps[n].min1 = val;
						}
						else
						{
							if( val < tmps[n].min2 )
								tmps[n].min2 = val;
						}

					}  
				}


				memOffset += c_ldpc;
			}

			for( k = 0; k < c_ldpc; k++ )
				dcs[stateOffset+k] = tmps[k];

		}

		for( parity = 0, i = 0; i < r_ldpc; i++ )
			parity |= synd[i];


		if( parity == 0 )
			break; 

	} // stage

	return parity ? -iter : iter+1; 
}

#else
int min_sum_decod_qc_lm( DEC_STATE* st, double y[], double decword[], int maxsteps, int decision, double alpha )
{
	int   i, j, k, n;
	int iter;
	int parity;
	short **matr       = st->hd;
	short *synd        = st->syndr; 
	MS_DATA *soft      = st->ms_soft;
	short *BnNS        = st->ms_BnNS;  
	MS_DATA *buffer    = st->ms_buffer;
	MS_DATA *rbuffer   = st->ms_rbuffer;
	MS_DATA *rsoft     = st->ms_rsoft;
	MS_DEC_STATE *dcs  = st->ms_dcs;		
	MS_DEC_STATE *tmps = st->ms_tmps;		


	int m  = st->m;
	int rh = st->rh;
	int nh = st->nh;
	
	int c_ldpc = m;
	int r_ldpc = m * rh;
	int n_ldpc = m * nh;


	for( i = 0; i < r_ldpc; i++ )
	{
		dcs[i].min1 = 0;
		dcs[i].min2 = 0;
		dcs[i].pos  = 0;
		dcs[i].sign = 0;
	}
/*	
	for( i = 0; i < n_ldpc; i++ )
	{
		double val = abs(y[i]) * 4;
		int sign = y[i] < 0;
		val = (int)(val + 0.5);
		y[i] = sign ? -val : val;
	}
*/

	memset( BnNS, 0, rh*n_ldpc*sizeof(BnNS[0]) );




	// check input codeword
	for( j = 0; j < rh; j++ )
	{
		for( i = 0; i < nh; i++ ) 
		{
			int pos_r = j * m;
			int pos_n = i * m;
			int circ = matr[j][i];

			if( circ != - 1 )
			{
				rotate( &y[pos_n], rsoft, circ, sizeof(y[0]), m );

				for( k = 0; k < m; k++ )
					synd[pos_r + k] ^= rsoft[k] < 0;
			}
		}
	}

	parity = 0;
	for( i = 0; i < r_ldpc; i++ )
		parity |= synd[i];


	for( iter = 0; iter < maxsteps; iter++ )
	{

		int memOffset;


		// INIT_STAGE
		memset( synd, 0, r_ldpc*sizeof(synd[0]) );
		memset( soft, 0, n_ldpc*sizeof(soft[0]) );
		memOffset = 0;

		// STATE 1: compute sum
		for( j = 0; j < rh; j++ ) 
		{
			int stateOffset = j*c_ldpc;

			for( k = 0; k < nh; k++ )
			{
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					for( n = 0; n < c_ldpc; n++ )
					{
						MS_DATA tmp = dcs[stateOffset+n].pos == k ? dcs[stateOffset+n].min2 : dcs[stateOffset+n].min1;

						buffer[n] = (BnNS[memOffset+n] ^ dcs[stateOffset+n].sign) ? -tmp : tmp;
					}

					rotate( buffer, rbuffer, c_ldpc - circ, sizeof(buffer[0]), c_ldpc ); 

					for( n = 0; n < c_ldpc; n++ )
					{
						MS_DATA t = soft[k * c_ldpc + n] + rbuffer[n];
						soft[k * c_ldpc + n] = t;//limit_val( t, MAX_VAL );
					}

				}

				memOffset += c_ldpc;
			}

		} 

		// STATE 2: 
		if( decision )
		{
			for( k = 0; k < n_ldpc; k++ )
			{
				soft[k] = y[k] + soft[k] * alpha; //alpha-normalization
				decword[k] = soft[k];
			}
		}
		else
		{
			for( k = 0; k < n_ldpc; k++ )
			{
				soft[k] = y[k] + soft[k] * alpha; //alpha-normalization
				decword[k] = soft[k] < 0;
			}
		}


		// STATE 3: update statistic
		memOffset = 0;
		for( j = 0; j < rh; j++ )
		{
			int stateOffset = j*c_ldpc;

			for( k = 0; k < c_ldpc; k++ )
			{
				tmps[k].min1 = MAX_VAL;
				tmps[k].min2 = MAX_VAL;
				tmps[k].pos  = 0;
				tmps[k].sign = 0;
			}

			for( k = 0; k < nh; k++ )
			{
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					rotate( &soft[k*c_ldpc], rsoft, circ, sizeof(soft[0]), c_ldpc ); 

					for( n = 0; n < c_ldpc; n++ )
						synd[stateOffset+n] ^= rsoft[n] < 0;

					for( n = 0; n < c_ldpc; n++ )
						buffer[n] = dcs[stateOffset+n].pos == k ? dcs[stateOffset+n].min2 : dcs[stateOffset+n].min1;

					for( n = 0; n < c_ldpc; n++ )
					{
						short sign = BnNS[memOffset+n] ^ dcs[stateOffset+n].sign;
						MS_DATA val  = buffer[n] * alpha;
						MS_DATA t    = (BnNS[memOffset+n] ^ dcs[stateOffset+n].sign) ? -val : val;

						t  = rsoft[n] - t;

						sign = t < 0;
						BnNS[memOffset+n] = sign;
						tmps[n].sign ^= sign;

//						val = abs( t );	// incorrect for double t
						val = t < 0.0 ? -t : t;
						val = (val > MAX_VAL) ? MAX_VAL : val;


						if( val < tmps[n].min1 )      
						{
							tmps[n].pos = k;
							tmps[n].min2 = tmps[n].min1; 
							tmps[n].min1 = val;
						}
						else
						{
							if( val < tmps[n].min2 )
								tmps[n].min2 = val;
						}

					}  
				}


				memOffset += c_ldpc;
			}

			for( k = 0; k < c_ldpc; k++ )
				dcs[stateOffset+k] = tmps[k];

		}

		for( parity = 0, i = 0; i < r_ldpc; i++ )
			parity |= synd[i];


		if( parity == 0 )
			break; 

	} // stage

	return parity ? -iter : iter+1; 
}
#endif


int imin_sum_decod_qc_lm( DEC_STATE* st, double y[], double decword[], int maxsteps, int decision, double alpha, double thr, int qbits, int dbits )
{
	int   i, j, k, n;
	int iter;
	int parity;
	short **matr        = st->hd;
	short *synd         = st->syndr; 
	IMS_DATA *soft      = st->ims_soft;
	short *BnNS         = st->ims_BnNS;  
	IMS_DATA *buffer    = st->ims_buffer;
	IMS_DATA *rbuffer   = st->ims_rbuffer;
	IMS_DATA *rsoft     = st->ims_rsoft;
	IMS_DEC_STATE *dcs  = st->ims_dcs;		
	IMS_DEC_STATE *tmps = st->ims_tmps;		
	IMS_DATA *iy        = st->ims_y;
	IMS_DATA max_data   = (IMS_DATA)((1L << (dbits-1)) - 1);
	IMS_DATA max_quant  = (IMS_DATA)((1L << (qbits-1)) - 1);


	int m  = st->m;
	int rh = st->rh;
	int nh = st->nh;
	
	int c_ldpc = m;
	int r_ldpc = m * rh;
	int n_ldpc = m * nh;

#ifdef MS_MUL_CORRECTION
	int ialpha = (int)(alpha * (1L << MS_ALPHA_FPP));
#else
	IMS_DATA delta = (IMS_DATA)(max_data * alpha);
#endif

	for( i = 0; i < r_ldpc; i++ )
	{
		dcs[i].min1 = 0;
		dcs[i].min2 = 0;
		dcs[i].pos  = 0;
		dcs[i].sign = 0;
	}


	{
		double coef;
		double en = 0;
	
		for( i = 0; i < n_ldpc; i++ )
			en += y[i] * y[i];

		coef = sqrt(n_ldpc / en);

		for( i = 0; i < n_ldpc; i++ )
		{
			int ival;
			double val = y[i];
			int   sign = 0;

			if( val < 0 )
			{
				val = -val;
				sign = 1;
			}

			val *= coef;
			if( val > thr ) 
				val = thr;
            ival  = (short)floor( val * max_quant / thr + 0.5 );

			iy[i] = sign ? -ival : ival;
		}
	}

	memset( BnNS, 0, rh*n_ldpc*sizeof(BnNS[0]) );
/*
	// check input codeword
	for( j = 0; j < rh; j++ )
	{
		for( i = 0; i < nh; i++ ) 
		{
			int pos_r = j * m;
			int pos_n = i * m;
			int circ = matr[j][i];

			if( circ != - 1 )
			{
				rotate( &y[pos_n], rsoft, circ, sizeof(y[0]), m );

				for( k = 0; k < m; k++ )
					synd[pos_r + k] ^= rsoft[k] < 0;
			}
		}
	}

	parity = 0;
	for( i = 0; i < r_ldpc; i++ )
		parity |= synd[i];
*/

	for( iter = 0; iter < maxsteps; iter++ )
	{

		int memOffset;


		// INIT_STAGE
		memset( synd, 0, r_ldpc*sizeof(synd[0]) );
		memset( soft, 0, n_ldpc*sizeof(soft[0]) );
		memOffset = 0;

		// STATE 1: compute sum
		for( j = 0; j < rh; j++ ) 
		{
			int stateOffset = j*c_ldpc;

			for( k = 0; k < nh; k++ )
			{
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					for( n = 0; n < c_ldpc; n++ )
					{
						IMS_DATA tmp = dcs[stateOffset+n].pos == k ? dcs[stateOffset+n].min2 : dcs[stateOffset+n].min1;
#ifdef MS_MUL_CORRECTION
						tmp = (IMS_DATA)((tmp * ialpha) >> MS_ALPHA_FPP);
#else
						tmp -= delta;
						if( tmp < 0 ) tmp = 0;
#endif

						buffer[n] = (BnNS[memOffset+n] ^ dcs[stateOffset+n].sign) ? -tmp : tmp;
					}

					rotate( buffer, rbuffer, c_ldpc - circ, sizeof(buffer[0]), c_ldpc ); 

					for( n = 0; n < c_ldpc; n++ )
					{
						IMS_DATA t = soft[k * c_ldpc + n] + rbuffer[n];
						soft[k * c_ldpc + n] = limit_val( t, max_data );
					}

				}

				memOffset += c_ldpc;
			}

		} 

		// STATE 2: 
		if( decision )
		{
			for( k = 0; k < n_ldpc; k++ )
			{
#if 0
				soft[k] = iy[k] + (IMS_DATA)((soft[k] * ialpha) >> MS_ALPHA_FPP); //alpha-normalization
#else
				soft[k] = iy[k] + soft[k]; 
#endif
				soft[k] = limit_val( soft[k], max_data );
				decword[k] = soft[k];
			}
		}
		else
		{
			for( k = 0; k < n_ldpc; k++ )
			{
#if 0
				soft[k] = iy[k] + (IMS_DATA)((soft[k] * ialpha) >> MS_ALPHA_FPP); //alpha-normalization
#else
				soft[k] = iy[k] + soft[k]; 
#endif
				soft[k] = limit_val( soft[k], max_data );
				decword[k]   = soft[k] < 0;
			}
		}



		// STATE 3: update statistic
		memOffset = 0;
		for( j = 0; j < rh; j++ )
		{
			int stateOffset = j*c_ldpc;

			for( k = 0; k < c_ldpc; k++ )
			{
				tmps[k].min1 = max_data;
				tmps[k].min2 = max_data;
				tmps[k].pos  = 0;
				tmps[k].sign = 0;
			}

			for( k = 0; k < nh; k++ )
			{
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					rotate( &soft[k*c_ldpc], rsoft, circ, sizeof(soft[0]), c_ldpc ); 

					for( n = 0; n < c_ldpc; n++ )
						synd[stateOffset+n] ^= rsoft[n] < 0;

					for( n = 0; n < c_ldpc; n++ )
						buffer[n] = dcs[stateOffset+n].pos == k ? dcs[stateOffset+n].min2 : dcs[stateOffset+n].min1;

					for( n = 0; n < c_ldpc; n++ )
					{
						short sign   = BnNS[memOffset+n] ^ dcs[stateOffset+n].sign;
#ifdef MS_MUL_CORRECTION
						IMS_DATA val = (IMS_DATA)((buffer[n] * ialpha) >> MS_ALPHA_FPP);
#else
						IMS_DATA tmp = buffer[n] - delta;
						IMS_DATA val = tmp < 0 ? 0 : tmp;
#endif
						IMS_DATA t   = sign ? -val : val;
                        IMS_DATA v2c_msg = rsoft[n] - t;
						short v2c_sign = v2c_msg < 0;
						
						BnNS[memOffset+n] = v2c_sign;
						tmps[n].sign     ^= v2c_sign;

						val = v2c_msg < 0 ? -v2c_msg : v2c_msg;
						val = (val > max_data) ? max_data : val;


						if( val < tmps[n].min1 )      
						{
							tmps[n].pos = k;
							tmps[n].min2 = tmps[n].min1; 
							tmps[n].min1 = val;
						}
						else
						{
							if( val < tmps[n].min2 )
								tmps[n].min2 = val;
						}

					}  
				}


				memOffset += c_ldpc;
			}

			for( k = 0; k < c_ldpc; k++ )
				dcs[stateOffset+k] = tmps[k];

		}

		for( parity = 0, i = 0; i < r_ldpc; i++ )
			parity |= synd[i];


		if( parity == 0 )
			break; 

	} // stage

	return parity ? -iter : iter+1; 
}




int syndrome_graph_cycle( short **hb, int rh, int nh, int m, short qhard[], short synd[], short buf[], short **hc_rl, short **mul/*, int rweight[]*/ )
{
	// Computes syndrome using precomputed multiplication table
	int i, j, k;
	int jj;

	for( i = 0; i < rh; i++ )
	{
		for( k = 0; k < m; k++ )
			synd[k] = 0;

		jj = 0;
		for( j = 0; j < nh; j++ )
		{
			int circ = hb[i][j];

			if( circ != -1 )
			{
				int symb  = hc_rl[i][jj];
				int pos_n = j * m;			
				
				rotate( &qhard[pos_n], buf, circ, sizeof(qhard[0]), m );

				for( k = 0; k < m; k++ )
#ifdef ORIG_TABLES
					synd[k] ^= mul[ buf[k] ][symb];
#else
					synd[k] ^= mul[symb][ buf[k] ];
#endif
				jj++;
			}
		}

		for( k = 0; k < m; k++ )
		{
			if( synd[k] )
				return i*m+k+1; 
		}
	}

	return 0;
}




#if 01
void fast_hadamar_C_16( double x[], double y[], int q, int c )
{
	double t;
	memcpy( y, x, q * sizeof(x[0]) );
	
	t = y[ 0];  y[ 0] += y[ 1];  y[ 1] = t - y[ 1];
	t = y[ 2];  y[ 2] += y[ 3];  y[ 3] = t - y[ 3];
	t = y[ 4];  y[ 4] += y[ 5];  y[ 5] = t - y[ 5];
	t = y[ 6];  y[ 6] += y[ 7];  y[ 7] = t - y[ 7];
	t = y[ 8];  y[ 8] += y[ 9];  y[ 9] = t - y[ 9];
	t = y[10];  y[10] += y[11];  y[11] = t - y[11];
	t = y[12];  y[12] += y[13];  y[13] = t - y[13];
	t = y[14];  y[14] += y[15];  y[15] = t - y[15];
	
	t = y[ 0];  y[ 0] += y[ 2];  y[ 2] = t - y[ 2];
	t = y[ 1];  y[ 1] += y[ 3];  y[ 3] = t - y[ 3];
	t = y[ 4];  y[ 4] += y[ 6];  y[ 6] = t - y[ 6];
	t = y[ 5];  y[ 5] += y[ 7];  y[ 7] = t - y[ 7];
	t = y[ 8];  y[ 8] += y[10];  y[10] = t - y[10];
	t = y[ 9];  y[ 9] += y[11];  y[11] = t - y[11];
	t = y[12];  y[12] += y[14];  y[14] = t - y[14];
	t = y[13];  y[13] += y[15];  y[15] = t - y[15];
	
	t = y[ 0];  y[ 0] += y[ 4];  y[ 4] = t - y[ 4];
	t = y[ 1];  y[ 1] += y[ 5];  y[ 5] = t - y[ 5];
	t = y[ 2];  y[ 2] += y[ 6];  y[ 6] = t - y[ 6];
	t = y[ 3];  y[ 3] += y[ 7];  y[ 7] = t - y[ 7];
	t = y[ 8];  y[ 8] += y[12];  y[12] = t - y[12];
	t = y[ 9];  y[ 9] += y[13];  y[13] = t - y[13];
	t = y[10];  y[10] += y[14];  y[14] = t - y[14];
	t = y[11];  y[11] += y[15];  y[15] = t - y[15];
	
	t = y[ 0];  y[ 0] += y[ 8];  y[ 8] = t - y[ 8];
	t = y[ 1];  y[ 1] += y[ 9];  y[ 9] = t - y[ 9];
	t = y[ 2];  y[ 2] += y[10];  y[10] = t - y[10];
	t = y[ 3];  y[ 3] += y[11];  y[11] = t - y[11];
	t = y[ 4];  y[ 4] += y[12];  y[12] = t - y[12];
	t = y[ 5];  y[ 5] += y[13];  y[13] = t - y[13];
	t = y[ 6];  y[ 6] += y[14];  y[14] = t - y[14];
	t = y[ 7];  y[ 7] += y[15];  y[15] = t - y[15];
}


void fast_hadamar_C_64( double x[], double y[], int q, int c )
{
	double t;
	memcpy( y, x, q * sizeof(x[0]) );
	t = y[ 0];  y[ 0] += y[ 1];  y[ 1] = t - y[ 1];
	t = y[ 2];  y[ 2] += y[ 3];  y[ 3] = t - y[ 3];
	t = y[ 4];  y[ 4] += y[ 5];  y[ 5] = t - y[ 5];
	t = y[ 6];  y[ 6] += y[ 7];  y[ 7] = t - y[ 7];
	t = y[ 8];  y[ 8] += y[ 9];  y[ 9] = t - y[ 9];
	t = y[10];  y[10] += y[11];  y[11] = t - y[11];
	t = y[12];  y[12] += y[13];  y[13] = t - y[13];
	t = y[14];  y[14] += y[15];  y[15] = t - y[15];
	t = y[16];  y[16] += y[17];  y[17] = t - y[17];
	t = y[18];  y[18] += y[19];  y[19] = t - y[19];
	t = y[20];  y[20] += y[21];  y[21] = t - y[21];
	t = y[22];  y[22] += y[23];  y[23] = t - y[23];
	t = y[24];  y[24] += y[25];  y[25] = t - y[25];
	t = y[26];  y[26] += y[27];  y[27] = t - y[27];
	t = y[28];  y[28] += y[29];  y[29] = t - y[29];
	t = y[30];  y[30] += y[31];  y[31] = t - y[31];
	t = y[32];  y[32] += y[33];  y[33] = t - y[33];
	t = y[34];  y[34] += y[35];  y[35] = t - y[35];
	t = y[36];  y[36] += y[37];  y[37] = t - y[37];
	t = y[38];  y[38] += y[39];  y[39] = t - y[39];
	t = y[40];  y[40] += y[41];  y[41] = t - y[41];
	t = y[42];  y[42] += y[43];  y[43] = t - y[43];
	t = y[44];  y[44] += y[45];  y[45] = t - y[45];
	t = y[46];  y[46] += y[47];  y[47] = t - y[47];
	t = y[48];  y[48] += y[49];  y[49] = t - y[49];
	t = y[50];  y[50] += y[51];  y[51] = t - y[51];
	t = y[52];  y[52] += y[53];  y[53] = t - y[53];
	t = y[54];  y[54] += y[55];  y[55] = t - y[55];
	t = y[56];  y[56] += y[57];  y[57] = t - y[57];
	t = y[58];  y[58] += y[59];  y[59] = t - y[59];
	t = y[60];  y[60] += y[61];  y[61] = t - y[61];
	t = y[62];  y[62] += y[63];  y[63] = t - y[63];
	t = y[ 0];  y[ 0] += y[ 2];  y[ 2] = t - y[ 2];
	t = y[ 1];  y[ 1] += y[ 3];  y[ 3] = t - y[ 3];
	t = y[ 4];  y[ 4] += y[ 6];  y[ 6] = t - y[ 6];
	t = y[ 5];  y[ 5] += y[ 7];  y[ 7] = t - y[ 7];
	t = y[ 8];  y[ 8] += y[10];  y[10] = t - y[10];
	t = y[ 9];  y[ 9] += y[11];  y[11] = t - y[11];
	t = y[12];  y[12] += y[14];  y[14] = t - y[14];
	t = y[13];  y[13] += y[15];  y[15] = t - y[15];
	t = y[16];  y[16] += y[18];  y[18] = t - y[18];
	t = y[17];  y[17] += y[19];  y[19] = t - y[19];
	t = y[20];  y[20] += y[22];  y[22] = t - y[22];
	t = y[21];  y[21] += y[23];  y[23] = t - y[23];
	t = y[24];  y[24] += y[26];  y[26] = t - y[26];
	t = y[25];  y[25] += y[27];  y[27] = t - y[27];
	t = y[28];  y[28] += y[30];  y[30] = t - y[30];
	t = y[29];  y[29] += y[31];  y[31] = t - y[31];
	t = y[32];  y[32] += y[34];  y[34] = t - y[34];
	t = y[33];  y[33] += y[35];  y[35] = t - y[35];
	t = y[36];  y[36] += y[38];  y[38] = t - y[38];
	t = y[37];  y[37] += y[39];  y[39] = t - y[39];
	t = y[40];  y[40] += y[42];  y[42] = t - y[42];
	t = y[41];  y[41] += y[43];  y[43] = t - y[43];
	t = y[44];  y[44] += y[46];  y[46] = t - y[46];
	t = y[45];  y[45] += y[47];  y[47] = t - y[47];
	t = y[48];  y[48] += y[50];  y[50] = t - y[50];
	t = y[49];  y[49] += y[51];  y[51] = t - y[51];
	t = y[52];  y[52] += y[54];  y[54] = t - y[54];
	t = y[53];  y[53] += y[55];  y[55] = t - y[55];
	t = y[56];  y[56] += y[58];  y[58] = t - y[58];
	t = y[57];  y[57] += y[59];  y[59] = t - y[59];
	t = y[60];  y[60] += y[62];  y[62] = t - y[62];
	t = y[61];  y[61] += y[63];  y[63] = t - y[63];
	t = y[ 0];  y[ 0] += y[ 4];  y[ 4] = t - y[ 4];
	t = y[ 1];  y[ 1] += y[ 5];  y[ 5] = t - y[ 5];
	t = y[ 2];  y[ 2] += y[ 6];  y[ 6] = t - y[ 6];
	t = y[ 3];  y[ 3] += y[ 7];  y[ 7] = t - y[ 7];
	t = y[ 8];  y[ 8] += y[12];  y[12] = t - y[12];
	t = y[ 9];  y[ 9] += y[13];  y[13] = t - y[13];
	t = y[10];  y[10] += y[14];  y[14] = t - y[14];
	t = y[11];  y[11] += y[15];  y[15] = t - y[15];
	t = y[16];  y[16] += y[20];  y[20] = t - y[20];
	t = y[17];  y[17] += y[21];  y[21] = t - y[21];
	t = y[18];  y[18] += y[22];  y[22] = t - y[22];
	t = y[19];  y[19] += y[23];  y[23] = t - y[23];
	t = y[24];  y[24] += y[28];  y[28] = t - y[28];
	t = y[25];  y[25] += y[29];  y[29] = t - y[29];
	t = y[26];  y[26] += y[30];  y[30] = t - y[30];
	t = y[27];  y[27] += y[31];  y[31] = t - y[31];
	t = y[32];  y[32] += y[36];  y[36] = t - y[36];
	t = y[33];  y[33] += y[37];  y[37] = t - y[37];
	t = y[34];  y[34] += y[38];  y[38] = t - y[38];
	t = y[35];  y[35] += y[39];  y[39] = t - y[39];
	t = y[40];  y[40] += y[44];  y[44] = t - y[44];
	t = y[41];  y[41] += y[45];  y[45] = t - y[45];
	t = y[42];  y[42] += y[46];  y[46] = t - y[46];
	t = y[43];  y[43] += y[47];  y[47] = t - y[47];
	t = y[48];  y[48] += y[52];  y[52] = t - y[52];
	t = y[49];  y[49] += y[53];  y[53] = t - y[53];
	t = y[50];  y[50] += y[54];  y[54] = t - y[54];
	t = y[51];  y[51] += y[55];  y[55] = t - y[55];
	t = y[56];  y[56] += y[60];  y[60] = t - y[60];
	t = y[57];  y[57] += y[61];  y[61] = t - y[61];
	t = y[58];  y[58] += y[62];  y[62] = t - y[62];
	t = y[59];  y[59] += y[63];  y[63] = t - y[63];
	t = y[ 0];  y[ 0] += y[ 8];  y[ 8] = t - y[ 8];
	t = y[ 1];  y[ 1] += y[ 9];  y[ 9] = t - y[ 9];
	t = y[ 2];  y[ 2] += y[10];  y[10] = t - y[10];
	t = y[ 3];  y[ 3] += y[11];  y[11] = t - y[11];
	t = y[ 4];  y[ 4] += y[12];  y[12] = t - y[12];
	t = y[ 5];  y[ 5] += y[13];  y[13] = t - y[13];
	t = y[ 6];  y[ 6] += y[14];  y[14] = t - y[14];
	t = y[ 7];  y[ 7] += y[15];  y[15] = t - y[15];
	t = y[16];  y[16] += y[24];  y[24] = t - y[24];
	t = y[17];  y[17] += y[25];  y[25] = t - y[25];
	t = y[18];  y[18] += y[26];  y[26] = t - y[26];
	t = y[19];  y[19] += y[27];  y[27] = t - y[27];
	t = y[20];  y[20] += y[28];  y[28] = t - y[28];
	t = y[21];  y[21] += y[29];  y[29] = t - y[29];
	t = y[22];  y[22] += y[30];  y[30] = t - y[30];
	t = y[23];  y[23] += y[31];  y[31] = t - y[31];
	t = y[32];  y[32] += y[40];  y[40] = t - y[40];
	t = y[33];  y[33] += y[41];  y[41] = t - y[41];
	t = y[34];  y[34] += y[42];  y[42] = t - y[42];
	t = y[35];  y[35] += y[43];  y[43] = t - y[43];
	t = y[36];  y[36] += y[44];  y[44] = t - y[44];
	t = y[37];  y[37] += y[45];  y[45] = t - y[45];
	t = y[38];  y[38] += y[46];  y[46] = t - y[46];
	t = y[39];  y[39] += y[47];  y[47] = t - y[47];
	t = y[48];  y[48] += y[56];  y[56] = t - y[56];
	t = y[49];  y[49] += y[57];  y[57] = t - y[57];
	t = y[50];  y[50] += y[58];  y[58] = t - y[58];
	t = y[51];  y[51] += y[59];  y[59] = t - y[59];
	t = y[52];  y[52] += y[60];  y[60] = t - y[60];
	t = y[53];  y[53] += y[61];  y[61] = t - y[61];
	t = y[54];  y[54] += y[62];  y[62] = t - y[62];
	t = y[55];  y[55] += y[63];  y[63] = t - y[63];
	t = y[ 0];  y[ 0] += y[16];  y[16] = t - y[16];
	t = y[ 1];  y[ 1] += y[17];  y[17] = t - y[17];
	t = y[ 2];  y[ 2] += y[18];  y[18] = t - y[18];
	t = y[ 3];  y[ 3] += y[19];  y[19] = t - y[19];
	t = y[ 4];  y[ 4] += y[20];  y[20] = t - y[20];
	t = y[ 5];  y[ 5] += y[21];  y[21] = t - y[21];
	t = y[ 6];  y[ 6] += y[22];  y[22] = t - y[22];
	t = y[ 7];  y[ 7] += y[23];  y[23] = t - y[23];
	t = y[ 8];  y[ 8] += y[24];  y[24] = t - y[24];
	t = y[ 9];  y[ 9] += y[25];  y[25] = t - y[25];
	t = y[10];  y[10] += y[26];  y[26] = t - y[26];
	t = y[11];  y[11] += y[27];  y[27] = t - y[27];
	t = y[12];  y[12] += y[28];  y[28] = t - y[28];
	t = y[13];  y[13] += y[29];  y[29] = t - y[29];
	t = y[14];  y[14] += y[30];  y[30] = t - y[30];
	t = y[15];  y[15] += y[31];  y[31] = t - y[31];
	t = y[32];  y[32] += y[48];  y[48] = t - y[48];
	t = y[33];  y[33] += y[49];  y[49] = t - y[49];
	t = y[34];  y[34] += y[50];  y[50] = t - y[50];
	t = y[35];  y[35] += y[51];  y[51] = t - y[51];
	t = y[36];  y[36] += y[52];  y[52] = t - y[52];
	t = y[37];  y[37] += y[53];  y[53] = t - y[53];
	t = y[38];  y[38] += y[54];  y[54] = t - y[54];
	t = y[39];  y[39] += y[55];  y[55] = t - y[55];
	t = y[40];  y[40] += y[56];  y[56] = t - y[56];
	t = y[41];  y[41] += y[57];  y[57] = t - y[57];
	t = y[42];  y[42] += y[58];  y[58] = t - y[58];
	t = y[43];  y[43] += y[59];  y[59] = t - y[59];
	t = y[44];  y[44] += y[60];  y[60] = t - y[60];
	t = y[45];  y[45] += y[61];  y[61] = t - y[61];
	t = y[46];  y[46] += y[62];  y[62] = t - y[62];
	t = y[47];  y[47] += y[63];  y[63] = t - y[63];
	t = y[ 0];  y[ 0] += y[32];  y[32] = t - y[32];
	t = y[ 1];  y[ 1] += y[33];  y[33] = t - y[33];
	t = y[ 2];  y[ 2] += y[34];  y[34] = t - y[34];
	t = y[ 3];  y[ 3] += y[35];  y[35] = t - y[35];
	t = y[ 4];  y[ 4] += y[36];  y[36] = t - y[36];
	t = y[ 5];  y[ 5] += y[37];  y[37] = t - y[37];
	t = y[ 6];  y[ 6] += y[38];  y[38] = t - y[38];
	t = y[ 7];  y[ 7] += y[39];  y[39] = t - y[39];
	t = y[ 8];  y[ 8] += y[40];  y[40] = t - y[40];
	t = y[ 9];  y[ 9] += y[41];  y[41] = t - y[41];
	t = y[10];  y[10] += y[42];  y[42] = t - y[42];
	t = y[11];  y[11] += y[43];  y[43] = t - y[43];
	t = y[12];  y[12] += y[44];  y[44] = t - y[44];
	t = y[13];  y[13] += y[45];  y[45] = t - y[45];
	t = y[14];  y[14] += y[46];  y[46] = t - y[46];
	t = y[15];  y[15] += y[47];  y[47] = t - y[47];
	t = y[16];  y[16] += y[48];  y[48] = t - y[48];
	t = y[17];  y[17] += y[49];  y[49] = t - y[49];
	t = y[18];  y[18] += y[50];  y[50] = t - y[50];
	t = y[19];  y[19] += y[51];  y[51] = t - y[51];
	t = y[20];  y[20] += y[52];  y[52] = t - y[52];
	t = y[21];  y[21] += y[53];  y[53] = t - y[53];
	t = y[22];  y[22] += y[54];  y[54] = t - y[54];
	t = y[23];  y[23] += y[55];  y[55] = t - y[55];
	t = y[24];  y[24] += y[56];  y[56] = t - y[56];
	t = y[25];  y[25] += y[57];  y[57] = t - y[57];
	t = y[26];  y[26] += y[58];  y[58] = t - y[58];
	t = y[27];  y[27] += y[59];  y[59] = t - y[59];
	t = y[28];  y[28] += y[60];  y[60] = t - y[60];
	t = y[29];  y[29] += y[61];  y[61] = t - y[61];
	t = y[30];  y[30] += y[62];  y[62] = t - y[62];
	t = y[31];  y[31] += y[63];  y[63] = t - y[63];
}
#else
void fast_hadamar_C_16( double x[], double y[], int q, int c )
{
	int i;
	double t;
	int indx[] =
	{
		0,  1,		2,  3, 		4,  5,		6,  7,
		8,  9,		10, 11,		12, 13,		14, 15,
		
		0,  2,		1,  3,		4,  6,		5,  7,
		8, 10,		9, 11,		12, 14,		13, 15,
		
		0,  4,		1,  5,		2,  6,		3,  7,
		8, 12,		9, 13,		10, 14,		11, 15,

		0,  8,		1,  9,		2, 10,		3, 11,
		4, 12,		5, 13,		6, 14,		7, 15,
	};
	int steps = sizeof( indx) /  sizeof(indx[0]) / 2;
	int j0, j1, k;

	memcpy( y, x, q * sizeof(x[0]) );

	for( i = k = 0; i < steps; i++, k+= 2 )
	{
		j0 = indx[k]; j1 = indx[k+1];
		t = y[j0];              //t=y(j0+1);
		y[j0] += y[j1];         //y(j0+1)=y(j0+1)+y(j1+1);
		y[j1] = t - y[j1];      //y(j1+1)=t-y(j1+1);
	}
}


void fast_hadamar_C_64( double x[], double y[], int q, int c )
{
	int i;
	double t;
	int indx[] =
	{
		0,  1,		  2,  3,		  4,  5,		  6,  7,
		8,  9,		 10, 11,		 12, 13,		 14, 15,
		16, 17,		 18, 19,		 20, 21,		 22, 23,
		24, 25,		 26, 27,		 28, 29,		 30, 31,
		32, 33,		 34, 35,		 36, 37,		 38, 39,
		40, 41,		 42, 43,		 44, 45,		 46, 47,
		48, 49,		 50, 51,		 52, 53,		 54, 55,
		56, 57,		 58, 59,		 60, 61,		 62, 63,
		0,  2,		  1,  3,		  4,  6,		  5,  7,
		8, 10,		  9, 11,		 12, 14,		 13, 15,
		16, 18,		 17, 19,		 20, 22,		 21, 23,
		24, 26,		 25, 27,		 28, 30,		 29, 31,
		32, 34,		 33, 35,		 36, 38,		 37, 39,
		40, 42,		 41, 43,		 44, 46,		 45, 47,
		48, 50,		 49, 51,		 52, 54,		 53, 55,
		56, 58,		 57, 59,		 60, 62,		 61, 63,
		0,  4,		  1,  5,		  2,  6,		  3,  7,
		8, 12,		  9, 13,		 10, 14,		 11, 15,
		16, 20,		 17, 21,		 18, 22,		 19, 23,
		24, 28,		 25, 29,		 26, 30,		 27, 31,
		32, 36,		 33, 37,		 34, 38,		 35, 39,
		40, 44,		 41, 45,		 42, 46,		 43, 47,
		48, 52,		 49, 53,		 50, 54,		 51, 55,
		56, 60,		 57, 61,		 58, 62,		 59, 63,
		0,  8,		  1,  9,		  2, 10,		  3, 11,
		4, 12,		  5, 13,		  6, 14,		  7, 15,
		16, 24,		 17, 25,		 18, 26,		 19, 27,
		20, 28,		 21, 29,		 22, 30,		 23, 31,
		32, 40,		 33, 41,		 34, 42,		 35, 43,
		36, 44,		 37, 45,		 38, 46,		 39, 47,
		48, 56,		 49, 57,		 50, 58,		 51, 59,
		52, 60,		 53, 61,		 54, 62,		 55, 63,
		0, 16,		  1, 17,		  2, 18,		  3, 19,
		4, 20,		  5, 21,		  6, 22,		  7, 23,
		8, 24,		  9, 25,		 10, 26,		 11, 27,
		12, 28,		 13, 29,		 14, 30,		 15, 31,
		32, 48,		 33, 49,		 34, 50,		 35, 51,
		36, 52,		 37, 53,		 38, 54,		 39, 55,
		40, 56,		 41, 57,		 42, 58,		 43, 59,
		44, 60,		 45, 61,		 46, 62,		 47, 63,
		0, 32,		  1, 33,		  2, 34,		  3, 35,
		4, 36,		  5, 37,		  6, 38,		  7, 39,
		8, 40,		  9, 41,		 10, 42,		 11, 43,
		12, 44,		 13, 45,		 14, 46,		 15, 47,
		16, 48,		 17, 49,		 18, 50,		 19, 51,
		20, 52,		 21, 53,		 22, 54,		 23, 55,
		24, 56,		 25, 57,		 26, 58,		 27, 59,
		28, 60,		 29, 61,		 30, 62,		 31, 63,
	};
	int steps = sizeof( indx) /  sizeof(indx[0]) / 2;
	int j0, j1, k;

	memcpy( y, x, q * sizeof(x[0]) );

	for( i = k = 0; i < steps; i++, k+= 2 )
	{
		j0 = indx[k]; j1 = indx[k+1];
		t = y[j0];              //t=y(j0+1);
		y[j0] += y[j1];         //y(j0+1)=y(j0+1)+y(j1+1);
		y[j1] = t - y[j1];      //y(j1+1)=t-y(j1+1);
	}
}
#endif

void fast_hadamar_C( double x[], double y[], int q, int c )
{
	int f, i, j;
	int lim1 = c;
	int lim2 = q/2;
	int jh, j1, j0;
	double t;

	f = 1;
	memcpy( y, x, q * sizeof(x[0]) );

	for( i = 0; i < c; i++ )
	{
		for( j = 0; j < lim2; j++ )
		{
			jh = j >> i;            //jh=bitshift(j,-i);
			j1 = j & ( f-1 );       //jl=bitand(j,f-1);
			j0 = (jh << (i+1)) + j1;  //j0=bitshift(jh,i+1)+jl;
			j1 = j0 + f;            //j1=j0+f;
			t = y[j0];              //t=y(j0+1);
			y[j0] += y[j1];         //y(j0+1)=y(j0+1)+y(j1+1);
			y[j1] = t - y[j1];      //y(j1+1)=t-y(j1+1);
		}
		f += f;                     //f = f+f;
	}
}


void hadamar( double x[], double y[], int q, int c )
{
	switch( q )
	{
	case 16 : fast_hadamar_C_16( x, y, q, c ); break;
	case 64 : fast_hadamar_C_64( x, y, q, c ); break;
	default:  fast_hadamar_C( x, y, q, c ); 
	}
}

void double_int_double( double x[], int n, int fpp )
{
	int i;
	int one = 1 << fpp;
	for( i = 0; i < n; i++ )
	{
		int    sign = x[i] < 0 ? 1 : 0;
		double val  = x[i] < 0 ? -x[i] : x[i];
		int t = (int)(val * one + 0.5);
//		if( t == 0 )
//			t = 1;
		t = sign ? -t : t;
		x[i] = (double)t / (double)one;
	}
}

double log_to_double( i32 slog, int fpp ) 
{
	double res;
	double one;
	if( slog < 0 )
		return 0;
	else
	{
		i32  shift = slog - fpp;

		if( shift > 0 )
			return 1 < shift;

		res = 1.0;
		one = 1.0 / (double)(1L << 16);
		
		while( shift < 0 )
		{
			shift += 16;
			res *= one;
		}
		res *= 1 << shift;

		return res;
	}
}

typedef struct  
{
	int val;
	int shift;
} FHT_DATA;

FHT_DATA double2data( double x, int fpp )
{
	FHT_DATA res;
	int sign, absval, shift, size;
	int one = 1 << fpp;

	if( x < 0 )
	{
		sign = 1;
		absval = (int)(-x * one);
	}
	else
	{
		sign = 0;
		absval = (int)(x * one);
	}

	if( absval )
	{
		size = get_bitsize( absval );
		shift = fpp - size;
		absval <<= shift;
	}
	else
	{
		shift = 0;
	}

	res.val = sign ? -absval : absval;
	res.shift = -shift;

	return res;
}

FHT_DATA mul_fht_data( FHT_DATA x, FHT_DATA y, int fpp )
{
	FHT_DATA res;

	res.val = (int)((((i64)x.val * (i64)y.val) << 1) >> fpp);
	res.shift = x.shift + y.shift - 1;
	return res;
}

FHT_DATA mul_msb_fht_data( FHT_DATA x, FHT_DATA y, int fpp, int msb )
{
	FHT_DATA res;
	int shift;
	i64 r;
	int a = x.val;
	int b = y.val;
	
	shift = msb < fpp ? fpp - msb : 0;
	
	a >>= shift;
	b >>= shift;

	res.shift = x.shift + y.shift - 1;

	r = ((i64)a * (i64)b) << 1;
	r <<= (shift+shift);
	res.val = (int)(r >> fpp);

	return res;
}

double fht_data2double( FHT_DATA x, int fpp )
{
	int shift;
	double one = (double)(1 << fpp );
	double res = x.val / one;
	if( x.shift < 0 )
	{
		shift = -x.shift;
		while(shift > 16 )
		{
			res /= (double)(1 << 16);
			shift -= 16;
		}
		res /= (double)(1 << shift);

	}
	else
	{
		shift = x.shift;
		while(shift > 16 )
		{
			res *= (double)(1 << 16);
			shift -= 16;
		}
		res *= (double)(1 << shift);
	}

	return res;
}


void  map_graph
( 
	double **soft_in, 
	double **soft_outs, 
	int pos_r, 
	int step, 
	short *rl, 
	int q_bits, 
	int rw, 
	short **mul, 
	short **div, 
	short *mask,
	short *hard,
	double *HAD
)
{

	int i, j, k;
	// ARRAYS
	static double S_probH[RWMAX][QMAX];
	static double Sigma_forwardH[RWMAX][QMAX];
	static double Sigma_backH[RWMAX][QMAX]; // Hadamar domain
	static FHT_DATA S_prob[RWMAX][QMAX];
	static FHT_DATA Sigma_f[RWMAX][QMAX];
	static FHT_DATA Sigma_b[RWMAX][QMAX];;
	static double s[QMAX];
	int binlogq = q_bits;
	int q = 1 << q_bits;
	int idx;
	double qinv = 1.0 / (double)q;


	// Permutation
	for( j = 0; j < rw; j++ )
	{
		if( mask[j] == 0 )
		{
			for( i = 0; i < q; i++ )
#ifdef ORIG_TABLES
				s[i] = soft_in[ div[i][rl[j]] ][j*step + pos_r];  
#else
				s[i] = soft_in[ div[rl[j]][i] ][j*step + pos_r];  
#endif

#if 1
			normalize_abs( s, q, 1,  HAD_FPP );
			double_int_double( s, q, HAD_FPP );
#endif
			hadamar( s, S_probH[j], q, binlogq );
		}
		else
		{
#ifdef ORIG_TABLES
			short t = mul[hard[j]][rl[j]];
#else
			short t = mul[rl[j]][hard[j]];
#endif
			for( i = 0; i < q; i++ )
				S_probH[j][i] = HAD[t * q + i];
		}

		normalize_abs( S_probH[j], q, 1, MAP_FPP );
		double_int_double( S_probH[j], q, MAP_FPP );
	}



#ifndef MAP_GRAPH_USE_INT
	// Alpha
	for( i = 0; i < q; i++ )
		Sigma_forwardH[0][i] = S_probH[0][i];

	for( k = 1; k < rw-1; k++ )     //loop over nonzero symbols of the check                                
		for( i = 0; i < q; i++ )
			Sigma_forwardH[k][i] = S_probH[k][i] * Sigma_forwardH[k-1][i];

	// Beta
	for( i = 0; i < q; i++ )
		Sigma_backH[rw-1][i] = S_probH[rw-1][i];

	for( k = rw-2; k > 0; k-- ) //loop over nonzero symbols of the check
		for( i = 0; i < q; i++ )
			Sigma_backH[k][i] = S_probH[k][i] * Sigma_backH[k+1][i];

#else
	for( j = 0; j < rw; j++ )
		for( i = 0; i < q; i++ )
			S_prob[j][i] = double2data( S_probH[j][i], MAP_FPP );

	for( i = 0; i < q; i++ )
		Sigma_f[0][i] = double2data( S_probH[0][i], MAP_FPP );

	for( k = 1; k < rw-1; k++ )     //loop over nonzero symbols of the check                                
		for( i = 0; i < q; i++ )
			Sigma_f[k][i] = mul_msb_fht_data( S_prob[k][i], Sigma_f[k-1][i], MAP_FPP, MAP_MULT );

	for( i = 0; i < q; i++ )
		Sigma_b[rw-1][i] = double2data( S_probH[rw-1][i], MAP_FPP );

	for( k = rw-2; k > 0; k-- ) //loop over nonzero symbols of the check
		for( i = 0; i < q; i++ )
			Sigma_b[k][i] = mul_msb_fht_data( S_prob[k][i], Sigma_b[k+1][i], MAP_FPP, MAP_MULT );

	for( k = 0; k < rw-1; k++ )     //loop over nonzero symbols of the check                                
		for( i = 0; i < q; i++ )
			Sigma_forwardH[k][i] = fht_data2double( Sigma_f[k][i], MAP_FPP );

	for( k = 1; k < rw; k++ )     //loop over nonzero symbols of the check                                
		for( i = 0; i < q; i++ )
			Sigma_backH[k][i] = fht_data2double( Sigma_b[k][i], MAP_FPP );
#endif


    
	// Sigma
	idx = pos_r;
	for( j = 0; j < rw; j++ )
	{
		if( mask[j] == 0 )
		{
			double *zptr;
			double Z[QMAX];
			
			if( j == 0 )
				zptr = &Sigma_backH[1][0];
			else
			{
				if( j == rw-1 )
					zptr = &Sigma_forwardH[rw-2][0];
				else
				{
					for( i = 0; i < q; i++ )
					{
#ifndef MAP_GRAPH_USE_INT
						Z[i] = Sigma_forwardH[j-1][i] * Sigma_backH[j+1][i];
#else
						FHT_DATA t01 = mul_msb_fht_data( Sigma_f[j-1][i], Sigma_b[j+1][i], MAP_FPP, MAP_MULT );
						Z[i] = fht_data2double( t01, MAP_FPP );
#endif
					}

					zptr = Z;
				}
			}

#ifdef SCALABLE
			normalize_abs( zptr, q, 1, HAD_FPP );
			double_int_double( zptr, q, HAD_FPP );
#endif

			hadamar( zptr, s, q, binlogq );


#if 01
			for( i = 0; i < q; i++ )
				if( s[i] < 0 )
					s[i] = 0;

			normalize( s, q, 1, PROB_FPP );
			double_int_double( s, q, PROB_FPP );

			for( i = 0; i < q; i++ )
			{
#ifdef ORIG_TABLES
				soft_outs[i][idx] = s[ mul[i][rl[j]] ];
#else
				soft_outs[i][idx] = s[ mul[rl[j]][i] ];
#endif
			}


#else
			for( i = 0; i < q; i++ )
			{
#ifdef ORIG_TABLES
				soft_outs[i][idx] = s[ mul[i][rl[j]] ] / q;
#else
				soft_outs[i][idx] = s[ mul[rl[j]][i] ] * qinv;
#endif
			}
#endif

		}

		idx += step;
	}
}


int plength( int pol[] )
{
	int i = 0;
	do
	i++;
	while( pol[i] );

	return i;
}

int pol_bank( int c[], int *pol, int *M )
{
	static int prim_polinomial[10][100] =
	{
		{03, 0},
		{07, 0},
		{015, 0},
		{023,037, 0},  
		{045, 075, 067, 0},   
		{0103,0147,0155, 0}, 
		{0203, 0211, 0217, 0235, 0367, 0277, 0325,  0313, 0345, 0},   
		{0435, 0551, 0453, 0545, 0543, 0537, 0703, 0747, 0}, 
		{
			01021, 01131, 01461, 01423, 01055, 01167, 01541, 01333, 
			01605, 01751, 01743, 01617, 01553, 01157, 01715, 01563, 
			01713, 01175, 01725, 01225, 01275, 01773, 01425, 01267, 0
		},
		{
			02011, 02415, 03771, 02157, 03515, 02773, 02033, 02443, 
			02461, 03023, 03543, 02745, 02431, 03177, 03525, 02617, 
			03471, 03323, 03507, 03623, 02707, 02327, 03265, 02055, 
			03575, 03171, 02047, 03025, 03337, 03211, 0
		}
	};	
	int *prim;
	if( c[0] > 10 )
	{
		*pol = 1;
		*M = 1;
		return 0;
	}

	prim = prim_polinomial[c[0] - 1];
	*M = plength( prim );

	*pol = prim[0];
	if( c[1] != 0  && c[1] < *M )
		*pol = prim[c[1] - 1];

	return 1;
}

int gf2( int m, int p, short *gf_log, short *gf_alog )
{
	// v is field arrray of gf(2^m) over prim polynom p
	int i;
	int m2 = 1 << m;
	int element;
/*	
	if( nargin == 1 )
	{
		P=[3,7,11,19,37,67,137, 285, 529, 1033];
		p=P(m);
	}
*/

	element = 1;
	gf_log[0] = -1; 
	gf_alog[0] = 1;

	for( i = 1; i < m2; i++ )
	{
		gf_alog[i-1] = element;
		gf_log[element] = i-1;

		element <<= 1;
		if( element >= m2 )
			element ^= p;

		if( element == 1 && i != (m2-1))
		{
			printf("p is not primitive\n");
			return -1;
		}
	}

	return 0;
}

static void p2table( short **mul, short **div, short *gf_log, short *gf_alog, int list[], int q_bits, int n_used )
{
	int i, j, res;
	int n = 1 << q_bits;
	int mod = n - 1;
	for( i = 0; i < n; i++ )
	{
		for( j = 0; j < n_used; j++ )
		{

			if( i == 0 )
#ifdef ORIG_TABLES
				mul[i][j] = 0;
#else
				mul[j][i] = 0;
#endif
			else
			{
				int x = gf_log[i];
				int y = gf_log[list[j]];

				res = x + y;
				res = res < mod ? res : res - mod;

#ifdef ORIG_TABLES
				mul[i][j] = gf_alog[res];
#else 
				mul[j][i] = gf_alog[res];
#endif


				res = x - y;
				res = res < 0 ? res + mod : res;
#ifdef ORIG_TABLES
				div[i][j] = gf_alog[res];
#else
				div[j][i] = gf_alog[res];
#endif
			}
/*
			if( i == 0 )
#ifdef ORIG_TABLES
				mul[i][j] = 0;
#else
				mul[j][i] = 0;
#endif
			else
			{
				int x = i;
				int y = gf_log[list[j]];

				res = x + y;
				res = res < mod ? res : res - mod;

#ifdef ORIG_TABLES
				mul[i][j] = gf_alog[res];
#else 
				mul[j][i] = res;
#endif


				res = x - y;
				res = res < 0 ? res + mod : res;
#ifdef ORIG_TABLES
				div[i][j] = gf_alog[res];
#else
				div[j][i] = res;
#endif
			}
*/
		}
	}
}

//#define SOFT_IN_OUTS
//#define MAKE_DUMP



double get_minMain( double x[], int flag[], int n )
{
	int i;
	double m;
	
	m = x[0];

	for( i = 0; i < n; i++ )
	{
		if( flag[i] == 0 && m > x[i] )
			m = x[i];
	}

	return m;
}


void sort( ELEMENT x[], int n, int k )
{
	int i, j;

	if( k >= n - 1 )
		k = n - 1;

	for( i = 0; i < k; i++ )
	{
		for( j = 0; j < n-1; j++ )
		{
			if( x[j].val > x[j+1].val )
			{
				ELEMENT t = x[j];
				x[j] = x[j+1];
				x[j+1] = t;
			}
		}
	}
}

void sort_ui32( ui32 x[], int n, int k )
{
	int i, j;

	if( k >= n - 1 )
		k = n - 1;

	for( i = 0; i < k; i++ )
	{
		for( j = 0; j < n-1; j++ )
		{
			if( x[j] > x[j+1] )
			{
				i32 t = x[j];
				x[j] = x[j+1];
				x[j+1] = t;
			}
		}
	}
}

#if 0 //#ifndef _WIN64
ui32 get_bitsize( ui32 d )
{
	ui32 res;
	if( d == 0 )
		return 0;

	__asm
	{ 
		mov edx, d;
		bsr eax, edx;
		mov res, eax;
	}

	return res + 1;
}
#else
ui32 get_bitsize( ui32 d )
{
	ui32 shift = 0;

	if( d == 0 )
		return 0;

	while( d  )
	{
		shift++;
		d >>= 1;
	}

	return shift;
}
#endif

#if defined( NORM_MAX )
// simple normalization by max value
void normalize( double *x, int n, int step )
{
	double c;
	int j;
	c = 0.0;
	n *= step;

	for( j = 0; j < n; j += step )
	{
		if( c < x[j] )
			c = x[j];
	}

	c *= 1.05;
	c = 1.0 / c;

	for( j = 0; j < n; j += step )
		x[j] *= c;
}
#elif defined( NORM_SHIFT )

void normalize( double *x, int n, int step, int fpp )
{

	double c;
	int j;
    int shift;
	int nshift;
    ui32 d = 1;
    ui32 dd;
	int one = 1 << fpp;

	c = 0.0;
	n *= step;

	for( j = 0; j < n; j += step )
	{
		ui32 a = (ui32)(x[j] * one);
        d |= a;
//        d |= ix[j];
	}
   
    dd = d;
    shift = 0;
    
 
	shift = get_bitsize( d );

    if( shift > 30 )
       shift = 30;
    
  
   
   nshift =  (fpp - shift);
   if( nshift < 0 )
   {
       c = 1.0 / (1 << -nshift);

	   for( j = 0; j < n; j += step )
		   x[j] *= c;
   }
   else
   {
       c = 1 << nshift;

	   for( j = 0; j < n; j += step )
		   x[j] *= c;
   }
}

void normalize_abs( double *x, int n, int step, int fpp )
{

	double c;
	int j;
	int shift;
	int nshift;
	ui32 d = 1;
	ui32 dd;
	int one = 1 << fpp;

	c = 0.0;
	n *= step;

	for( j = 0; j < n; j += step )
	{
		double val = x[j] < 0 ? -x[j] : x[j];
		ui32 a = (ui32)(val * one);
		d |= a;
		//        d |= ix[j];
	}

	dd = d;
	shift = 0;


	shift = get_bitsize( d );

	if( shift > 30 )
		shift = 30;



	nshift =  (fpp - shift);
	if( nshift < 0 )
	{
		c = 1.0 / (1 << -nshift);

		for( j = 0; j < n; j += step )
			x[j] *= c;
	}
	else
	{
		c = 1 << nshift;

		for( j = 0; j < n; j += step )
			x[j] *= c;
	}
}
/*
void inormalize( ui32 *x, int n, int step )
{

	double c;
	int j;
	int shift;
	int nshift;
	ui32 d = 1;
	ui32 dd;
	static int cnt=0;
	cnt++;

	n *= step;

	for( j = 0; j < n; j += step )
		d |= x[j];

	dd = d;
	shift = 0;


	shift = get_bitsize( d );

	if( shift > 30 )
		shift = 30;

	d = 1L << shift;


	nshift =  (PROB_FPP - shift);
	if( nshift < 0 )
	{
		for( j = 0; j < n; j += step )
			x[j] = (x[j] >> nshift) | 1;
	}
	else
	{
		for( j = 0; j < n; j += step )
			x[j] = (x[j] << nshift) | 1;
	}

}
*/
#else
// true normalization by sum
void normalize( double *x, int n, int step )
{
	double c;
	int j;
	c = 0.0;
	n *= step;
/*
	for( j = 0; j < n; j += step )
	{
		if( x[j] > 1.0 )
			j = j;
	}
*/
	for( j = 0; j < n; j += step )
		c += x[j];


	c = 1.0 / c;
	for( j = 0; j < n; j += step )
		x[j] *= c;
}
#endif

int sum_prod_gfq_decod_lm( DEC_STATE* st, double *soft[], short *qhard, double *decword[], int maxiter, double p_thr )
{
	int i, j, k;
	int iter;
	int syn;

	short *posh = st->fht_pos;
	short **hc = st->hc;
	short **hb = st->hb;
	short **hc_rl = st->fht_hc_rl;
	short *synd = st->syndr;
	short *buf = st->fht_buf;
	short **hb_ri = st->fht_hb_ri; 
	short **hb_ci = st->fht_hb_ci; 
	short **hb_cj = st->fht_hb_cj;
	short **mul_tab = st->fht_t;
	short **div_tab = st->fht_ti;
	double **soft_in  = st->fht_soft_in;
	double **soft_outs= st->fht_soft_outs;
	double **soft_out = st->fht_soft_out;
	double **buf0 = st->fht_buf0;
	double **buf1 = st->fht_buf1;
	double *HAD = st->fht_HAD;
	short *smask = st->fht_smask;
	int *rweight = st->fht_rw;
	int *list  = st->fht_list;
	int cw2 = st->fht_all_cw_2;
	int q_bits = st->q_bits;
	int max_rw = st->max_rw;
	int q  = st->q;
	int nh = st->nh;
	int rh = st->rh;
	int m  = st->m;
	int rg = rh * m;
	int ng = nh * m;
	int r = rg;
	int use_p_thr = p_thr == 0 ? 0 : 1;

	if( cw2 == 0 )
		return -1;


#ifdef SCALABLE 
	for( i = 0; i < ng; i++ )
	{
		double tmp[QMAX];
		for( j = 0; j < q; j++ )
			tmp[j] = soft[j][i];

		normalize( tmp, q, 1, INP_FPP );
		double_int_double( tmp, q, INP_FPP );

		for( j = 0; j < q; j++ )
			soft[j][i] = tmp[j];
	}
#endif

	// Initialize each constituent code symbol by its input soft data 

	for( j = 0; j < rh; j++ )
	{
		int cnt = 0;

		for( i = 0; i < nh; i++ )
		{
			int pos_n = i * m;
			int pos_r = j * m;
			int circ = hb[j][i];

			if( circ != -1 )
			{
				for( k = 0; k < q; k++ )
					rotate( &soft[k][pos_n], &soft_in[k][cnt*r+pos_r], circ, sizeof(soft[0][0]), m );
				cnt += 1;
			}
		}
	}


	// just to compute syndrome before iterations
	for( i = 0; i < q; i++ )
		for( j = 0; j < ng; j++ )
			soft_out[i][j] = soft[i][j];

	if( use_p_thr )
	{
		for( i = 0; i < ng; i++ )
			smask[i] = 0;
	}

 
	for( iter = 0; iter < maxiter; iter++ )
	{
		int dcmp = 0;
		int flag_cnt = 0;


#ifdef MAKE_DUMP
		FILE *fdump = fopen("d:\\huawei\\FHTdecoder\\C\\dump_c.txt", "wt");
		fprintf(fdump, "iter: %3d\n", iter);
#endif

		for( i = 0; i < ng; i++ )
		{
			double max;
			int pos;

			if( use_p_thr & smask[i] )	continue;

			max = 0.0;
			pos = 0;

			for( j = 0; j < q; j++ )
			{
				if( max < soft_out[j][i] )
				{
					max = soft_out[j][i];
					pos = j;
				}
			}

			qhard[i] = pos;
			if( use_p_thr && (1 - max < p_thr) )
			{
				int ih = i / m;
				int k  = i % m;

				smask[i] = 1;

				for( j = 0; j < q; j++ )
					soft[j][i] = 0.0; 
				
				soft[pos][i] = 1.0;
				/*
				// useless cycle
				for( j = 0; j < 2; j++ )
				{
					int row   = hb_ci[ih][j];
					int circ  = hb[row][ih];
					int ph    = hb_cj[ih][j];
					int pos_r = row*m;

					soft_in[pos][ ph*r + pos_r + (m - circ + k) % m ] = 1;
				}
				*/
				flag_cnt++;
			}
		}

		if( iter > 0 )
		{
#if 0
			FILE *fid = fopen("d:\\huawei\\FHTdecoder\\C\\hardq_c.txt", "wt");
			for( j = 0; j < ng; j++ )
				fprintf(fid, "%d\n", qhard[j]);
			fclose(fid);    
#endif

#if 0
			FILE *fid = fopen("d:\\huawei\\FHTdecoder\\C\\soft_out_c.txt", "wt");
			for( i = 0; i < q; i++ )
			{
				for( j = 0; j < ng; j++ )
					fprintf(fid, "%6.4f ", soft_out[i][j]);
				fprintf(fid, "\n");
			}
			fclose(fid);    
#endif
		}
#ifdef MAKE_DUMP
		{
			int iii, jjj;
			fprintf(fdump, "soft\n");
			for( iii = 0; iii < q; iii++ )
			{
				for( jjj = 0; jjj < ng; jjj++ )
					fprintf(fdump, "%8.6f ",  soft[iii][jjj] );
				fprintf(fdump, "\n");
			}
		}
#endif
		// check synfrome

		syn = syndrome_graph_cycle( hb, rh, nh, m, qhard, synd, buf, hc_rl, mul_tab/*, rweight*/ );

		if( syn == 0 )
		{
#ifdef MAKE_DUMP
			fclose( fdump );
#endif
			return iter; //break;	// Decoding is finished
		}
		
		// check nodes processing
		for( j = 0; j < rh; j++ )  //loop over checks
		{
			int pos_r = j * m;
			int weight = rweight[j];
			static short mask[1024] = {0};
			static short hard[1024] = {0};

			for( k = 0; k < m; k++ )
			{
				// decode constituent code of each row
				int ii;
				i = pos_r + k;
//				if( use_p_thr )
				{
					for( ii = 0; ii < weight; ii++ )
					{
						short col  = hb_ri[j][ii];
						short circ = hb[j][col];
						int jj     = (circ + k) % m;
						int pos_n  = col * m + jj;
						if( use_p_thr )
							mask[ii] = smask[pos_n];
						hard[ii] = qhard[pos_n];
					}
				}

#ifndef SOFT_IN_OUTS
				map_graph( soft_in, soft_outs, i, r, hc_rl[j], q_bits, weight, mul_tab, div_tab, mask, hard, HAD );
#else
				map_graph( soft_in, soft_in, i, r, hc_rl[j], q_bits, weight, mul_tab, div_tab, mask, hard, HAD );
#endif
#ifdef MAKE_DUMP
				{
					int iii, jjj;
					fprintf(fdump, "%4d soft_in\n", i);

					for( iii = 0; iii < q; iii++ )
					{
						for( jjj = 0; jjj < rweight[j]; jjj++ )
							fprintf(fdump, "%8.6f ",  soft_in[iii][jjj*r + i] );
						fprintf(fdump, "\n");
					}
/*
					fprintf(fdump, "mask\n", i);
					for( jjj = 0; jjj < weight; jjj++ )
						fprintf(fdump, "%1d ",  mask[jjj] );
					fprintf(fdump, "\n");
*/
/*
					fprintf(fdump, "hard\n", i);
					for( jjj = 0; jjj < weight; jjj++ )
						fprintf(fdump, "%2d ",  hard[jjj] );
					fprintf(fdump, "\n");
*/
					fprintf(fdump, "%4d soft_outs\n", i);
					for( iii = 0; iii < q; iii++ )
					{
						for( jjj = 0; jjj < rweight[j]; jjj++ )
							fprintf(fdump, "%8.6f ",  soft_outs[iii][jjj*r + i] );
						fprintf(fdump, "\n");
					}
				}
#endif
			}
		}

#ifdef MAKE_DUMP
		fprintf(fdump, "all soft_outs\n");

		for( i = 0; i < q; i++ )
		{
			for( j = 0; j < 4; j++ )
			{
				fprintf(fdump, "i = %4d, j = %4d\n", i+1, j+1 );
				for( k = 0; k < rg; k++ )
				{
					fprintf(fdump, "%8.6f ",  soft_outs[i][j*rg + k] );
				}
				fprintf(fdump, "\n");
			}
			fprintf(fdump, "\n");
		}
#endif

		// symbol nodes

		for( i = 0; i < rh; i++ )
			posh[i] = 0;



#ifdef COLUMN_BY_COLUMN


#ifdef MAKE_LIST
		for( i = 0; i < nh; i++ )
		{
			static double bf0[QMAX];
			static double bf1[QMAX];
			int pos_n = i * m;
			int row0 = hb_ci[i][0];
			int row1 = hb_ci[i][1];
			int circ0 = hb[row0][i];
			int circ1 = hb[row1][i];
			int ph0 = posh[row0];
			int ph1 = posh[row1];
			int pos_r0 = row0*m;
			int pos_r1 = row1*m;
			int index0 = ph0*r + pos_r0;
			int index1 = ph1*r + pos_r1;

			int kkk0 = (m - circ0) % m;
			int kkk1 = (m - circ1) % m;

			for( k = 0; k < m; k++ )
			{
				int ii = pos_n + k;
				double a;//, b, c;

				if( kkk0 == m )
					kkk0 = 0;

				if( kkk1 == m )
					kkk1 = 0;

				if( use_p_thr & smask[ii] )
				{
					for( j = 0; j < q; j++ )	soft_out[j][ii] = soft[j][ii];
				}
				else
				{
					static int flag[QMAX];
					static ELEMENT t[QMAX];
#ifdef SCALABLE
					static double sft[QMAX];
					static ui32 soft_log[QMAX];
					static ui32 bf0_log[QMAX];
					static ui32 bf1_log[QMAX];
					static ui32 ibf0[QMAX];
					static ui32 ibf1[QMAX];
#endif


#ifdef  SCALABLE
#ifndef SOFT_IN_OUTS
					for( j = 0; j < q; j++ )
					{
						bf0[j] = soft_outs[j][index0 + kkk0];
						bf1[j] = soft_outs[j][index1 + kkk1];
					}
#else
					for( j = 0; j < q; j++ )
					{
						bf0[j] = soft_in[j][index0 + kkk0];
						bf1[j] = soft_in[j][index1 + kkk1];
					}
#endif
					for( j = 0; j < q; j++ )
						sft[j] = soft[j][ii];

					//	rough calculation

					for( j = 0; j < q; j++ )
					{
						int xlog    = get_bitsize( (ui32)(sft[j]  * ONE_INP) );
						int b0log  = get_bitsize( (ui32)(bf0[j] * ONE_TMP));
						int b1log  = get_bitsize( (ui32)(bf1[j] * ONE_TMP));

						bf0_log[j]  = xlog + b1log;
						bf1_log[j]  = xlog + b0log;
						soft_log[j] = xlog + b0log + b1log;
					}

					// find best elements

					for( j = 0; j < q; j++ )
					{
						t[j].val = soft_log[j];
						t[j].pos = j;
					}

#else //SCALABLE
					for( j = 0; j < q; j++ )
					{
#ifndef SOFT_IN_OUTS
						double b0 = soft_outs[j][index0 + kkk0];
						double b1 = soft_outs[j][index1 + kkk1];
#else
						double b0 = soft_in[j][index0 + kkk0];
						double b1 = soft_in[j][index1 + kkk1];
#endif
						double x = soft[j][ii];
						double y0 = x * b0;
						double y1 = x * b1;

						bf0[j] = y1;
						bf1[j] = y0;
						soft_out[j][ii] = y1 * b0;
					}

					for( j = 0; j < q; j++ )
					{
						t[j].val = soft_out[j][ii];
						t[j].pos = j;
					}
#endif	// SCALABLE
					
					sort( t, q, NMAIN );
                   
					memset( flag, 0, sizeof(flag) );
					for( j = q - NMAIN; j < q; j++ )
						flag[t[j].pos] = 1;


#ifdef SCALABLE
					for( j = 0; j < q; j++ )
					{
						if( flag[j] )
						{
							// precision computation for strong elements
#if 0
							double b0 = bf0[j];
							double b1 = bf1[j];
							double x  = sft[j];
#else
							double b0 = (double)( (ui32)(bf0[j] * ONE_TMP) | 1 ) / (double)ONE_TMP;
							double b1 = (double)( (ui32)(bf1[j] * ONE_TMP) | 1) / (double)ONE_TMP;
							double x  = (double)( (ui32)(sft[j] * ONE_INP) | 1) / (double)ONE_INP;
#endif
							double y0 = x * b0;
							double y1 = x * b1;

							bf0[j] = y1;
							bf1[j] = y0;
							
							soft_out[j][ii] = y1 * b0;
						}
						else
						{
							// put rough results for weak elements
							//- 3 - EUG trick
							soft_out[j][ii] = log_to_double( soft_log[j] - 3, INP_FPP + TMP_FPP + TMP_FPP ); 

							bf0[j] = log_to_double( bf0_log[j] - 2, INP_FPP+ TMP_FPP );  

							bf1[j] = log_to_double( bf1_log[j] - 2, INP_FPP+ TMP_FPP );  
						}
					}

					// normalization
					a = 0.0;
					for( j = 0; j < q; j++ )
						a += soft_out[j][ii];

					a = 1.0 / a;
					for( j = 0; j < q; j++ )
						soft_out[j][ii] *= a;


					normalize_abs( bf0, q, 1, PROB_FPP );
					normalize_abs( bf1, q, 1, PROB_FPP );

					double_int_double( bf0, q, PROB_FPP );
					double_int_double( bf1, q, PROB_FPP );

					for( j = 0; j < q; j++ )
					{
						soft_in[j][index0 + kkk0] = bf0[j];
						soft_in[j][index1 + kkk1] = bf1[j];
					}

#else //SCALABLE
					// put average for weak elements
					b = c = 0;
					for( j = 0; j < q; j++ )
					{
						if( flag[j] == 0 )
						{
							b += bf0[j];
							c += bf1[j];
						}
					}

                    b /= q;
                    c /= q;
                    
					for( j = 0; j < q; j++ )
					{
						if( flag[j] == 0 )
						{
							bf0[j] = b;
							bf1[j] = c;
						}
					}

					// normalization
					a = 0.0;
					for( j = 0; j < q; j++ )
						a += soft_out[j][ii];

					a = 1.0 / a;
					for( j = 0; j < q; j++ )
						soft_out[j][ii] *= a;

					normalize( bf0, q, 1 );
					normalize( bf1, q, 1 );

					for( j = 0; j < q; j++ )
					{
						soft_in[j][index0 + kkk0] = bf0[j];
						soft_in[j][index1 + kkk1] = bf1[j];
					}
#endif  //SCALABLE
				}

				kkk0++;
				kkk1++;

				dcmp++;
			}

			posh[row0] += 1;
			posh[row1] += 1;

		}


#else //MAKE_LIST
		for( i = 0; i < nh; i++ )
		{
			static double bf0[1024];
			static double bf1[1024];
			int pos_n = i * m;
			int row0 = hb_ci[i][0];
			int row1 = hb_ci[i][1];
			int circ0 = hb[row0][i];
			int circ1 = hb[row1][i];
			int ph0 = posh[row0];
			int ph1 = posh[row1];
			int pos_r0 = row0*m;
			int pos_r1 = row1*m;
			int index0 = ph0*r + pos_r0;
			int index1 = ph1*r + pos_r1;

			int kkk0 = (m - circ0) % m;
			int kkk1 = (m - circ1) % m;

			for( k = 0; k < m; k++ )
			{
				int ii = pos_n + k;
				double a, b, c;

#ifdef MAKE_DUMP
//				fprintf(fdump, "%4d mask = %1d\n", dcmp, smask[ii] );				
				fprintf(fdump, "soft_out %8d\n", ii+1);
#endif

				if( kkk0 == m )
					kkk0 = 0;

				if( kkk1 == m )
					kkk1 = 0;

				if( use_p_thr & smask[ii] )
				{
					for( j = 0; j < q; j++ )	soft_out[j][ii] = soft[j][ii];
#ifdef MAKE_DUMP
					for( j = 0; j < q; j++ )
						fprintf(fdump, "%8.6f ", soft_out[j][ii] );
					fprintf(fdump, "\n");
#endif
				}
				else
				{
					for( j = 0; j < q; j++ )
					{
#ifndef SOFT_IN_OUTS
						double b0 = soft_outs[j][index0 + kkk0];
						double b1 = soft_outs[j][index1 + kkk1];
#else
						double b0 = soft_in[j][index0 + kkk0];
						double b1 = soft_in[j][index1 + kkk1];
#endif
						double x = soft[j][ii];
						double y0 = x * b0;
						double y1 = x * b1;
							
						bf0[j] = y1;
						bf1[j] = y0;
						soft_out[j][ii] = y1 * b0;
					}

					a = 0;

					for( j = 0; j < q; j++ )
						a += soft_out[j][ii];
#ifdef MAKE_DUMP
					for( j = 0; j < q; j++ )
						fprintf(fdump, "%8.6f ", soft_out[j][ii] );
					fprintf(fdump, "\n");
#endif
					a = 1.0 / a;
					for( j = 0; j < q; j++ )
						soft_out[j][ii] *= a;

					c = b = 0;

					for( j = 0; j < q; j++ )
					{
						b += bf0[j];
						c += bf1[j];
					}
#ifdef MAKE_DUMP
					fprintf(fdump, "soft_in\n");
					for( j = 0; j < q; j++ )
						fprintf(fdump, "%8.6f ", bf0[j] );
					fprintf(fdump, "\n");

					fprintf(fdump, "soft_in\n");
					for( j = 0; j < q; j++ )
						fprintf(fdump, "%8.6f ", bf1[j] );
					fprintf(fdump, "\n");
#endif
					b = 1.0 / b;
					c = 1.0 / c;
					for( j = 0; j < q; j++ )
					{
						bf0[j] *= b;
						bf1[j] *= c;
					}

					for( j = 0; j < q; j++ )
					{
						soft_in[j][index0 + kkk0] = bf0[j];
						soft_in[j][index1 + kkk1] = bf1[j];
					}
				}

				kkk0++;
				kkk1++;

				dcmp++;
			}

			posh[row0] += 1;
			posh[row1] += 1;

		}
#endif //MAKE_LIST

#else   //COLUMN_BY_COLUMN
		for( i = 0; i < nh; i++ )
		{
			int pos_n = i * m;
			int r0 = hb_ci[i][0];
			int r1 = hb_ci[i][1];
			int circ0 = hb[r0][i];
			int circ1 = hb[r1][i];
			int ph0 = posh[r0];
			int ph1 = posh[r1];
			int pos_r0 = r0*m;
			int pos_r1 = r1*m;
			int index0 = ph0*r + pos_r0;
			int index1 = ph1*r + pos_r1;

			for( j = 0; j < q; j++ )
			{
#ifndef SOFT_IN_OUTS
				rotate( &soft_outs[j][index0], &buf0[j][0], m - circ0, sizeof(soft_in[0][0]), m );
				rotate( &soft_outs[j][index1], &buf1[j][0], m - circ1, sizeof(soft_in[0][0]), m );
#else
				rotate( &soft_in[j][index0], &buf0[j][0], m - circ0, sizeof(soft_in[0][0]), m );
				rotate( &soft_in[j][index1], &buf1[j][0], m - circ1, sizeof(soft_in[0][0]), m );
#endif
			}

			for( j = 0; j < q; j++ )
			{
				for( k = 0; k < m; k++ )
				{
					int ii = pos_n + k;
					double x = soft[j][ii];
					double b0 = buf0[j][k];
					double b1 = buf1[j][k];
					double y0 = x * b0;
					double y1 = x * b1;

					buf0[j][k] = y1;
					buf1[j][k] = y0;
					soft_out[j][ii] = y1 * b0;
				}
			}


			for( k = 0; k < m; k++ )
			{
				int ii = pos_n + k;
				double a = 0;

				for( j = 0; j < q; j++ )
					a += soft_out[j][ii];

				a = 1.0 / a;
				for( j = 0; j < q; j++ )
					soft_out[j][ii] *= a;
			}


			for( k = 0; k < m; k++ )
			{
				double a = 0;
				double b = 0;
				int ii = pos_n + k;
			
				for( j = 0; j < q; j++ )
				{
					a += buf0[j][k];
					b += buf1[j][k];
				}

				a = 1.0 / a;
				b = 1.0 / b;
				for( j = 0; j < q; j++ )
				{
					buf0[j][k] *= a;
					buf1[j][k] *= b;
				}
			}

			for( j = 0; j < q; j++ )
			{
				rotate( &buf0[j][0], &soft_in[j][index0], circ0, sizeof(soft_in[0][0]), m );
				rotate( &buf1[j][0], &soft_in[j][index1], circ1, sizeof(soft_in[0][0]), m );
			}

			posh[r0] += 1;
			posh[r1] += 1;

		}
#endif    //  COLUMN_BY_COLUMN
#ifdef MAKE_DUMP
		fclose( fdump );
#endif
	}

	return -iter;	// errors detected but not corrected
};


#ifndef SKIP_MEX
// %
// % decoder interface:
// % 1. SETUP ( before simulation )
// % ========
// % result = decoder( decoder_type, c, HB, HC, M );
// %    decodec_type - defined decoder:
// %        0 - bp
// %        1 - sum prod
// %		2 - advanced sum prod
// %		3 - min-sum
// %		4 - int min-sum
// %		5 - int advanced sum prod
// %		6 - fht sum prod
// %        7 - int fht sum prod
// %    c  - size of symbol (in bits)
// %    HB - circulant matrix
// %    HC - symbol matrix
// %    M  - matrix extension

//
// %  result == 1, if OK
// %
// % 2. DECODE
// % =========
// %  [it, soft, decword] = decoder( chan_out, maxiter );
// %	 chan_out -  channel output
// %	 maxiter  -  max number of iterations
// %	 
// %     it       -  number of made iterations 
// %     soft     -  decoded codeword ( soft decision )
// %     decword  -  decoded codeword ( hard decision )
// %
// %
// % 3. CLEAR ( free memory allocations )
// % ========
// %  decoder();     
// % 


void mexFunction(int nOut, mxArray *pOut[], int nInp, const mxArray *pInp[])
{
    static int nh, mh, M;
    static int maxiter; // = mxGetPr(prhs[3])[0];
    static int N, R;
    static int dectype;
    static int BBsize;
    int i, iter;
	static double	ms_alpha	= MS_ALPHA;
	static double	ms_thr		= MS_THR;
	static int		ms_qbits	= MS_QBITS;
	static int		ms_dbits	= MS_DBITS;
	static int		decision	= 0;
    static int      q_bits;
    static int     q;
    static int     ng;
	static double p_thr;
	/*
	#define MS_ALPHA 0.9//1.0
	#define MS_THR   1.4	
	#define MS_QBITS 5
	#define MS_DBITS 16
	*/
    
    // этот вариант переделан только для q-ичного кода
    
    static DEC_STATE* state = NULL;

     switch( nInp )
     {
     case 0:
           // close decoder
           // mexPrintf("close the decoder\n");
            if( state ) 
                decod_close( state );
			state = NULL;
            return;
             
     case 2:
     case 3:
     {
         double *p;
         double **pp;
           // run decoder
           // mexPrintf("run the decoder\n");
          // [iter_C, softq_C, hardq_C] =  decode( soft, maxsteps );            
            maxiter = (int)mxGetPr(pInp[1])[0];
            q   = (int)mxGetM(pInp[0]);
            ng = (int)mxGetN(pInp[0]);
           
			p_thr = nInp == 2 ? 0.0 : (double)mxGetPr(pInp[2])[0];
		   	
			// mexPrintf("nOut: %d\n", nOut);
            // mexPrintf("maxiter: %d, q: %d, state_q: %d, n: %d, thr: %f\n", maxiter, q, state->q, ng, p_thr);
            

           if( state->q == q )
            {
                pOut[0]=mxCreateDoubleMatrix(q,ng,mxREAL);
                pOut[1]=mxCreateDoubleMatrix(1,1,mxREAL); 
           
	        if( nOut == 3 )
                     pOut[2]=mxCreateDoubleMatrix(1,ng,mxREAL);

                if( q == 1 )  
		    unpackRow(mxGetPr(pInp[0]), 1, ng, state->y );
                else
                    unpackMatrix(mxGetPr(pInp[0]), q, ng, state->qy );


                decision = 0;   
                
                if( q == 1 )
                {
		    if( nInp == 3 )
                        decision = (int)mxGetPr(pInp[2])[0];
		    else
                        decision = 0;
                }

                switch( dectype )
                {
                   case 0: iter = bp_decod_qc_lm( state, state->y, state->decword, maxiter, decision);			    break;
                   case 1: iter = sum_prod_decod_qc_lm( state, state->y, state->decword, maxiter, decision);	            break;
                   case 2: iter = sum_prod_gf2_decod_qc_lm( state, state->y, state->decword, maxiter, decision);            break;
                   case 3: iter = min_sum_decod_qc_lm( state, state->y, state->decword, maxiter, decision, ms_alpha );      break;
                   case 4: iter = imin_sum_decod_qc_lm( state, state->y, state->decword, maxiter, decision, ms_alpha, ms_thr, ms_qbits, ms_dbits); break;
                   case 5: iter = isum_prod_gf2_decod_qc_lm( state, state->y, state->decword, maxiter, decision);           break;
	           case 6: iter = sum_prod_gfq_decod_lm( state, state->qy, state->qhard, state->qdecword, maxiter, p_thr ); break;
                   case 7: iter = tdmp_sum_prod_gf2_decod_qc_lm( state, state->y, state->decword, maxiter, decision);       break;
                }



   
		pOut[1] =  mxCreateDoubleScalar(iter);

               

                p = mxGetPr(pOut[0]);
               
                if( q == 1 )
                    for(  i = 0;  i < ng;  i++ )
                        p[i] = state->decword[i];
                else
                    for(  i = 0;  i < ng;  i++ )
                        p[i] = state->qhard[i];
         
            }
            else
            {
                mexErrMsgTxt("Allocation error");
                *(mxGetPr(pOut[0])) = maxiter;
            }
            return;
     }
     default:
         // open decoder
         //mexPrintf("open the decoder\n");

         dectype = (int)mxGetPr(pInp[0])[0];
         switch( dectype )
         {
             case FHT_DEC:
                  // Open qDEc  decoder_type, c, HB, HC, M
              q_bits = (int)mxGetPr(pInp[1])[0];
              mh = (int)mxGetM(pInp[2]);
              nh = (int)mxGetN(pInp[2]);
                        
              M = (int)mxGetPr(pInp[4])[0];
              N = nh * M;
              R = mh * M;

              if( state )
                  decod_close( state );
              
              state = decod_open( dectype, q_bits, mh, nh, M );
              if( state == NULL )
              {
                 mexErrMsgTxt("Allocation error");
                 return;
              }
              
              unpackMatrix_double2short(mxGetPr(pInp[2]), mh, nh, state->hb );
              unpackMatrix_double2short(mxGetPr(pInp[3]), mh, nh, state->hc );
            
              decod_init( state );
              pOut[0]=mxCreateDoubleMatrix(1,1,mxREAL); 
              *(mxGetPr(pOut[0])) = 1;
 //             mexPrintf("decoder open: OK\n");
              return;

	     default:
                  // Open qDEc  decoder_type, c, HB, M
	      //mexPrintf("open bin decoder\n");

              q_bits = (int)mxGetPr(pInp[1])[0];
              mh = (int)mxGetM(pInp[2]);
              nh = (int)mxGetN(pInp[2]);
                        
              M = (int)mxGetPr(pInp[3])[0];
              N = nh * M;
              R = mh * M;

	      //mexPrintf("qbits = %d, mh = %d, nh = %d, M = %d\n", q_bits, mh, nh, M);


              if( state )
                  decod_close( state );
              
              state = decod_open( dectype, q_bits, mh, nh, M );
              if( state == NULL )
              {
                 mexErrMsgTxt("Allocation error");
                 return;
              }

	     // mexPrintf("Allocation OK\n");
              
              unpackMatrix_double2short(mxGetPr(pInp[2]), mh, nh, state->hd );

	      //mexPrintf("unpacking OK\n");
            
              decod_init( state );

	      //mexPrintf("init OK\n");

              pOut[0]=mxCreateDoubleMatrix(1,1,mxREAL); 
              *(mxGetPr(pOut[0])) = 1;
              mexPrintf("decoder open: OK\n");
              return;

         }
     }
     
#if 0
    else
    //if( nInp >= 4 )
    {
        // SETUP :result = bp_decod_qc_lm( HC, M, niter, dectype, decision ); OR result = bp_decod_qc_lm( HC, M, niter, dectype, decision, alpha );( for dectype == 3 )
        if( nOut != 1 )
        {
            mexErrMsgTxt("Only one output argument allowed for SETUP");
        }

        M = (int)mxGetPr(pInp[1])[0];
        maxiter = (int)mxGetPr(pInp[2])[0];
        dectype = (int)mxGetPr(pInp[3])[0];
		if( nInp >= 5 )
			ms_alpha = (double)mxGetPr(pInp[4])[0];
//        if( dectype < 0 || dectype > 1 )
//            mexErrMsgTxt("Illegal decoder type");
		if( nInp >= 6 )
			ms_thr = (double)mxGetPr(pInp[5])[0];
		if( nInp >= 7 )
			ms_qbits = (int)mxGetPr(pInp[6])[0];
		if( nInp == 8 )
			ms_dbits = (int)mxGetPr(pInp[7])[0];

        
        mh = (int)mxGetM(pInp[0]);
        nh = (int)mxGetN(pInp[0]);
        N = nh * M;
        R = mh * M;
        BBsize = __max(mh,nh-mh);
        mexPrintf("nh = %d mh = %d N = %d maxiter = %d\n", nh, mh, N, maxiter);
        
        // Check if dynamic arrays were allocated
        if( state )
        {
             decod_close( state, state->rh, state->nh, state->m, state->n, __max(state->rh,state->nh-state->rh), state->codec_id );
        }
        state = decod_open( dectype, mh, nh, N, M, BBsize, R );
        if( state == NULL )
             mexErrMsgTxt("Allocation error");
        unpackMatrix_double2short(mxGetPr(pInp[0]), mh, nh, state->hd );
        
        pOut[0]=mxCreateDoubleMatrix(1,1,mxREAL); 
        *(mxGetPr(pOut[0])) = 1;
        return;
    }
    
    // decoding
    if( nInp == 1 || nInp == 2 )
    {
        double *p;
        if( state == NULL )
            mexErrMsgTxt("Decoder not opend"); 
        
        if( nOut != 2 )
            mexErrMsgTxt("Only two output argument allowed for DECODING");
        // Decoding, only 1 input param - input vector
//        unpackRow(mxGetPr(pInp[0]), 1, N, y );
        pOut[0]=mxCreateDoubleMatrix(1,N,mxREAL);
        p = mxGetPr(pOut[0]);

		unpackRow(mxGetPr(pInp[0]), 1, N, state->y );
		if( nInp == 2 )
			decision = (int)mxGetPr(pInp[1])[0];
		else
			decision = 0;

        switch( dectype )
        {
            case 0: iter = bp_decod_qc_lm( state, state->y, state->decword, maxiter, decision);			 break;
            case 1: iter = sum_prod_decod_qc_lm( state, state->y, state->decword, maxiter, decision);	 break;
            case 2: iter = sum_prod_gf2_decod_qc_lm( state, state->y, state->decword, maxiter, decision); break;
			case 3: iter = min_sum_decod_qc_lm( state, state->y, state->decword, maxiter, decision, ms_alpha ); break;
			case 4: iter = imin_sum_decod_qc_lm( state, state->y, state->decword, maxiter, decision, ms_alpha, ms_thr, ms_qbits, ms_dbits); break;
            case 5: iter = isum_prod_gf2_decod_qc_lm( state, state->y, state->decword, maxiter, decision); break;

        }

		for( i = 0; i < N; i++ )
			p[i] = state->decword[i];

		pOut[1] =  mxCreateDoubleScalar(iter);
        //pOut[0]=mxCreateDoubleMatrix(1,1,mxREAL); 
        //*(mxGetPr(pOut[0])) = iter;
        return;
    }
    
    if( nInp == 0 && state != NULL )
    {
        // Check if dynamic arrays were allocated
        decod_close( state, state->rh, state->nh, state->m, state->n, __max(state->rh,state->nh-state->rh), state->codec_id );
        state = NULL;
    }
#endif  //0
}
#endif


