#include <stdio.h>
#include <malloc.h>
#include <stdarg.h>


#include "3d_array.h"

void put_addr( ARRAY *x, void* addr ){ x->addr    = (int**)addr; }
void put_ncol( ARRAY *x, int n )     { x->size[0] = n;    }
void put_nrow(  ARRAY *x, int n )    { x->size[1] = n;    }
void put_nsheet(  ARRAY *x, int n )  { x->size[2] = n;    }
void put_ndim( ARRAY *x, int n )     { x->ndim    = n;    }

void* Alloc1d_int( int b )
{
	return (int*)malloc( b * sizeof(int) );
}



void** Alloc2d_int( int b, int c )
{
	void **p;
	int i;

	p = (void**)malloc( b * sizeof(void*) );
	for( i = 0; i < b; i++ )
		p[i] = Alloc1d_int( c );

	return p;
}

void*** Alloc3d_int( int b, int c, int d )
{
	void ***p;
	int i;

	p = (void***)malloc( d * sizeof(void**) );
	for( i = 0; i < d; i++ )
		p[i] = Alloc2d_int( b, c );

	return p;
}



ARRAY zeros( int ndim, ... ) 
{
	ARRAY x;
	int *addr1;
	int **addr2;
	int ***addr3;
	int i, j, k;
	int ncol, nrow, nsheet;
	va_list ptr;
	
	va_start( ptr, ndim );

	switch( ndim )
	{
	case 1:
		ncol = va_arg( ptr, int );
		put_addr( &x, Alloc1d_int( ncol ) );
		put_ncol( &x, ncol );
		addr1 = (int*)get_addr( x );
		for( i = 0; i < ncol; i++ )
				addr1[i] = 0;
		break;
	case 2:
		ncol = va_arg( ptr, int );
		nrow = va_arg( ptr, int );

		put_addr( &x, Alloc2d_int( nrow, ncol ) );
		put_ncol( &x, ncol );
		put_nrow( &x, nrow );
		put_ndim( &x, ndim );

		addr2 = (int**)get_addr( x );
		for( i = 0; i < nrow; i++ )
			for( j = 0; j < ncol; j++ )
				addr2[i][j] = 0;
		break;

	case 3:
		ncol   = va_arg( ptr, int );
		nrow   = va_arg( ptr, int );
		nsheet = va_arg( ptr, int );
		
		put_addr( &x, Alloc3d_int( nrow, ncol, nsheet ) );
		put_ncol( &x, ncol );
		put_nrow( &x, nrow );
		put_nsheet( &x, nsheet );
		put_ndim( &x, ndim );

		addr3 = (int***)get_addr( x );
		for( i = 0; i < nsheet; i++ )
			for( j = 0; j < nrow; j++ )
				for( k = 0; k < ncol; k++ )
					addr3[i][j][k] = 0;
		break;

	default:
		printf("zeros: unknown dimension\n" );
	}

	va_end( ptr );

	return x;
}

ARRAY ones( int ndim, int n, int k )
{
	ARRAY x;
	int **addr;
	int i, j;

	switch( ndim )
	{
	case 2:
		put_addr( &x, Alloc2d_int( n, k ) );
		put_ncol( &x, k );
		put_nrow( &x, n );
		put_ndim( &x, ndim );

		addr = (int**)get_addr( x );
		for( i = 0; i < n; i++ )
			for( j = 0; j < k; j++ )
				addr[i][j] = 1;
		break;

	default:
		printf("ones: unknown dimension\n" );
	}
	return x;
}

ARRAY minus_ones( int ndim, int n, int k )
{
	ARRAY x;
	int **addr;
	int i, j;

	switch( ndim )
	{
	case 2:
		put_addr( &x, Alloc2d_int( k, n ) );
		put_ncol( &x, n );
		put_nrow( &x, k );
		put_ndim( &x, ndim );

		addr = (int**)get_addr( x );
		for( i = 0; i < k; i++ )
			for( j = 0; j < n; j++ )
				addr[i][j] = -1;
		break;

	default:
		printf("minus_ones: unknown dimension\n" );
	}
	return x;
}


void free_ARRAY( ARRAY x )
{
	switch( x.ndim )
	{
	case 1:
		{
			int ncol = get_ncol( x );
			int *matr  = (int*)get_addr( x );
			free( matr );
			break;
		}

	case 2:
		{
			int nrow = get_nrow( x );
			int **matr  = (int**)get_addr( x );
			int i;
			for( i = 0; i < nrow; i++ )
				free( matr[i] );
			free( matr );
			break;
		}

	case 3:
		{
			int nrow   = get_nrow( x );
			int nsheet = get_nsheet( x );
			int ***y   = (int***)get_addr( x );
			int i, k;
			for( k = 0; k < nsheet; k++ )
			{
				int **matr = y[k];
				for( i = 0; i < nrow; i++ )
					free( matr[i] );
				free( matr );
			}
			free( y );
			break;
		}
	default: printf("free_ARRAY: unknown dimension\n" );
	}
}

int save_ARRAY( char *fileName, ARRAY x )
{
	int i, j;
	FILE *fp = fopen(fileName, "wt");
	int ncol   = get_ncol( x );
	int nrow   = get_nrow( x );
	int **y = (int**)get_addr( x );

	if( fp == NULL ){ printf("file %s: BAD\n", fileName ); return -1;}

	for( i = 0; i < nrow; i++ )
	{
		for( j = 0; j < ncol; j++ )
			fprintf( fp, "%4d", y[i][j] );
		fprintf( fp, "\n" );
	}

	fclose( fp );
	return 0;
}

int find( int **B, int pos, int n, int mode, int val, int *res )
{
	int i;
	int k;

	k = 0;
	if( mode == 1 )
	{
		// find in column
		k = 0;
		for( i = 0; i < n; i++ )
		{
			if( B[i][pos] >= val )
				res[k++] = i;
		}
	}
	else
	{
		// find in row
		k = 0;
		for( i = 0; i < n; i++ )
		{
			if( B[pos][i] >= val )
				res[k++] = i;
		}
	}

	return k;
}
