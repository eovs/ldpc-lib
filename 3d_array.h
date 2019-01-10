#ifndef _3D_ARRAY_H_
#define _3D_ARRAY_H_

typedef struct
{
	int **addr;
	int ndim;
	int size[10];
} ARRAY;

#define get_ncol( x )   (x.size[0])
#define get_nrow( x )   (x.size[1])
#define get_nsheet( x ) (x.size[2])
#define get_addr( x )   ((void*)x.addr)

void put_addr( ARRAY *x, void* addr );
void put_ncol( ARRAY *x, int n );
void put_nrow(  ARRAY *x, int n );
void put_nsheet(  ARRAY *x, int n );
void put_ndim( ARRAY *x, int n );

void* Alloc1d_int( int b );
void** Alloc2d_int( int b, int c );
void*** Alloc3d_int( int b, int c, int d );

ARRAY zeros( int ndim, ... );
ARRAY ones( int ndim, int n, int k );
ARRAY minus_ones( int ndim, int n, int k );
void free_ARRAY( ARRAY x );

int save_ARRAY( char *fileName, ARRAY x );
int find( int **B, int pos, int n, int mode, int val, int *res );

#endif //_3D_ARRAY_H_