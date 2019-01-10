#include <stdio.h>
#include <malloc.h>

#include "trace_pm.h"
#include "3d_array.h"


#define my_min( x, y ) ((x) < (y) ? (x) : (y))
#define my_max( x, y ) ((x) < (y) ? (y) : (x))


void sum_sp_pol_mon
( 
//int wx, int *dx, int *cx, int *xs, 
//int wy, int *dy, int *cy, int *ys,
//int *wz, int *dz, int *cz, int *zs 
	POLYNOM *acc,
	POLYNOM *tmp,
	POLYNOM *res
);

void copy_2D_to_3d( ARRAY dst, ARRAY src )
{
	int i, j;
	int ***y = (int***)get_addr( dst );
	int **x  = (int**)get_addr( src );
	int nrow_src = get_nrow( src );
	int ncol_src = get_ncol( src );
	int nrow_dst = get_nrow( dst );
	int ncol_dst = get_ncol( dst );

	if( nrow_src != nrow_dst || ncol_src != ncol_dst )
		return;

	for( i = 0; i < nrow_src; i++ )
		for( j = 0; j < ncol_src; j++ )
			y[0][i][j] = x[i][j];
}

void copy_2D_to_3d_( ARRAY dst, ARRAY src )
{
	int i, j;
	int ***y = (int***)get_addr( dst );
	int **x  = (int**)get_addr( src );
	int nrow_src = get_nrow( src );
	int ncol_src = get_ncol( src );
	int nrow_dst = get_nrow( dst );
	int nsheet_dst = get_nsheet( dst );

	if( nrow_src != nsheet_dst || ncol_src != nrow_dst )
		return;

	for( i = 0; i < nrow_src; i++ )
		for( j = 0; j < ncol_src; j++ )
			y[i][j][0] = x[i][j];
}

void trace_bound_pol_mon_pm( ARRAY HD, int M, int gmax, int gtarget, int *S, int *SA )
{
	int i, j;
	int N;
	int d;
#ifdef NATURAL_POLYNOM
	ARRAY polB, aceB;
	int *buffer;
	int *pd;
	int *pc;
	int *ps;
	int *ps_org;
	int w;
#else
	ARRAY WB, DB, CB, BS;
#endif
	TANNER_MON_RES tmr;
	HP2A_MON_RES hmr;
	MUL_MAT_MAT_MON_RES mmmmr;
	ARRAY HP, HA;
	ARRAY A, AC;
	int cnt = 0;

	for( i = 0; i < gmax; i++ ) S[i] = 0;	// Spectrum
	for( i = 0; i < gmax; i++ ) SA[i] = 0;	// Spectrum
	
	tmr = tanner_mon( HD );
	HP = tmr.array1;
	HA = tmr.array2;


	hmr = hp2a_mon_pm( HP, HA, M );
	A  = hmr.array1;
	AC = hmr.array2;

	//save_ARRAY( "A.txt", A );
	//save_ARRAY( "AC.txt", AC );
	N = get_ncol( A );

#ifdef NATURAL_POLYNOM
	polB = zeros( 3, M, N, N );
	aceB = zeros( 3, M, N, N );
	buffer = (int*)malloc( M * 4 * sizeof( buffer[0] ) );
	pd = buffer;
	pc = pd + M;
	ps = pc + M;
	ps_org = ps + M;


	{
		int  **a = (int**)get_addr( A );
		int  **ac = (int**)get_addr( AC );

		int ***ppol = (int***)get_addr( polB );
		int ***pace = (int***)get_addr( aceB );

		for( i = 0; i < N; i++ )
		{
			for( j = 0; j < N; j++ )
			{
				if( a[i][j] >= 0 )
				{
					int pos = a[i][j];
					ppol[i][j][pos] = 1;
					pace[i][j][pos] = ac[i][j];
				}

			}
		}
	}
#else
	//	WB=double(A>=0);	//% B  is accumulator, WB are Hamming weights, 2D, integer 
	WB = zeros( 2, N, N );
	{
		int **wb = (int**)get_addr( WB );
		int **a  = (int**)get_addr( A );
		for( i = 0; i < N; i++ )
			for( j = 0; j < N; j++ )
				wb[i][j] = a[i][j] >= 0;
	}

	//	DB=zeros(N,N,2);    //% degrees,3D
	DB = zeros( 3, 2, N, N );
	//	DB(:,:,1)=A;        //% degrees, 3D
	copy_2D_to_3d_( DB, A );
	
	//	CB=zeros(N,N,2);    //% coefficients, 3D
	CB = zeros( 3, 2, N, N );
	
	//	CB(:,:,1)=WB;       //%  coefs, 3D
	copy_2D_to_3d_( CB, WB );

	//	BS=zeros(N,N,2); 
	BS = zeros( 3, 2, N, N );
	//	BS(:,:,1)=AC;       //% Accumullator for ACE, 3D
	copy_2D_to_3d_( BS, AC );
#endif

	for( d = 3; d < gmax; d += 2 )
	{
		//printf("d=%d\n", d);
		//[WB,DB,CB,BS]=mul_mat_mat_mon(WB,DB,CB,BS,A,AC,M);
#ifdef NATURAL_POLYNOM
		mmmmr = mul_mat_mat_mon( N, polB, aceB, A, AC, M );

		free_ARRAY( polB );
		free_ARRAY( aceB );
		polB = mmmmr.pol;
		aceB = mmmmr.ase;

#else
		mmmmr = mul_mat_mat_mon( WB, DB, CB, BS, A, AC, M );
		free_ARRAY( WB );
		free_ARRAY( DB );
		free_ARRAY( CB );
		free_ARRAY( BS );
		WB = mmmmr.ret1;
		DB = mmmmr.ret2;
		CB = mmmmr.ret3;
		BS = mmmmr.ret4;
#endif

		for( i = 0; i < N; i++ )
		{
#ifdef NATURAL_POLYNOM
			int ***ppol = (int***)get_addr( polB );
			int ***pace = (int***)get_addr( aceB );

			if( ppol[i][i][0] != 0  )
			{
				//% Increase spectra
				S[d] = S[d] + ppol[i][i][0];

				if( SA[d] > 0 )
					SA[d] = my_min( SA[d], pace[i][i][0] );
				else
					SA[d] = pace[i][i][0];

				//% Subtract constant
				
#if 1
				// know how!!!
				//%BS(I,I,1:wb)=BS(I,I,2:wb+1);
				w = pack_pol( ppol[i][i], pace[i][i], M, pd, pc, ps_org );
				ppol[i][i][0] = 0;
				w =  pack_pol( ppol[i][i], pace[i][i], M, pd, pc, ps );
				unpack_pol( w, pd, pc, ps_org, M, ppol[i][i], pace[i][i] );
#else
				ppol[i][i][0] = 0;
				//BS(I,I,1:wb)=BS(I,I,2:wb+1);
				pace[i][i][0] = 0;
#endif
			}
#else
			int ***pCB = (int***)get_addr( CB );
			int ***pBS = (int***)get_addr( BS );
			int ***pDB = (int***)get_addr( DB );
			int **pWB  = (int**)get_addr( WB );
			int wb = pWB[i][i];

			if( wb > 0 && pDB[i][i][0] == 0 )
			{
				int k;
				//% Increase spectra
				S[d] = S[d] + pCB[i][i][0];
				if( SA[d] > 0 )
					SA[d] = my_min( SA[d], pBS[i][i][0] );
				else
					SA[d] = pBS[i][i][0];
				
				//% Subtract constant
				wb = wb-1;
				pWB[i][i] = wb;
				//DB(I,I,1:wb)=DB(I,I,2:wb+1);
				for( k = 0; k < wb; k++ )
					pDB[i][i][k] = pDB[i][i][k+1];

				//CB(I,I,1:wb)=CB(I,I,2:wb+1);            
				for( k = 0; k < wb; k++ )
					pCB[i][i][k] = pCB[i][i][k+1];
				//%BS(I,I,1:wb)=BS(I,I,2:wb+1);
				//for( k = 0; k < wb; k++ )
				//	pBS[i][i][k] = pBS[i][i][k+1];
			}
#endif
		}

		if( S[d] > 0 )
		{
			cnt++;
			if( cnt == gtarget )
				break;
		}
	}

	for( i = 3; i < gmax; i += 2 )
		S[i] = S[i] / (i+1);
	
	for( i = 3; i < gmax; i += 2 )
		SA[i] /= 2;


#ifdef NATURAL_POLYNOM
	free_ARRAY( polB );
	free_ARRAY( aceB );
	free( buffer );
#else
	free_ARRAY( WB );
	free_ARRAY( BS );
	free_ARRAY( CB );
	free_ARRAY( DB );
#endif
	free_ARRAY( HP );
	free_ARRAY( HA );
	free_ARRAY( A );
	free_ARRAY( AC );
}

TANNER_MON_RES tanner_mon( ARRAY hd )
{
	int cw[1000];
	int i, j;
	int num;
	int IC;
	int **T, **W;
	TANNER_MON_RES res;
	int **HD = (int**)get_addr( hd );
	int nrow = get_nrow( hd );
	int ncol = get_ncol( hd );
	
	//cw - column weights
	num = 0;
	for( i = 0; i < ncol; i++ )
	{
		int s = 0;
		for( j = 0; j < nrow; j++ )
			s += HD[j][i] >= 0;
		cw[i] = s;
		num += cw[i]; 
	}

	res.array1 = minus_ones( 2, num, nrow+ncol );
	T = (int**)get_addr( res.array1 ); 
	
	res.array2 = zeros( 2, num, nrow+ncol );
	W = (int**)get_addr( res.array2 );


	IC = 0; // current column
	for( i = 0; i < ncol; i++ )
	{
		for( j = 0; j < nrow; j++ )
		{
			if( HD[j][i] >= 0 )
			{
				T[i][IC] = 1;
				T[ncol + j][IC] = HD[j][i];
				W[ncol + j][IC] = cw[i];

				IC=IC+1;
			}
		}
	}


	return res;
}


HP2A_MON_RES hp2a_mon_pm( ARRAY hd, ARRAY w, int M )
{
	HP2A_MON_RES res;
	int N    = get_ncol( hd );
	int nrow = get_nrow( hd );
	int **HD = (int**)get_addr( hd );
	int **W  = (int**)get_addr( w );
	int i, j, h;
	int I, J;
	int num;
	ARRAY a, wa, am, wam, ap, wap; 
	ARRAY e, x, s;
	int **A, **WA;
	int **Am, **WAm;
	int **Ap, **WAp;
	int **E, **X, **S;
	int **Apm, **WApm;
	int *f = (int*)malloc( N * sizeof(int) );
	
	a = minus_ones( 2, N*2, N*2 );	A = (int**)get_addr( a );
	wa = zeros( 2, N*2, N*2 );  	WA = (int**)get_addr( wa );
	
	// %minus		
	am = minus_ones( 2, N, N ); Am = (int**)get_addr( am );
	wam = zeros( 2, N, N );	WAm = (int**)get_addr( wam );

	// % plus
	ap = minus_ones( 2, N, N );	Ap = (int**)get_addr( ap );
	wap = zeros( 2, N, N );	WAp = (int**)get_addr( wap );


	//% Edges: loop over columns of Tanner
	e = zeros( 2, 2, 2*N );		E = (int**)get_addr( e );
	x = zeros( 2, 1, 2*N );		X = (int**)get_addr( x );
	//% column weights
	s = zeros( 2, 1, 2*N );		S = (int**)get_addr( s );
	


	for( i = 0; i < N; i++ )
	{
		for( j = 0, num = 0; j < nrow; j++ ){if( HD[j][i] >= 0 ){f[num++] = j;}}

		E[i][0] = f[0];   //% forward  part 1 to part 2
		E[i][1] = f[1]; 
		E[i+N][0] = f[1]; //% back  from 2 to 1;
		E[i+N][1] = f[0]; 
		X[i][0]   = HD[f[1]][i];
		X[i+N][0] = (M - X[i][0]) % M;
		S[i][0]   = W[f[1]][i];
		S[i+N][0] = W[f[1]][i];
	}

	//% Edge adjacency matrix (Aminus)
	for( i = 0; i < N; i++ )          //% i is head of the current edge
	{
		j = E[i][1];        //% j is tail of the current edge
		//% find all edges starting at j ending anywhere except i+N or i-N 
		//f=find( E(N+1:2*N,1)==j );
		for( h = N, num = 0; h < 2*N; h++ ){if( E[h][0] == j ){f[num++] = h;}}

		for( h = 0; h < num; h++ )
		{
			if( f[h] != i+N )
			{
				A[i][f[h]]     = X[f[h]][0];
				Am[i][f[h]-N]  = X[f[h]][0];
				WA[i][f[h]]    = S[f[h]][0];
				WAm[i][f[h]-N] = S[f[h]][0];
			}
		}
	}

	//%    (Aplus)
	for( i = N; i < 2*N; i++ ) //          % i is head of the current edge
	{
		j = E[i][1];        //% j is tail of the current edge
		//% find all edges starting at j ending anywhere except i+N or i-N 
		//% if i<=N, ib=i+N; else ib=i-N; end
		//f=find(E(1:N,1)==j); 
		for( h = 0, num = 0; h < N; h++ ){if( E[h][0] == j ){f[num++] = h;}}

		for( h = 0; h < num; h++ )
		{
			if( f[h] != i-N )
			{
				A[i][f[h]]     = X[f[h]][0];
				Ap[i-N][f[h]]  = X[f[h]][0];
				WA[i][f[h]]    = S[f[h]][0];
				WAp[i-N][f[h]] = S[f[h]][0];
			}
		}
	}


	//% Multiplication Am*Ap
	res.array1 = minus_ones( 2, N, N ); 	Apm  = (int**)get_addr( res.array1 );
	res.array2 = zeros( 2, N, N );			WApm = (int**)get_addr( res.array2 );
	
	for( I = 0; I < N; I++ )  //% rows Am
	{
		for( J = 0; J < N; J++ ) //% columns Ap
		{
			//%s=-1; ws=0;  //% accumulator
			for( j = 0; j < N; j++ ) // % along row
			{
				if( Am[I][j] >= 0 && Ap[j][J] >= 0 )
				{
					Apm[I][J]  =  (Am[I][j] +  Ap[j][J]) % M; 
					WApm[I][J] = WAm[I][j] + WAp[j][J];
				}
			}
		}
	}

	free_ARRAY( s );
	free_ARRAY( x );
	free_ARRAY( e );

	free_ARRAY( ap );
	free_ARRAY( wap );

	free_ARRAY( am );
	free_ARRAY( wam );

	free_ARRAY( a );
	free_ARRAY( wa );

	free( f );
	return res;
}


int find_first( int *x, int n, int lim )
{
	int i;
	for( i = 0; i < n; i++ )
		if( x[i] >= lim )
			break;
	return i;
}


void sum_sp_pol_mon
( 
	//int wx, int *dx, int *cx, int *xs, 
	//int wy, int *dy, int *cy, int *ys,
	//int *wz, int *dz, int *cz, int *zs 
	POLYNOM *acc, 
	POLYNOM *tmp,
	POLYNOM *res
)
{
	int cntX, cntY, cntZ;
	int wx  = acc->degree;
	int *dx = acc->pos;
	int *cx = acc->coef;
	int *xs = acc->ace; 

	int wy  = tmp->degree;
	int *dy = tmp->pos;
	int *cy = tmp->coef;
	int *ys = tmp->ace;

	int *wz = &res->degree;
	int *dz = res->pos;
	int *cz = res->coef;
	int *zs = res->ace; 
	
	*wz = 0;
	cntX = 0;
	cntY = 0;
	cntZ   = 0;

	while( 1 )
	{
		if( cntX == wx )
		{
			// dx - done
			while( cntY < wy )
			{
				dz[cntZ] = dy[cntY];
				cz[cntZ] = cy[cntY];
				zs[cntZ] = ys[cntY];
				cntY++;
				cntZ++;
			}
			*wz = cntZ;
			return; 
		}
		else
		{
			if( cntY == wy )
			{
				// dy - done;
				while( cntX < wx )
				{
					dz[cntZ] = dx[cntX];
					cz[cntZ] = cx[cntX];
					zs[cntZ] = xs[cntX];
					cntX++;
					cntZ++;
				}
				*wz = cntZ;
				return; 
			}
			else
			{
				if( dx[cntX] < dy[cntY] )
				{
					dz[cntZ] = dx[cntX];
					cz[cntZ] = cx[cntX];
					zs[cntZ] = xs[cntX];
					cntX++;
					cntZ++;
				}
				else
				{
					if( dx[cntX] > dy[cntY] )
					{
						dz[cntZ] = dy[cntY];
						cz[cntZ] = cy[cntY];
						zs[cntZ] = ys[cntY];
						cntY++;
						cntZ++;
					}
					else
					{
						int t = cx[cntX] + cy[cntY];
						if( t != 0 )
						{
							dz[cntZ] = dx[cntX];
							cz[cntZ] = t;
							zs[cntZ] = my_min( xs[cntX], ys[cntY] );
							cntX++;
							cntY++;
							cntZ++;
						}
					}
				}
			}

		}
	}
}

#ifdef NATURAL_POLYNOM
void unpack_pol( int wx, int dx[], int cx[], int xs[], int M, int pol[], int ase[] )
{
	int i;
	for( i = 0; i < M; i++ )	pol[i] = 0;
	for( i = 0; i < M; i++ )	ase[i] = 0;

	for( i = 0; i < wx; i++ )
	{
		pol[dx[i]] = cx[i];
		ase[dx[i]] = xs[i];
	}
}

int pack_pol( int pol[], int ase[], int M, int dx[], int cx[], int xs[] )
{
	int i, j;

	j = 0;
	for( i = 0; i < M; i++ )
	{
		if( pol[i] )
		{
			dx[j] = i;
			cx[j] = pol[i];
			xs[j] = ase[i];
			j++;
		}
	}

	return j;
}

void add_pol( int pol0[], int ace0[], int pol1[], int ace1[], int M )
{
	int i;
	
	for( i = 0; i < M; i++ )
	{
		if( pol0[i] == 0 )
		{
			pol0[i] = pol1[i];
			ace0[i] = ace1[i];
		}
		else
		{
			if( pol1[i] == 0 )
			{
				pol0[i] = pol0[i];
				ace0[i] = ace0[i];
			}
			else
			{
				pol0[i] += pol1[i];
				ace0[i] = my_min( ace0[i], ace1[i] );
			}
		}
	}
}

MUL_MAT_MAT_MON_RES mul_mat_mat_mon( int bb, ARRAY polB, ARRAY aceB, ARRAY A, ARRAY AS, int M)
{
	ARRAY polY, aceY;
	int I, J;
	MUL_MAT_MAT_MON_RES res;

	int **pA  = (int**)get_addr( A );
	int **pAS = (int**)get_addr( AS );
	int ***ppolB;
	int ***paceB;
	int ***ppolY;
	int ***paceY;

	int h, i, j;

	int *pol1 = (int*)malloc( M * sizeof(pol1[0]) );
	int *ace1 = (int*)malloc( M * sizeof(ace1[0]) );


	polY = zeros( 3, M, bb, bb );	
	aceY = zeros( 3, M, bb, bb );	


	//% main loop
	ppolB = (int***)get_addr( polB );
	paceB = (int***)get_addr( aceB );

	ppolY = (int***)get_addr( polY );
	paceY = (int***)get_addr( aceY );

	for( I=0; I < bb; I++ ) //% rows of B
	{

		for( J=0; J < bb; J++ ) //% columns of A (square marices)
		{

			int* pol0 = ppolY[I][J];
			int* ase0 = paceY[I][J];

			for( j = 0; j < bb; j++ ) //% loop along column
			{
				int a  =  pA[j][J];
				int as = pAS[j][J];

				if( a >= 0 )
				{
					int wb = 0;
					for( h = 0; h < M; h++ )
					{
						if( ppolB[I][j][h] )
						{
							wb = 1;
							break;
						}
					}

					if( wb )
					{
						for( h = 0, i = M - a; i < M; i++, h++ ) pol1[h] = ppolB[I][j][i];
						for( i = 0; i < M - a; i++, h++ )        pol1[h] = ppolB[I][j][i];

						for( h = 0, i = M - a; i < M; i++, h++ ) ace1[h] = paceB[I][j][i];
						for( i = 0; i < M - a; i++, h++ )        ace1[h] = paceB[I][j][i];

						for( h = 0; h < M; h++ )
						{
							if( pol1[h] )
								ace1[h] += as;
						}

						add_pol( pol0, ase0, pol1, ace1, M );
					}
				}
			}
		}
	}

	free( pol1 );
	free( ace1 );

	res.pol = polY;
	res.ase = aceY;
	return res;
}
#else
MUL_MAT_MAT_MON_RES mul_mat_mat_mon( ARRAY WB, ARRAY DB, ARRAY CB, ARRAY BS, ARRAY A, ARRAY AS, int M)
{
	int bb, LB, LX;
	ARRAY WX, CX, DX, XS;
	int I, J;
	int w;
	MUL_MAT_MAT_MON_RES res;

	int **pWB = (int**)get_addr( WB );
	int **pA = (int**)get_addr( A );
	int **pAS = (int**)get_addr( AS );
	int ***pDB = (int***)get_addr( DB );
	int ***pCB = (int***)get_addr( CB );
	int ***pBS = (int***)get_addr( BS );
	int **pWX;
	int ***pDX;
	int ***pCX;
	int ***pXS;
	int wx;
	int* db;
	int* cb;
	int* bs;
	int* dx;
	int* cx;
	int* xs;
	int* rdx = (int*)malloc(M*sizeof(int));
	int* rcx = (int*)malloc(M*sizeof(int));
	int* rxs = (int*)malloc(M*sizeof(int));
	int tpos;
	int h, i, j;

	//% A,B are degree matrices
	//	% B is 3D
	//	% A is 2D
	//	% M is modulo

	//	[bb,~,LB]=size(DB);
	bb = get_nrow( DB );
//	LB = get_nsheet( DB );
	LB = get_ncol( DB );

	LX=LB+1;
	//db=zeros(1,LB); % Column degrees 
	db = (int*)calloc( LB, sizeof(db[0]) );	
	//cb=zeros(1,LB); % coefs
	cb = (int*)calloc( LB, sizeof(cb[0]) );
	//bs=zeros(1,LB); % ace
	bs = (int*)calloc( LB, sizeof(bs[0]) );
 
	//% ARRAY for result
	//WX=zeros(bb,bb);    % weights
	WX = zeros( 2, bb, bb );	
	//DX=zeros(bb,bb,LX); % degrees
	DX = zeros( 3, M, bb, bb );
	//CX=zeros(bb,bb,LX); % coefs
	CX = zeros( 3, M, bb, bb );
	//XS=zeros(bb,bb,LX); % ACEs
	XS = zeros( 3, M, bb, bb );

	//% main loop
	w=0;  //% max weight

	pWX = (int**)get_addr( WX );
	pDX = (int***)get_addr( DX );
	pCX = (int***)get_addr( CX );
	pXS = (int***)get_addr( XS );

	for( I=0; I < bb; I++ ) //% rows of B
	{
		POLYNOM acc, tmp, res;

		for( J=0; J < bb; J++ ) //% columns of A (square marices)
		{
			dx = &pDX[I][J][0]; //% degrees
			cx = &pCX[I][J][0]; //% coefs
			xs = &pXS[I][J][0]; //% ACEs
			wx = 0; 

			acc.pos    = dx;
			acc.coef   = cx;
			acc.ace    = xs;

			tmp.pos    = db;
			tmp.coef   = cb;
			tmp.ace    = bs;

			res.pos    = rdx;
			res.coef   = rcx;
			res.ace    = rxs;

			for( j = 0; j < bb; j++ ) //% loop along column
			{
				int a = pA[j][J], as = pAS[j][J];
				int wb = pWB[I][j];

				if( wb > 0 && a >= 0 )
				{
					acc.degree = wx;
					tmp.degree = wb;

					tpos = find_first(  &pDB[I][j][0], wb, M-a );

					for( h = 0, i = tpos; i < wb; i++, h++ ) db[h] = pDB[I][j][i] + a - M;
					for( i = 0; i < tpos; i++, h++ )         db[h] = pDB[I][j][i] + a;

					for( h = 0, i = tpos; i < wb; i++, h++ ) cb[h] = pCB[I][j][i];
					for( i = 0; i < tpos; i++, h++ )         cb[h] = pCB[I][j][i];

					for( h = 0, i = tpos; i < wb; i++, h++ ) bs[h] = pBS[I][j][i] + as;
					for( i = 0; i < tpos; i++, h++ )         bs[h] = pBS[I][j][i] + as;
	

					//% Accumulation
					//[wx,dx,cx,xs]=sum_sp_pol_mon( wx,dx,cx,xs, wb,db,cb,bs );
					sum_sp_pol_mon(  &acc, &tmp, &res );
					wx = res.degree;
					for( h = 0; h < wx; h++ )	dx[h] = res.pos[h];
					for( h = 0; h < wx; h++ )	cx[h] = res.coef[h];
					for( h = 0; h < wx; h++ )	xs[h] = res.ace[h];
				}
			}

			pWX[I][J] = wx;    //% weights
			w = my_max( w, wx );  //% max weight
		}
	}

	if( w < my_max(LX,2) )
		LX = w;
	
	free( db );
	free( cb );
	free( bs );

	free( rdx );
	free( rcx );
	free( rxs );

	res.ret1 = WX;
	res.ret2 = DX;
	res.ret3 = CX;
	res.ret4 = XS;
	res.ret5 = LX;
	return res;
}
#endif
