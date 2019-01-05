#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "modulation.h"

#ifndef SKIP_MEX
#include <mex.h>
#endif  //SKIP_MEX
/*
 * Module for direct and inverse codeword permutation 
 */



#ifndef SKIP_MEX
void unpackMatrix(double mat[], int height, int width, double *result[] )
{
    int i,j;
	for( i = 0; i < height; ++i)
	{
		for( j = 0; j < width; ++j)
		{
			result[i][j] = mat[j * height + i];
		}
	}
}
void unpackMatrix_double2int(double mat[], int height, int width, int *result[] )
{
    int i,j;
	for( i = 0; i < height; ++i)
	{
		for( j = 0; j < width; ++j)
		{
			result[i][j] = (int)mat[j * height + i];
		}
	}
}


#endif  // SKIP_MEX
//unpackMatrix_double2short
static void unpackMatrix_double2short(double mat[], int height, int width, short *result[] )
{
    int i,j;
	for( i = 0; i < height; ++i)
	{
		for( j = 0; j < width; ++j)
		{
			result[i][j] = (short)mat[j * height + i];
		}
	}
}


static int** Alloc2d_int( int b, int c )
{
	int **p;
	int i;
	// MB: unused and gives a memory leak
	//int *buf = (int*)calloc( b*c, sizeof(int) );

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

static void free2d_int( int **p )
{
	free( p[0] );
	free( p );
}


static short** Alloc2d_short( int b, int c )
{
	short **p;
	int i;
	// MB: unused and gives a memory leak
	//int *buf = (int*)calloc( b*c, sizeof(int) );

	p = (short**)calloc( b, sizeof(short*) );
	assert(p);
	p[0] = (short*)calloc( b*c, sizeof(short) );
	for( i = 1; i < b; i++ )
	{
		p[i] = p[i-1] + c;
		assert( p[i] );
	}
	return p;
}

static void free2d_short( short **p )
{
	free( p[0] );
	free( p );
}

#define NOT_FOUND   -1
static int Qav[4] = { 4, 16, 64, 256 };
static int logQav[4] = { 2, 4, 6, 8 };
static int findQ( int q )
{
    for( int i = 0; i < 4; i++ )
    {
        if( q == Qav[i] )
            return i;
    }
    return NOT_FOUND;
                
}

int find_max( int a[], int size )
{
    int maxval = a[0];
    for( int i = 1; i < size; i++ )
        if( a[i] > maxval )
             maxval = a[i];
    return maxval;
}

static unsigned int rand_next = 1;
int myrand(void)
{
    rand_next = rand_next * 1103515245 + 12345;
    return (int)(rand_next & 0x3FFFFFFF );
}



PERMSTATE* Permutations_Open( int b, int c, int M, int QAM, int halfmlog, int perm_mode, int block_size, int step_size )
{
    PERMSTATE* st;
    //int n0, c0;
    //int i, j;
    //int n0;
    
    st = (PERMSTATE*)calloc(1, sizeof(PERMSTATE));
    if( !st ) 
        return NULL;
    
    st->b = b;
    st->c = c;
    st->M = M;
    st->QAM = QAM;
    st->halfmlog = halfmlog;
    st->mlog = halfmlog * 2;
    st->perm_mode = perm_mode;
    st->Ncode = c * M;
    st->Rcode = b * M;

    
    
    if( perm_mode == 3 )
    {
        st->block_size = block_size;
        st->Nblocks = st->Ncode / block_size;
		st->ShortBlockSize = st->Ncode % block_size;
    }
    if( perm_mode == 4 )
    {
        st->step_size = step_size;
        st->block_size = st->Ncode / step_size;
        st->Nblocks = st->Ncode / st->block_size;
		st->ShortBlockSize = st->Ncode % st->block_size;
    }
    rand_next = 1;  // initial value for internal random generator
    
	st->perm_short		= NULL;
	st->invperm_short	= NULL;
    
    st->HC = Alloc2d_int( b, c );
    if( st->HC==NULL)
        return NULL;
    

    int rem = st->Ncode % st->mlog;
    
    st->extra_bits = 0;
    if( rem )
    {
        st->extra_bits = st->mlog - rem;
    }
    
    st->n0 = st->Ncode / halfmlog;
    st->c0 = c / halfmlog;

    int guard = st->n0 ;

    st->cw = (int*)calloc(c, sizeof( st->cw[0] ));
    if( !st->cw )
        return NULL;

    st->V = (int*)calloc( c, sizeof( st->V[0] ));
    if( !st->V )
        return NULL;

    st->p = (int*)calloc( c, sizeof( st->p[0] ));
    if( !st->p )
        return NULL;

    st->q = (int*)calloc( c, sizeof( st->q[0] ));
    if( !st->q )
        return NULL;

    st->s = (int*)calloc( st->n0, sizeof( st->s[0] ));
    if( !st->s )
        return NULL;

    st->perm = (int*)calloc( st->Ncode+guard, sizeof( st->perm[0] ) );
    if( st->perm == NULL )
        return NULL;
    st->tmpbuf = (double*)calloc( st->Ncode+guard, sizeof( st->tmpbuf[0] ) );
    if( st->tmpbuf == NULL )
        return NULL;
    st->invperm = (int*)calloc( st->Ncode+guard, sizeof( st->invperm[0] ) );
    if( st->invperm == NULL )
        return NULL;
    
    st->temp = (PERM_FOR_ORDER*)calloc( st->Ncode+guard, sizeof(PERM_FOR_ORDER));
    if( !st->temp )
        return NULL;

	st->buffer = (double*)calloc( st->Ncode+guard, sizeof( st->buffer[0] ));
	if( !st->buffer )
		return NULL;

	st->perm_codeword = (double*)calloc( st->Ncode+guard, sizeof( st->perm_codeword[0] ));
	if( !st->perm_codeword )
		return NULL;

	if( st->ShortBlockSize )
	{
		st->perm_short		= (int*)calloc( st->ShortBlockSize, sizeof( st->perm_short[0] ));
		st->invperm_short	= (int*)calloc( st->ShortBlockSize, sizeof( st->invperm_short[0] ));
	}
    
    
    return st;

}

void Permutations_Close( PERMSTATE *st )
{
    if( st->HC ) free2d_int( st->HC );
    if( st->perm ) free( st->perm );
    if( st->invperm ) free( st->invperm );
    if( st->cw ) free( st->cw );
    if( st->V ) free( st->V );
    if( st->p ) free( st->p );
    if( st->q ) free( st->q );
    if( st->s ) free( st->s );
    if( st->temp ) free( st->temp );
    if( st->tmpbuf ) free( st->tmpbuf );
	if( st->buffer ) free( st->buffer );
	if( st->perm_codeword ) free( st->perm_codeword );
	if( st->perm_short ) free( st->perm_short );
	if( st->invperm_short ) free( st->invperm_short );
    if( st ) free( st );
    
}

int permcomp( const void* a, const void *b )
{
    PERM_FOR_ORDER *p1 = (PERM_FOR_ORDER*) a;
    PERM_FOR_ORDER *p2 = (PERM_FOR_ORDER*) b;
    
    if( p1->val == p2->val )
        return 0;
    else
        return p1->val - p2->val;
}

int check_rn( int rn, int a[], int size )
{
    int i;
    for( i = 0; i < size; i++ )
    {
        if( a[i] == rn )
            return 0;
    }
    return 1;
}

void random_perm_gen( int perm[], int size )
{
    //int i;
    int cnt = 0;
    
    //mexPrintf("Start generation size = %d\n", size);
    while( cnt < size )
    {
        int rn = myrand() % size;
        int ret = check_rn( rn, perm, cnt );
        if( ret == 0 )
            continue;
        perm[cnt++] = rn;
//        if( cnt % 1000 == 0 )
//            mexPrintf("cnt = %d\n", cnt );
    }
//    mexPrintf("Stop generation\n");

}
void Permutation_Init(  PERMSTATE* state, short **hc )
{
    int i, j, l;
    int c = state->c;
    int b = state->b;
    int mlog = state->mlog;
    int halfmlog = state->halfmlog;
    int M = state->M;
    
    for( i = 0; i < b; i++ )
    {
        for( j = 0; j < c; j++ )
        {
            state->HC[i][j] = hc[i][j];
        }
    }
    
    int n0 = state->Ncode / halfmlog;
    int c0 = c / halfmlog;
        
    int n0_mod = state->Ncode % halfmlog;
    int c0_mod = c % halfmlog;
    
    #ifndef SKIP_MEX
    mexPrintf("c0_mod = %d  n0_mod = %d\n", c0_mod, n0_mod );
    #endif
    memset( state->cw, 0 , state->c*sizeof( state->cw[0] ));

    for( j = 0; j < c; j++ )
    {
        for( i = 0; i < b; i++ )
        {
            if( state->HC[i][j] >= 0 )
                state->cw[j] += 1;
        }
    }

    int mw = find_max( state->cw, c );
    #ifndef SKIP_MEX
    mexPrintf("mw  = %d  halfmlog = %d\n", mw, halfmlog );
    for( j = 0; j < c; j++ )
        mexPrintf("%3d ", state->cw[j] );
    mexPrintf("\n");
    
    #endif
    
    int mrp = 0;
    int k = 0;

    for( j = 0; j < halfmlog; j++ )
    {
//        for( i = 0; i < c0 ; i++ )
        for( i = 0; i*halfmlog < c ; i++ )
        {
            state->V[j*c0 + i] = k + i*halfmlog;
        }
        k++;
    }

    #ifndef SKIP_MEX
    for( j = 0; j < c; j++ )
        mexPrintf("%3d ", state->V[j] );
    mexPrintf("\n");
    
    #endif
    
    if( c0_mod )
    {
        for( k = c0*halfmlog; k < c; k++ )
            state->V[k] = k;
    }

    for( i = 0; i < c; i++ )
    {
        if( state->cw[state->V[i]] == mw )
            mrp++;
    }

    
    j = l = 0;
    for( i = 0; i < c; i++ )
    {
        if( state->V[i] < mrp )
            state->p[j++] = i;
        else
            state->q[l++] = state->V[i];
    }

      
    j = 0;
    for( i = 0; i < c-mrp; i++ )
        state->V[j++] = state->q[i];    //state->V[state->q[i]];
    for( i = 0; i < mrp; i++ )
        state->V[c-mrp+i] = i;
#ifndef SKIP_MEX
    mexPrintf("mrp = %d\n",mrp);
    mexPrintf("V (q,p):\n");
    for( i = 0; i < c; i++ )
        mexPrintf("%d ",state->V[i]);
    mexPrintf("\n");
#endif
        
//h=0;
//s=zeros(1,n0);
    int h = 0;
    memset( state->s, 0, n0*sizeof(state->s[0]));
        
    #ifndef SKIP_MEX
    for( i = 0; i < state->Ncode; i++ )
        state->perm[i] = -1;        // for testing only
    #endif

    switch( state->perm_mode )
    {
    case 0:     // No permutation
            for( i = 0; i < state->Ncode; i++ )
                state->perm[i] = i;
            break;
            
    case 1:     // random perm.
                random_perm_gen( state->perm, state->Ncode );
#ifndef SKIP_MEX
                {
                    FILE *fp = fopen("rand_perm.txt","wt");
           
                    for( i = 0; i < state->Ncode; i+= state->halfmlog*2 )
                    {
                        int k;
                        for( k = 0; k < state->halfmlog; k++ )
                            fprintf(fp, "%5d ", state->perm[i+k] );
                        fprintf(fp, "  " );
                        for( k = state->halfmlog; k < state->mlog; k++ )
                            fprintf(fp, "%5d ", state->perm[i+k] );
                        fprintf(fp, "\n" );
                    }
                    fclose(fp);
                }       
#endif
// check distribution of 'bad'  positions
#ifndef SKIP_MEX
                {
                    int badfrq[9] = {0};
                    int sum = 0;
                    int badlim = mrp * state->M;
                    for( i = 0; i < state->Ncode; i+= state->halfmlog )
                    {
                        int ind = 0;
                        int k, mask = 1;;
                        for( k = 0; k < state->halfmlog; k++, mask <<= 1 )
                        {
                            if( state->perm[i+k] < badlim ) {ind |= mask; sum++; }
                        }
                        if( ind )
                            badfrq[ind-1]++;
                    }
                    mexPrintf("sum = %d\n", sum );
                    {
                        int k;
                        for( k = 0; k < state->mlog+1; k++ )
                            mexPrintf( "%5d ", badfrq[k] );
                        mexPrintf("\n" );
                        
                    }
                    //mexPrintf("frq:  %d %d %d %d %d %d %d\n", badfrq[0],  badfrq[1],  badfrq[2],  badfrq[3],  badfrq[4],  badfrq[5],  badfrq[6] ) ;        
                }    
#endif                
                break;

        case 2:     // deterministic perm.

            if( c0_mod == 0 && n0_mod == 0 )
            {  // simplest case
                for( i = 0; i < state->halfmlog; i++ )
                {
                    for( j = 0; j < state->c0; j++ )
                    {
                        for( k = 0; k < state->M; k++ )
                            state->s[j*state->M + k] = state->V[j+h]*state->M + k;
                    }
                    
                    for( j = k = 0; k < state->n0; k++, j += state->halfmlog )
                        state->perm[i+j] = state->s[k];
            
                    h += state->c0;
                    //mexPrintf(" h = %d\n", h );
                }
            }
            
            if( c0_mod != 0 )   //&& n0_mod == 0 )
            {
                int ibad = c - mrp;
                int imix1 = ibad - mrp;
                int imix2 = imix1 - mrp;
                int imix3 = imix2 - mrp;
                
                #ifndef SKIP_MEX
                mexPrintf("ibad = %d   imix1 = %d   imix2 = %d   imix3 = %d\n", ibad, imix1, imix2, imix3 ); 
                #endif
                int limit, len;
                switch( halfmlog )
                {
                    case 1:
                        limit = 0; // just to silent down uninitialized error messages
                        break;
                    case 2:
                        limit = imix1;
                        break;
                    case 3:
                        limit = imix2;
                        break;
                    case 4:
                        limit = imix3;
                        break;
                    default:
                        assert(!!"Unknown option for halfmlog!!");
                        exit(1);
                        break;
                }
                len = limit * M;

                for( k = 0; k < M; k++ )
                {
                    for( i = 0; i < limit; i++ )
                    {
                        state->perm[i*M + k] = state->V[i]*M + k;
                    }
                }
                
                
                int rem = len % halfmlog;
                int l = len;
                
                #ifndef SKIP_MEX
                mexPrintf("limit = %d  len = %d  rem = %d\n", limit,len,rem );
                #endif
                switch( halfmlog )
                {
                    case 1:
                        break;
                    case 2:
                        
                        for( i = 0; i < mrp; i++ )
                        {
                            switch( rem )
                            {
                            case 0:
                                for( k = 0; k < M; k++ )
                                {
                                    state->perm[l++] = state->V[imix1+i]*M + k;
                                    state->perm[l++] = state->V[ibad+i]*M + k;
                                }
                                break;
                            case 1:
                                for( k = 0; k < M; k++ )
                                {
                                    state->perm[l++] = state->V[ibad+i]*M + k;
                                    state->perm[l++] = state->V[imix1+i]*M + k;
                                }
                                break;
                            }
                        }
                        break;
                    case 3:
                        for( i = 0; i < mrp; i++ )
                        {
                            switch( rem )
                            {
                            case 0:
                                for( k = 0; k < M; k++ )
                                {
                                    state->perm[l++] = state->V[imix1+i]*M + k;
                                    state->perm[l++] = state->V[imix2+i]*M + k;
                                    state->perm[l++] = state->V[ibad+i]*M + k;
                                }
                                break;
                            case 1:
                                for( k = 0; k < M; k++ )
                                {
                                    state->perm[l++] = state->V[imix1+i]*M + k;
                                    state->perm[l++] = state->V[ibad+i]*M + k;
                                    state->perm[l++] = state->V[imix2+i]*M + k;
                                }
                                break;
                            case 2:
                                for( k = 0; k < M; k++ )
                                {
                                    state->perm[l++] = state->V[ibad+i]*M + k;
                                    state->perm[l++] = state->V[imix1+i]*M + k;
                                    state->perm[l++] = state->V[imix2+i]*M + k;
                                }
                                break;
                            }
                        }
                        break;
                    case 4:
                        //int imix3 = imix2 - M;
                        //int i = imix3 * M;
                        for( i = 0; i < mrp; i++ )
                        {
                            switch( rem )
                            {
                            case 0:
                                for( k = 0; k < M; k++ )
                                {
                                    state->perm[l++] = state->V[imix1+i]*M + k;
                                    state->perm[l++] = state->V[imix2+i]*M + k;
                                    state->perm[l++] = state->V[imix3+i]*M + k;
                                    state->perm[l++] = state->V[ibad+i]*M + k;
                                }
                                break;
                            case 1:
                                for( k = 0; k < M; k++ )
                                {
                                    state->perm[l++] = state->V[imix1+i]*M + k;
                                    state->perm[l++] = state->V[imix2+i]*M + k;
                                    state->perm[l++] = state->V[ibad+i]*M + k;
                                    state->perm[l++] = state->V[imix3+i]*M + k;
                                }
                                break;
                            case 2:
                                for( k = 0; k < M; k++ )
                                {
                                    state->perm[l++] = state->V[imix3+i]*M + k;
                                    state->perm[l++] = state->V[ibad+i]*M + k;
                                    state->perm[l++] = state->V[imix1+i]*M + k;
                                    state->perm[l++] = state->V[imix2+i]*M + k;
                                }
                                break;
                            case 3:
                                for( k = 0; k < M; k++ )
                                {
                                    state->perm[l++] = state->V[ibad+i]*M + k;
                                    state->perm[l++] = state->V[imix3+i]*M + k;
                                    state->perm[l++] = state->V[imix1+i]*M + k;
                                    state->perm[l++] = state->V[imix2+i]*M + k;
                                }
                                break;
                            }
                        }
                        
                        break;
                }
                
            }   //c0_mod != 0 && n0_mod == 0 
            
            #ifndef SKIP_MEX
            // testing
            mexPrintf("Test permutation, Ncode  = %d\n", state->Ncode);
            for( i = 0; i < state->Ncode; i++ )
            {
                if( state->perm[i] == -1 )
                    mexPrintf("i=%d\n",i);
            }
            
            {
                FILE *fp = fopen("perm.txt","wt");
           
                for( i = 0; i < state->Ncode; i+= halfmlog*2 )
                {
                    int k;
                    for( k = 0; k < halfmlog; k++ )
                        fprintf(fp, "%5d ", state->perm[i+k] );
                    fprintf( fp, "  " );
                    for( k = halfmlog; k < mlog; k++ )
                        fprintf(fp, "%5d ", state->perm[i+k] );
                    fprintf(fp, "\n");
                }
                fclose(fp);
            }           
            
// check distribution of 'bad'  positions
            {
                int badfrq[7] = {0};
                int sum = 0;
                int badlim = mrp * M;
                for( i = 0; i < state->Ncode; i+= 3 )
                {
                    int ind = 0;
                    if( state->perm[i] < badlim ) {ind |= 1; sum++; }
                    if( state->perm[i+1] < badlim ) { ind |= 2; sum++; }
                    if( state->perm[i+2] < badlim ) { ind |= 4; sum++; }
                    if( ind )
                        badfrq[ind-1]++;
                }
                mexPrintf("sum = %d\n", sum );
                mexPrintf("frq:  %d %d %d %d %d %d %d\n", badfrq[0],  badfrq[1],  badfrq[2],  badfrq[3],  badfrq[4],  badfrq[5],  badfrq[6] ) ;        
            }    
            
            #endif  //1
            
            break;

        case 3:     // block mode
            random_perm_gen( state->perm, state->block_size );
			if( state->ShortBlockSize )
			{
				random_perm_gen( state->perm_short, state->ShortBlockSize );
			}
            break;
        
        case 4:     // step mode
            random_perm_gen( state->perm, state->block_size );
			if( state->ShortBlockSize )
			{
				random_perm_gen( state->perm_short, state->ShortBlockSize );
			}
            break;
        }
        
        switch( state->perm_mode )
        {
            // create invperm
        case 0:
        case 1:
        case 2:
            for( i = 0; i < state->Ncode; i++ )
            {
                state->temp[i].val = state->perm[i];
                state->temp[i].ind = i;
            }
        
            qsort(  state->temp, state->Ncode, sizeof(PERM_FOR_ORDER), permcomp );
        
            for( i = 0; i < state->Ncode; i++)
                state->invperm[i] = state->temp[i].ind;
            break;
        case 3:
		case 4:
            for( i = 0; i < state->block_size; i++ )
            {
                state->temp[i].val = state->perm[i];
                state->temp[i].ind = i;
            }
        
            qsort(  state->temp, state->block_size, sizeof(PERM_FOR_ORDER), permcomp );
        
            for( i = 0; i < state->block_size; i++)
                state->invperm[i] = state->temp[i].ind;
			if( state->ShortBlockSize )
			{
	            for( i = 0; i < state->ShortBlockSize; i++ )
		        {
			        state->temp[i].val = state->perm_short[i];
				    state->temp[i].ind = i;
				}
        
				qsort(  state->temp, state->ShortBlockSize, sizeof(PERM_FOR_ORDER), permcomp );
        
				for( i = 0; i < state->ShortBlockSize; i++)
					state->invperm_short[i] = state->temp[i].ind;
			}

            break;
#if 0
		case 4:
            for( i = 0; i < state->block_size; i++ )
            {
                state->temp[i].val = state->perm[i];
                state->temp[i].ind = i;
            }
        
            qsort(  state->temp, state->block_size, sizeof(PERM_FOR_ORDER), permcomp );
        
            for( i = 0; i < state->block_size; i++)
                state->invperm[i] = state->temp[i].ind;
                break;
                    
#endif        
        }
    
}


void Permutation( PERMSTATE* state, int direction, double pinp[], double pout[] )
{
	int i;
        switch( state->perm_mode )
        {
        case 0:
        case 1:
        case 2:
            if( direction == 0 )
            {
                // direct permutation
            
                for( i = 0; i < state->Ncode; i++ )
                    pout[i] = pinp[ state->perm[i] ]; 
            
            }
            else
            {
                for( i = 0; i < state->Ncode; i++ )
                    pout[i] = pinp[ state->invperm[i] ]; 
            }
            break;
        case 3:   
            if( direction == 0)
            {
                for( int k = 0; k < state->Nblocks; k++ )
                {
                    int ind  = k * state->block_size;
                    memcpy( state->tmpbuf, pinp+ind, state->block_size * sizeof( pinp[0] ));
                    for( i = 0; i < state->block_size; i++ )
                    {
                        pout[ind + i] = state->tmpbuf[state->perm[i] ];
                    }
                }
				if( state->ShortBlockSize )
				{
					int ind  = state->Nblocks * state->block_size;
                    memcpy( state->tmpbuf, pinp+ind, state->ShortBlockSize * sizeof( pinp[0] ));

                    for( i = 0; i < state->ShortBlockSize; i++ )
                    {
                        pout[ind + i] = state->tmpbuf[state->perm_short[i] ];
                    }
				}
            }
            else
            {
                for( int k = 0; k < state->Nblocks; k++ )
                {
                    int ind  = k * state->block_size;
                    memcpy( state->tmpbuf, pinp+ind, state->block_size * sizeof( pinp[0] ));
                    for( i = 0; i < state->block_size; i++ )
                    {
                        pout[ind + i] = state->tmpbuf[state->invperm[i] ];
                    }
                }
				if( state->ShortBlockSize )
				{
					int ind  = state->Nblocks * state->block_size;
                    memcpy( state->tmpbuf, pinp+ind, state->ShortBlockSize * sizeof( pinp[0] ));

                    for( i = 0; i < state->ShortBlockSize; i++ )
                    {
                        pout[ind + i] = state->tmpbuf[state->invperm_short[i] ];
                    }
				}
            }
            break;
        case 4:    
            if( direction == 0 )
            {
                int j;
                for( int k = 0; k < state->Nblocks; k++ )
                {
                    for( i = 0, j = k; i < state->block_size; i++, j += state->step_size )
                        state->tmpbuf[i] = pinp[j];
                    
                    for( i = 0, j = k; i < state->block_size; i++, j += state->step_size )
                        pout[j] = state->tmpbuf[state->perm[i]];
                }
				if( state->ShortBlockSize )
				{

                    for( i = 0, j = state->Nblocks; i < state->ShortBlockSize; i++, j += state->step_size )
                        state->tmpbuf[i] = pinp[j];
                    
                    for( i = 0, j = state->Nblocks; i < state->ShortBlockSize; i++, j += state->step_size )
                        pout[j] = state->tmpbuf[state->perm_short[i]];
				}
            }
            else
            {
                int j;
                for( int k = 0; k < state->Nblocks; k++ )
                {
                    for( i = 0, j = k; i < state->block_size; i++, j += state->step_size )
                        state->tmpbuf[i] = pinp[j];
                    
                    for( i = 0, j = k; i < state->block_size; i++, j += state->step_size )
                        pout[j] = state->tmpbuf[state->invperm[i]];
                }
				if( state->ShortBlockSize )
				{

                    for( i = 0, j = state->Nblocks; i < state->ShortBlockSize; i++, j += state->step_size )
                        state->tmpbuf[i] = pinp[j];
                    
                    for( i = 0, j = state->Nblocks; i < state->ShortBlockSize; i++, j += state->step_size )
                        pout[j] = state->tmpbuf[state->invperm_short[i]];
				}
            }
            break;
        }   // end of switch


}




#ifndef SKIP_MEX
//
//  Call for SETUP
//  ==============
//  result = Direct_Inverse_Perm( HC, b, c, M, QAM, perm_mode );
//
//  result == 1, if OK
//
//  Call for direct permutation
//  ===========================
//  permcode = Direct_Inverse_Perm( codeword, 0 );
//  codeword - input
//  permcode - output
//
//  Call for inverse permutation
//  ============================
//  codeword = Direct_Inverse_Perm( permcode, 1 );
//  permcode - input
//  codeword - output
//
//  Call for Close module
//  =====================
//  Direct_Inverse_Perm();
//
void mexFunction(int nOut, mxArray *pOut[], int nInp, const mxArray *pInp[])
{
    static PERMSTATE* state = NULL;
    int QAM;
    int b,c,M;
    int halfmlog;
    int i,j,l;
    int perm_mode;
    int block_size = 0, step_size = 0;
    
    if( nInp >= 6 )
    {
        if( nOut != 1 )
        {
            mexErrMsgTxt("Only one output argument allowed for SETUP");
        }
        
        if( state )
        {
            mexPrintf("Close after restart\n");
            Permutations_Close( state );
        }
        // SETUP for permutation module
        QAM = (int)mxGetPr(pInp[4])[0];
        int indQ = findQ( QAM );
        if( indQ == NOT_FOUND )
        {
            mexErrMsgTxt("Only QAM-4, 16, 64, 256 are allowed");
        }
        int mlog = logQav[indQ];
        int halfmlog = mlog/2;
    
        b = (int)mxGetPr(pInp[1])[0];
        c = (int)mxGetPr(pInp[2])[0];
        M = (int)mxGetPr(pInp[3])[0];
        perm_mode = (int)mxGetPr(pInp[5])[0];
        
        if( perm_mode == 3 )
        {
            if( nInp < 8 )
                mexErrMsgTxt("Block_size is needed");
            block_size = (int)mxGetPr(pInp[6])[0];
            if( (c*M) % block_size )
                mexErrMsgTxt("Block_size must divide codelength");
        }   

        if( perm_mode == 4 )
        {
            if( nInp < 8 )
                mexErrMsgTxt("Step_size is needed");
            step_size = (int)mxGetPr(pInp[7])[0];
            if( (c*M) % step_size )
                mexErrMsgTxt("Step_size must divide codelength");
        }   
        
        mexPrintf("b = %d c = %d M=%d halfmlog=%d QAM = %d\n",b,c,M,halfmlog,QAM);
        
        state = Permutations_Open( b, c, M, QAM, halfmlog, perm_mode, block_size, step_size );
        if( !state ) 
             mexErrMsgTxt("Allocation error");
        
        short** hc = Alloc2d_short( b, c);
        if( !hc )
            mexErrMsgTxt("Allocation error");
        
        unpackMatrix_double2short(mxGetPr(pInp[0]), b, c, hc );
        
        
        Permutation_Init( state, hc );
       
        pOut[0]=mxCreateDoubleMatrix(1,1,mxREAL); 
        *(mxGetPr(pOut[0])) = 1;        // OK
        return;
        
    }
    
    if( nInp == 2 )
    {
        double *pout;
        double *pinp;
        int direction;
        
        if( state == NULL )
            mexErrMsgTxt("Decoder not opend"); 
        if( nOut != 1 )
            mexErrMsgTxt("Only one output argument allowed for permutations");
        
        pOut[0]=mxCreateDoubleMatrix(1,state->Ncode,mxREAL);
        pout = mxGetPr(pOut[0]);
        pinp = mxGetPr(pInp[0]);
        
        direction = (int)mxGetPr(pInp[1])[0];
        
        Permutation( state, direction, pinp, pout ); 
    }
    
    
    if( nInp == 0 && state != NULL )
    {
        // Check if dynamic arrays were allocated
        Permutations_Close( state );
        state = NULL;
    }
    

}
#endif
