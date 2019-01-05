#ifndef _DECODERS_H_
#define _DECODERS_H_


//#define ORIG_TABLES

#define DEC_DECISION 0	// 0 - hard decision, 1 - soft decision

//#define BP_USE_EPS	// be careful!!!

#define IASP_FIXED_POINT

#define MS_MUL_CORRECTION
#define MS_ALPHA_FPP 4

enum  DEC_ID
{
	BP_DEC, 
	SP_DEC, 
	ASP_DEC, 
	MS_DEC,
	IMS_DEC,
	IASP_DEC,
	FHT_DEC,
	TASP_DEC
};

extern char const * const DEC_FULL_NAME[];


#define SKIP  -1
#define TRUE_CIRCULANT


#ifdef MS_MUL_CORRECTION
#define MS_ALPHA 0.95
#else
#define MS_ALPHA 0.02//0.075//0.02
#endif

#define MS_THR   1.4
#define MS_QBITS 6//5//6
#define MS_DBITS (MS_QBITS + 2) //16

#define UFLT_MNT_16
#define FLT_MNT_16
//#define FLT_POW_16

#ifdef _MSC_VER
typedef unsigned __int64    ui64;
typedef __int64             i64;
#else
typedef unsigned long long  ui64;
typedef long long           i64;
#endif

typedef unsigned int        ui32;
typedef unsigned short      ui16;
typedef short               i16;
typedef int                 i32;
typedef unsigned char       ui8;
typedef signed char         i8;


#ifdef UFLT_MNT_16
typedef ui16	UFLT_MNT;
#else
typedef ui8	UFLT_MNT;
#endif

#ifdef FLT_MNT_16
typedef i16	FLT_MNT;
#else
typedef i8	FLT_MNT;
#endif

#ifdef FLT_POW_16
typedef i16		FLT_POW;
#else
typedef i8		FLT_POW;
#endif

typedef struct  
{
	UFLT_MNT m;
	FLT_POW p;
} UFLT16;

typedef struct  
{
	FLT_MNT m;
	FLT_POW p;
} FLT16;

#if 01
typedef double MS_DATA;
//typedef unsigned short UDATA;
#else
typedef char MS_DATA;
typedef unsigned char UDATA;
#endif

#if 01
typedef short IMS_DATA;
#else
typedef char IMS_DATA;
#endif


typedef struct  
{
	MS_DATA min1; 
	MS_DATA min2; 
	int pos;
	int sign;
} MS_DEC_STATE;

typedef struct  
{
	IMS_DATA min1; 
	IMS_DATA min2; 
	int pos;
	int sign;
} IMS_DEC_STATE;


typedef struct
{
	int q_bits;
	int q;
	int nh; 
	int rh;
	int m;
    int n;
	int maxiter;
    int codec_id;
	int bin_codec;
	
	int max_rw;

	short   **hd; 
	short   **hb; 
	short   **hc; 
	double  *y;
	double  *decword;
	short   *syndr;                     
	double  **qy;
	double  **qdecword;
	short   *qhard;

	// Belief Propagation Decoder
	double  *bp_data;      
	short   *bp_sign;     
	short   *bp_bs;        
	double  *bp_yd;                   
	double  *bp_s;                 
	double  **bp_ZZ;       
	short   **bp_BB;      

	// Sum-Prod Decoder
	double  *sp_data;
	double  *sp_yd;  
	double  *sp_s;   
	double  **sp_ZZ; 
	double  **sp_ZZ0;
	double  *sp_AA;  

	// Advanced Sum-Prod Decoder
	double  *asp_data0;
	double	*asp_data1;
	double  *asp_p0;
	double	*asp_p1;
	double  *asp_soft_out;
	int		*asp_posh;
	int     *asp_rw;
	int		**asp_hc_ri;
	double  **asp_state;
	int		asp_all_cw_2;

	// Integer Advanced Sum-Prod Decoder
#ifdef IASP_FIXED_POINT	
	ui16 *iasp_y;
	ui16 *iasp_data0;
	ui16 *iasp_data1;
	ui16 **iasp_state;
	ui32 *iasp_p0;
	ui32 *iasp_p1;
#else
	UFLT16 *iasp_y;
	ui16  *iasp_data;
	UFLT16 *iasp_data0;
	UFLT16 *iasp_data1;
	UFLT16 **iasp_state;
	UFLT16 *iasp_p0;
	UFLT16 *iasp_p1;
#endif
	ui16 *iasp_soft_out;
	int	*iasp_posh;
	int *iasp_rw;
	int	**iasp_hc_ri;
	int	iasp_all_cw_2;

	// Min-Sum Decoder
	MS_DATA *ms_soft;
	short	*ms_BnNS;
	MS_DEC_STATE *ms_dcs;
	MS_DEC_STATE *ms_tmps;
	MS_DATA *ms_buffer;
	MS_DATA *ms_rbuffer;
	MS_DATA *ms_rsoft;

	// Integer Min-Sum Decoder
	IMS_DATA *ims_soft;
	short	 *ims_BnNS;
	IMS_DATA *ims_y;
	IMS_DEC_STATE *ims_dcs;
	IMS_DEC_STATE *ims_tmps;
	IMS_DATA *ims_buffer;
	IMS_DATA *ims_rbuffer;
	IMS_DATA *ims_rsoft;

	// FHT Sum-Product
	int *fht_rw;
	double *fht_HAD;
	short *fht_smask;
	int	fht_all_cw_2;
	double fht_p_thr;
	int fht_ilist_size;
	short **fht_hb_ci;
	short **fht_hb_ri;
	short **fht_hc_rl;
	short **fht_hb_cj;
	short **fht_t;
	short **fht_ti;
	int *fht_list;
	int *fht_ilist;
	short *fht_buf;
	short *fht_pos;
	double **fht_buf0;
	double **fht_buf1;
	double **fht_soft_in;
	double **fht_soft_outs;
	double **fht_soft_out;

	// TDMP Advanced Sum-Prod Decoder
	double  *tasp_data0;
	double	*tasp_tmp;
//	double  *tasp_p0;
//	double	*tasp_p1;
	double  *tasp_soft_out;
//	int		*tasp_posh;
//	int     *tasp_rw;
//	int		**tasp_hc_ri;
	double  **tasp_state;
//	int		tasp_all_cw_2;

}DEC_STATE;


// functions prototypes
DEC_STATE* decod_open( int decoder_id, int q_bits, int mh, int nh, int M );
int decod_init( void* st );
void decod_close( DEC_STATE* st );
int bp_decod_qc_lm( DEC_STATE* st, double soft[], double decword[], int maxiter, int decision);
int sum_prod_decod_qc_lm( DEC_STATE* st, double soft[], double decword[], int maxiter, int decision );    
int sum_prod_gf2_decod_qc_lm( DEC_STATE* st, double soft[], double decword[], int maxiter, int decision );    
int min_sum_decod_qc_lm( DEC_STATE* st, double soft[], double decword[], int maxiter, int decision, double alpha );    
int imin_sum_decod_qc_lm( DEC_STATE* st, double soft[], double decword[], int maxiter, int decision, double alpha, double thr, int qbits, int dbits );    
int isum_prod_gf2_decod_qc_lm( DEC_STATE* st, double soft[], double decword[], int maxiter, int decision );    
int sum_prod_gfq_decod_lm( DEC_STATE* st, double *soft[], short *qhard, double *decword[], int maxiter, double p_thr );
int tdmp_sum_prod_gf2_decod_qc_lm( DEC_STATE* st, double soft[], double decword[], int maxiter, int decision );    




#endif	//_DECODERS_H_
