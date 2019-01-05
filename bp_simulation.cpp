#include <algorithm>
#include <cstring>
#include <cmath>
#include <exception>
#include <vector>

#include "bp_simulation.h"
#include "commons_portable.h"
#include "decoders.h"
#include "modulation.h"

using std::pair;
using std::vector;
using std::make_pair;
using std::log;
using std::exp;
using std::memset;

class interruption_exception : public std::exception {};

int qc_encode(
    matrix< int > const &mx,
    int tailbite_length,
    vector< bit > &cword
) {
    int b = mx.n_rows();
    int c = mx.n_cols();
    int r = b * tailbite_length;
    int n = c * tailbite_length;

    bool is_single_diagonal = mx(1, 0) < 0;

    // 1. Find a positive degree in (b-1)-th column.
    int p = 0;
    while (p < b && mx(p, b - 1) <= 0) {
        ++p;
    }
    if (p >= b && !is_single_diagonal) {
        return -1; // no positive degree
    }

    if ((int) cword.size() != n) {
        die("Input word has the length %d but should have been %d", (int) cword.size(), n);
    }

    // 2. Computing partial syndromes.
    vector< bit > synd(r), sumsynd(tailbite_length);
    for (int i = 0; i < b; ++i) {
        for (int j = b; j < c; ++j) {
            if (mx(i, j) >= 0) {
                int offset = mx(i, j);
                for (int h = 0; h < tailbite_length; ++h) {
                    bool cword_bit = cword[j * tailbite_length + (h + offset) % tailbite_length];
                    synd[i * tailbite_length + h] ^= cword_bit;
                }
            }
        }
        for (int h = 0; h < tailbite_length; ++h) {
            sumsynd[h] ^= synd[i * tailbite_length + h];
        }
    }

    if (is_single_diagonal) {
        // 3a. Single diagonal, copying the partial syndrome onto the codeword
        for (int i = 0; i < r; ++i) {
            cword[i] = synd[i];
        }
    } else {
        // 3b. Building up the rest of the codeword
        for (int h = 0; h < tailbite_length; ++h) {
            bool xh = sumsynd[(h + tailbite_length - mx(p, b - 1)) % tailbite_length];
            cword[(b - 1) * tailbite_length + h] = xh;
            cword[h] = synd[h];
            if (mx(0, b - 1) == 0) cword[h] ^= xh;
            if (mx(0, b - 1) >  0) cword[h] ^= sumsynd[h];
            for (int i = 1; i < b - 1; ++i) {
                int idx = i * tailbite_length + h;
                cword[idx] = synd[idx] ^ cword[idx - tailbite_length];
                if (mx(i, b - 1) == 0) cword[idx] ^= cword[(b - 1) * tailbite_length + h];
                if (mx(i, b - 1) >  0) cword[idx] ^= sumsynd[h];
            }
        }
    }

    // 5. Check for being a codeword.
    for (int i = 0; i < r; ++i) {
        synd[i] = 0;
    }
    for (int i = 0; i < b; ++i) {
        for (int j = 0; j < c; ++j) {
            if (mx(i, j) >= 0) {
                int offset = mx(i, j);
                for (int h = 0; h < tailbite_length; ++h) {
                    bool cword_bit = cword[j * tailbite_length + (h + offset) % tailbite_length];
                    synd[i * tailbite_length + h] ^= cword_bit;
                }
            }
        }
    }

    for (int i = 0; i < r; ++i) {
        if (synd[i]) {
            int jj = r;
            FILE *f_cw = fopen_or_die("cw.txt", "wt", "cw.txt cannot be opened");
            for (int j = b; j < c; ++j) {
                for (int i = 0; i < tailbite_length; ++i) {
                    fprintf(f_cw, "%2d", (bool) (cword[jj++]));
                }
                fprintf(f_cw, "\n");
            }
            fclose(f_cw);
            return 1;
        }
    }
    return 0;
}

void independent_validation(matrix< int > const &mx, int tailbite_length, vector< bit > const &codeword) {
    int rows = mx.n_rows(), cols = mx.n_cols();
	int syndrome_size = rows * tailbite_length;
    vector< bit > syndrome(syndrome_size);
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            if (mx(r, c) >= 0) {
                int offset = mx(r, c);
                for (int t = 0; t < tailbite_length; ++t) {
                    bool cword_bit = codeword[c * tailbite_length + (t + offset) % tailbite_length];
                    syndrome[r * tailbite_length + t] ^= cword_bit;
                }
            }
        }
    }
	// fixed bug: c < cols ===> c < syndrome_size
    for (int c = 0; c < syndrome_size; ++c) {
        if (syndrome[c]) {
            die("Codeword is bad [independent validation]");
        }
    }
}

int random_codeword(matrix< int > const &mx, int tailbite_length, vector< bit > &codeword) {
    vector< int > break_items;
    break_items.push_back(0);
    for (int i = 1, i_max = mx.n_rows(); i + 1 < i_max; ++i) {
        // The bidiagonal block begins
        if (mx(i, i) >= 0 && mx(i + 1, i) >= 0 && mx(i, i - 1) < 0) {
            break_items.push_back(i);
        }
        // The unidiagonal block begins
        if (i > 1 && mx(i, i) >= 0 && mx(i + 1, i) < 0 && mx(i, i - 1) < 0 && mx(i - 1, i - 1) >= 0 && mx(i - 1, i - 2) >= 0) {
            break_items.push_back(i);
        }
    }
    break_items.push_back(mx.n_rows());

    int n = mx.n_cols() * tailbite_length;

    vector< bit > cword(n);
    for (int i = break_items.back() * tailbite_length; i < n; ++i) {
        cword[i] = next_random_int(0, 2) == 1;
    }

    for (int i = (int) (break_items.size()) - 2; i >= 0; --i) {
        int curr_length = break_items[i + 1];
        int next_length = break_items[i];
        int offset = next_length;
        matrix< int > mx2(curr_length - next_length, mx.n_cols() - offset);
        for (int r = 0, r_max = mx2.n_rows(); r < r_max; ++r) {
            for (int c = 0, c_max = mx2.n_cols(); c < c_max; ++c) {
                mx2(r, c) = mx(r + offset, c + offset);
            }
        }
        vector< bit > local_cword;
        for (int z = offset * tailbite_length, j = z; j < n; ++j) {
            local_cword.push_back(cword[j]);
        }
        int qc_result = qc_encode(mx2, tailbite_length, local_cword);
        if (qc_result != 0) {
            return qc_result;
        }
        for (int z = offset * tailbite_length, j = z; j < n; ++j) {
            cword[j] = local_cword[j - z];
        }
    }
    codeword = cword;

    independent_validation(mx, tailbite_length, codeword);

    return 0;
}

#if 1
int bp_decod_lm(
    vector< double > &soft,
    matrix< int > const &C,
    matrix< int > const &V,
    matrix< int > const &VI,
    vector< int > const &col_weights,
    vector< int > const &row_weights,
    int max_iterations
) {
    int r = V.n_rows();
    int n = C.n_rows();
    const double eps = 2e-5;

    // to be initialized to the appropriate size lazily
    vector< double > y;
    vector< double > s;
    vector<  bit   > bs;
    matrix< double > Z;
    matrix<  bit   > B;

    for (int iter = 0; iter < max_iterations; ++iter) {
        int synd = 0;
        for (int i = 0; i < r; ++i) {
            synd = 0;
            for (int j = 0; j < row_weights[i]; ++j) {
                synd ^= soft[V(i, j)] < 0;
            }
            if (synd) {
                break;
            }
        }
        if (!synd) {
            return iter;
        }
        // First time there, so init arrays
        if (iter == 0) {
            int max_col_weight = 0;
            for (int i = 0; i < n; ++i) {
                soft[i] = std::max(std::min(soft[i], 20.0), -20.0);
                max_col_weight = std::max(max_col_weight, col_weights[i]);
            }
            y = soft;
            s.resize(r);
            bs.resize(r);
            Z = matrix< double >(n, max_col_weight);
            B = matrix<  bit   >(n, max_col_weight);
        }
        // Variable-node activation step
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < col_weights[i]; ++j) {
                double a = exp(soft[i] - Z(i, j));
                B(i, j) = a < 1;
                double x = log(std::abs((a - 1) / (a + 1)));
                if (std::abs(x) < eps) {
                    x = x <= 0 ? -eps : eps;
                }
                Z(i, j) = x;
            }
        }
        // Check-node activation step
        for (int i = 0; i < r; ++i) {
            s[i] = 0;
            bs[i] = false;
            for (int j = 0; j < row_weights[i]; ++j) {
                s[i]  += Z(V(i, j), VI(i, j));
                bs[i] ^= B(V(i, j), VI(i, j));
            }
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < col_weights[i]; ++j) {
                double a = exp(s[C(i, j)] - Z(i, j));
                bool b = bs[C(i, j)] ^ B(i, j);
                Z(i, j) = std::max(-19.07, std::min(19.07, (1 - 2 * b) * log((1 + a) / (1 - a))));
            }
        }
        // Output of an iteration
        for (int i = 0; i < n; ++i) {
            double s = y[i];
            for (int j = 0; j < col_weights[i]; ++j) {
                s += Z(i, j);
            }
            soft[i] = s;
        }
    }
    return -max_iterations;
}
#endif

pair<double, double> bp_simulation(
    matrix< int > const &code_generating_matrix,
    int tailbite_length,
    int max_iterations,
    int n_frame_errors,
    int n_experiments,
    double snr,
    double reference_frame_error,
	int decoder_type,
	int modulation_type,
	int permutation_type,
	int permutation_block,
	int permutation_inter,
	int punctured_blocks,
	int show_process
) {
    int b = code_generating_matrix.n_rows();
    int c = code_generating_matrix.n_cols();
    int r = b * tailbite_length;
    int n = c * tailbite_length;
	int M = tailbite_length;
    int i, j;

	    // 2. Channel model.
	int QAM;				// config
	int block_size = 1;		// config
	int step_size = 1;		// config
	int halfmlog;			// 0.5*log2(QAM )
	int extra_bits = 0;		// if n%(2*halfmlog) != 0
	int n_ext = n;			//
	double T = 26.0;	//?????????????????????????????????????????????

    static DEC_STATE *dec_state =  NULL;
	static PERMSTATE* perm_state = NULL;

	if( dec_state )
		decod_close( dec_state );

	dec_state = decod_open( decoder_type, 1, b, c, tailbite_length );
	if( !dec_state )
		die("Can't open decode");

	for( i = 0; i < b; i++ )
        for( j = 0; j < c; j++ )
            dec_state->hd[i][j] = code_generating_matrix(i, j);

	decod_init( dec_state );


	//modulation_type = MODULATION_QAM16;

	switch( modulation_type )
	{
	case MODULATION_SKIP:	QAM = 1;	halfmlog = 1;	break;
	case MODULATION_QAM4:	QAM = 4;	halfmlog = 1;	break;
	case MODULATION_QAM16:	QAM = 16;	halfmlog = 2;	break;
	case MODULATION_QAM64:	QAM = 64;	halfmlog = 3;	break;
	case MODULATION_QAM256:	QAM = 256;	halfmlog = 4;	break;
	default:                QAM = -1;   halfmlog = -1; die("Unknown modulation type: %d", modulation_type); break;
	}

	if( perm_state )
		Permutations_Close( perm_state );
	perm_state = Permutations_Open( b, c, M, QAM, halfmlog, permutation_type, permutation_block, permutation_inter );
	if( !perm_state )
		die("Can't open permutations");
	
	Permutation_Init( perm_state, dec_state->hd );

	int moddim = 2 * halfmlog;

	if(  n % moddim )
		extra_bits = moddim - n%moddim;   
    n_ext = n + extra_bits;


	static QAM_MODULATOR_STATE* qam_mod_st = NULL;	//QAM_modulator_open( QAM, n_ext, moddim );
	if( qam_mod_st )
		QAM_modulator_close(qam_mod_st);
	qam_mod_st = QAM_modulator_open( QAM, n_ext, moddim );
	if( !qam_mod_st )
		die("Can't open QAM modulator");


    //double bitrate = (double) (c - b) / c;
	double bitrate = (double) (c - b) / (c - punctured_blocks);
	double sigma = sqrt(pow(10, -snr / 10) / 2 / bitrate);

	double norm_factor = 2.0*(QAM-1.0)/3.0;

	double sigmaQAM = sqrt(pow(10., -snr / 10.) /(2 * bitrate * halfmlog*2) * norm_factor);
	int nQAMSignals = n_ext / moddim;
	int out_type = 0;	// 0 for BP_decoder and minsum decoder, 1 - for sumprod

    switch( decoder_type )
    {
    case BP_DEC:    out_type = 0;		break;
    case SP_DEC:    out_type = 1;		break;
    case ASP_DEC:   out_type = 1;		break;
    case MS_DEC:    out_type = 0;		break;
    case IMS_DEC:   out_type = 0;		break;
    case IASP_DEC:  out_type = 1;		break;
	case TASP_DEC:  out_type = 1;		break;
	case LMS_DEC:   out_type = 0;		break;
	case LCHE_DEC:  out_type = 1;		break;
    default: out_type = -1; die("Unknown decoder type: %d", decoder_type); break;
    }

	static QAM_DEMODULATOR_STATE* qam_demod_st = NULL;	//QAM_demodulator_open( T, sigmaQAM, QAM, n, tailbite_length, nQAMSignals, out_type );
	if( qam_demod_st )
		QAM_demodulator_close(qam_demod_st );
	qam_demod_st = QAM_demodulator_open( T, sigmaQAM, QAM, n, tailbite_length, nQAMSignals, out_type );
	if( !qam_demod_st )
		die("Can't open QAM demodulator");


#if 0

    // 1. Generating matrices.
    vector< int > col_weights(n), row_weights(r);
    matrix< int >  C(n, b);
    matrix< int >  V(r, c);
    matrix< int > VI(r, c);

    for (int i = 0; i < b; ++i) {
        for (int j = 0; j < c; ++j) {
            if (code_generating_matrix(i, j) >= 0) {
                for (int h = 0; h < tailbite_length; ++h) {
                    int r_ind = i * tailbite_length + h;
                    int c_ind =
                        j * tailbite_length +
                        (code_generating_matrix(i, j) + h) % tailbite_length;
                    C(c_ind, col_weights[c_ind]) = r_ind;
                    V(r_ind, row_weights[r_ind]) = c_ind;
                    VI(r_ind, row_weights[r_ind]) = col_weights[c_ind];
                    ++col_weights[c_ind];
                    ++row_weights[r_ind];
                }
            }
        }
    }
#endif




    // 3. Generating random codeword.
    vector< bit > codeword;
	bool zero_codeword = false;
    int codeword_exitcode = random_codeword(code_generating_matrix, tailbite_length, codeword);
    if (codeword_exitcode < 0) {
        // MB: shouldn't we die there?
        printf("Bad matrix: zero codeword will be used\n");
        //return make_pair(-10.0, -10.0);
		zero_codeword = true;
    } else if (codeword_exitcode > 0) {
        // MB: shouldn't we die there?
        printf("Bad encoding: zero codeword will be used\n");
        //return make_pair(10.0, 10.0);
		zero_codeword = true;
    }

	if( zero_codeword )
	{
		for( i = 0; i < n; i++ )
			codeword.push_back(0);
	}


	//EUG
	for( i = 0; i < n; i++ ) codeword[i] = 0;

    // 4. Simulation
	for( i = 0; i < n; i++ )
		perm_state->buffer[i] = codeword[i];
	
	Permutation(perm_state, 0, perm_state->buffer, perm_state->perm_codeword );		//direct permutation

	memset( perm_state->perm_codeword+n, 0, extra_bits * sizeof( perm_state->perm_codeword[0] ) );


    int nse = 0, nue = 0, nde = 0;
    int experiment = 0;
	
	

	QAM_modulator( qam_mod_st, perm_state->perm_codeword, qam_mod_st->dx  );

    try {
        console_exception_hook x_hook('x', "interrupts current code simulation", interruption_exception());
        while (nde < n_frame_errors && experiment <= n_experiments) {
            console_check_hooks();
            ++experiment;

			switch( modulation_type )
			{
			case MODULATION_SKIP:
				for ( i = 0; i < n; ++i) {
					double noise = next_random_gaussian();
					perm_state->buffer[i] = -2.0 * (sigma * noise + 2.0 * perm_state->perm_codeword[i] - 1.0) / (sigma * sigma);
				}
				break;

			case MODULATION_QAM4:
				for ( i = 0; i < n; ++i) {
					double noise = next_random_gaussian();
					perm_state->buffer[i] = -2.0 * (sigmaQAM * noise + 2.0 * perm_state->perm_codeword[i] - 1.0) / (sigmaQAM * sigmaQAM);
				}
				break;

			case MODULATION_QAM16:
				//break;

			case MODULATION_QAM64:
				//break;

			case MODULATION_QAM256:
				for( i = 0; i < 2*nQAMSignals; i++ )
				{
					double noise = next_random_gaussian();
					qam_mod_st->dx[i] += noise; 
				}
				Demodulate( qam_demod_st, qam_mod_st->dx, perm_state->buffer );
				for( i = 0; i < n; i++ )
					perm_state->buffer[i] = -perm_state->buffer[i] ;

				break;
			}

			//memcpy( perm_state->buffer, dec_state->y, n*sizeof(dec_state->y[0]));

			Permutation(perm_state, 1, perm_state->buffer, dec_state->y );		// Inverse permutation

#if 0
            int iter = bp_decod_lm( received, C, V, VI, col_weights, row_weights, max_iterations );
            int curr_nse = 0;
            for (int i = 0; i < n; ++i) {
                curr_nse += (received[i] < 0) != codeword[i];
            }
#else
            int iter, curr_nse = 0, curr_nse_info = 0;

			{
				double init_val = qam_demod_st->DemodOutType == 1 ? 0 : 0.5;

				int plen = dec_state->m * punctured_blocks;
				int pstart = n - plen;
//				int pstart = 0;
//				int pstart = dec_state->rh*dec_state->m;
				int pstop  = pstart + plen;

				for( i = pstart; i < pstop; i++ )
					dec_state->y[i] = init_val;
			}
/*
			for( i = 0; i < n; i++ )
			{
				dec_state->y[i] /= 10;
				if( dec_state->y[i] <= 0 )
					i = i;
			}

			dec_state->y[] = -0.1;
*/
            switch( decoder_type )
            {
            case BP_DEC:    iter = bp_decod_qc_lm( dec_state, dec_state->y, dec_state->decword, max_iterations, DEC_DECISION);				 break;
            case SP_DEC:    iter = sum_prod_decod_qc_lm( dec_state, dec_state->y, dec_state->decword, max_iterations, DEC_DECISION);			 break;
            case ASP_DEC:   iter = sum_prod_gf2_decod_qc_lm( dec_state, dec_state->y, dec_state->decword, max_iterations, DEC_DECISION);		 break;
            case MS_DEC:    iter = min_sum_decod_qc_lm( dec_state, dec_state->y, dec_state->decword, max_iterations, DEC_DECISION, MS_ALPHA); break;
            case IMS_DEC:   iter = imin_sum_decod_qc_lm( dec_state, dec_state->y, dec_state->decword, max_iterations, DEC_DECISION, MS_ALPHA, MS_THR, MS_QBITS, MS_DBITS); break;
            case IASP_DEC:  iter = isum_prod_gf2_decod_qc_lm( dec_state, dec_state->y, dec_state->decword, max_iterations, DEC_DECISION);		 break;
			case TASP_DEC:  iter = tdmp_sum_prod_gf2_decod_qc_lm( dec_state, dec_state->y, dec_state->decword, max_iterations, DEC_DECISION);		 break;
			case LMS_DEC:   iter = lmin_sum_decod_qc_lm( dec_state, dec_state->y, dec_state->decword, max_iterations, DEC_DECISION, MS_ALPHA, MS_BETA); break;
			case LCHE_DEC:  iter = lche_decod( dec_state, dec_state->y, dec_state->decword, max_iterations, DEC_DECISION);		 break;
            default: iter = -1; die("Unknown decoder type: %d", decoder_type); break;
            }

            if( DEC_DECISION == 0 )
            {
                for (i = 0; i < n; i++) {
                    if (dec_state->decword[i] != (short)codeword[i]) {
                        ++curr_nse;
                        if (i >= r) {
                            ++curr_nse_info;
                        }
                    }
                }
				if( curr_nse && iter < max_iterations )
					i = 0;
            }
            else
		    {
                switch( decoder_type )
                {
                case BP_DEC:
                case MS_DEC:
                case IMS_DEC:
				case LMS_DEC:
				case LCHE_DEC:
                    for (i = 0; i < n; i++) {
                        if ((dec_state->decword[i] < 0.0) != codeword[i]) {
                            ++curr_nse;
                            if (i >= r) {
                                ++curr_nse_info;
                            }
                        }
                    }
                    break;
                case SP_DEC:
                    for (i = 0; i < n; i++) {
                        if ((dec_state->decword[i] < 0.5) != codeword[i]) {
                            ++curr_nse;
                            if (i >= r) {
                                ++curr_nse_info;
                            }
                        }
                    }
                    break;
                case ASP_DEC:
                    for (i = 0; i < n; i++) {
                        if ((dec_state->decword[i] > 0.5) != codeword[i]) {
                            ++curr_nse;
                            if (i >= r) {
                                ++curr_nse_info;
                            }
                        }
                    }
                    break;
                case IASP_DEC:
//                    for (i = 0; i < n; i++) curr_nse += (dec_state->decword[i] > 0.5) != codeword[i];
                    break;
                }
            }

#endif
            if (curr_nse > 0) {
                nse += curr_nse_info;
                ++nde;
                if (iter >= 0) {
                    ++nue;
                }

				if( show_process )
				{
					printf("SNR=%5.3lf,step=%4d,s_ers=%d,f_ers=%d,u_ers=%d,BER=%5.3le,FER=%5.3le\n",
								snr, experiment, nse, nde, nue,
								(double) nse / experiment / (n - r),
								(double) nde / experiment);
				}

                if (nde >= 10 && (double) nde / experiment > 2.5 * reference_frame_error) {
                    break;
                }
            }
        }
    } catch (interruption_exception &) {
        // A cancelled code cannot be good
        return make_pair(-1.0, -1.0);
    }


	decod_close( dec_state );
	dec_state = NULL;
	Permutations_Close( perm_state );
	perm_state = NULL;
	QAM_modulator_close( qam_mod_st);
	qam_mod_st = NULL;
	QAM_demodulator_close( qam_demod_st );
	qam_demod_st = NULL;

    return make_pair((double) nse / experiment / (n - r), (double) nde / experiment);
}
