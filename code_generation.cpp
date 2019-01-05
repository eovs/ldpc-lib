#include <stdio.h>
#include <exception>

#include "commons_portable.h"
#include "code_generation.h"

using std::vector;

class interruption_exception : public std::exception {};

bool check_same_colums( matrix< int > const &code_matrix )
{
    int R = code_matrix.n_rows();
    int N = code_matrix.n_cols();
	

	for( int j = 0; j < N; j++ )
	{
		// compare j-th column with all others from j+1 to N
		int i;
		for( int jc = j+1; jc < N; jc ++ )
		{
			for( i = 0; i < R; i++ )
			{
				if( code_matrix(i,j) != code_matrix(i,jc) )
					break;
			}
			if( i == R )
				return true;
		}
	}
	return false;
}

bool generate_code(
    matrix< int > const &code_matrix,
    int target_girth,
    int max_module,
    int max_codes_to_test,
    matrix< int > &result_matrix
) {
    int R = code_matrix.n_rows();
    int N = code_matrix.n_cols();

    graph< vector_bag > g;
    vector< bit > set_zeros;
    int voltage_size = 0;

    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < R; ++i) {
            if (code_matrix(i, j) > 0) {
                vector_bag fw(voltage_size++);
                g.add_bidi_edge(j, N + i, fw, -fw);
                set_zeros.push_back(code_matrix(i, j) != 1);
            }
        }
    }

    vector< int > voltages(voltage_size);

    printf("\n n=%d, r=%d, g=%d, M=%d\n", N, R, target_girth, max_module);
    printf("\n Tanner graph size: nt=%d, rt=%d\n", voltage_size, g.n_vertices());

    printf("Generating first equations\n");
    equation_builder b(g, target_girth);
    if (b.tree_girth() < target_girth) {
        return false;
    }
    printf("ok: g=%d, eqs=%d, nodes=%d\n", b.tree_girth(), b.n_equations(), b.n_vertices());

    printf("Creating random %d columns\n", N);

    try {
        console_exception_hook x_hook('x', "interrupts current code generation", interruption_exception());
        for (int attempt = 0; attempt < max_codes_to_test; ++attempt) {
            for (int i = 0; i < voltage_size; ++i) {
                voltages[i] = set_zeros[i] ? 0 : next_random_int(0, max_module);
            }
            bool ok = b.check_voltages(voltages).check_modulo(max_module);
            if (ok) {
//                printf("\n  1 codes selected, %d tested\n", attempt + 1);
                result_matrix = code_matrix;
                for (int col = 0, vp = 0; col < N; ++col) {
                    for (int row = 0; row < R; ++row) {
                        if (code_matrix(row, col) == 0) {
                            result_matrix(row, col) = -1;
                        } else {
                            result_matrix(row, col) = voltages[vp++] % max_module;
                        }
                    }
                }
				// checking the same column
				bool same_cols = check_same_colums( result_matrix );
				if( same_cols )
				{
					printf("=============================> matrix with same columns\n" );
					continue;
				}
				else
				{
	                printf("\n  1 codes selected, %d tested\n", attempt + 1);
					return true;
				}
            }

            console_message_hook i_hook('i', "prints how much codes is tested", format_to_string("    n = %d, %d tested, 0 found\n", N, attempt + 1));
            console_check_hooks();

            if ((attempt + 1) % 100000 == 0) {
                printf("    %d tested\n", attempt + 1);
            }
        }
    } catch (interruption_exception &) {}
    return false;
}

