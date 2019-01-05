#ifndef __COMMONS_PORTABLE_H
#define __COMMONS_PORTABLE_H

#include <stdio.h>
#include <stdlib.h>

/*
 * This header file contains different utilities
 * which need to be portable but often are not.
 *
 * For each such utility, a distinct name is picked up
 * that is different from the implementations,
 * and then implemented in the best way
 * according to the platform/compiler.
 */

// Meta
void commons_portable_print_info();

// Randomness
extern int initial_random_seed;
int next_random_int(int minInclusive, int maxExclusive);
double next_random_normal();

// Console interactive input
bool is_keyboard_hit();
char get_character();

// File operations
FILE *fopen_or_die(char const *filename, char const *mode, char const *error_message);
void skip_until_newline(FILE *handle);

// Macros

// fscanf_1 does the same fscanf with one target argument does
// additionally checks if this argument was read successfully
// and fails with a meaningful message if not.
#define fscanf_1(handle, fmt_string, param_1)                               \
    if (fscanf(handle, fmt_string, &param_1) != 1) {                        \
        fprintf(                                                            \
            stderr,                                                         \
            "Scanning parameter " #param_1 " failed, file %s, line %d\n",   \
            __FILE__, __LINE__                                              \
        );                                                                  \
        exit(1);                                                            \
    }
// "" to shutdown buggy highlighting in Nano


#endif
