#ifndef __COMMONS_PORTABLE_H
#define __COMMONS_PORTABLE_H

#include <stdio.h>
#include <istream>
#include <string>
#include <exception>

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

#if __cplusplus >= 201103L
[[noreturn]]
void die(char const *format, ...);
#else
void die(char const *format, ...);
#endif

// Format-to-string
std::string format_to_string(char const *format, ...);

// Randomness
extern int initial_random_seed;
void reset_random();
void ensure_random_is_initialized();
int next_random_int(int minInclusive, int maxExclusive);
double next_random_01();
double next_random_gaussian();

// Console interactive input
//// Direct access
bool is_keyboard_hit();
char get_character();

//// Hooks
struct exception_wrapper {
    virtual void execute() const = 0;
    virtual ~exception_wrapper();
};

template<typename exception_type>
class actual_exception_wrapper : public exception_wrapper {
    exception_type my_exception;
public:
    actual_exception_wrapper(exception_type const &exception) : my_exception(exception) {}
    virtual void execute() const { throw my_exception; }
    virtual ~actual_exception_wrapper() {}
};

void console_digest_exception_wrapper(char symbol, std::string const &explanation, exception_wrapper *wrapper);

class console_message_hook {
    char symbol;
public:
    console_message_hook(char symbol, std::string const &explanation, std::string const &message);
    ~console_message_hook();
};
class console_exception_hook {
    char symbol;
public:
    template<typename exception_type>
    console_exception_hook(char symbol, std::string const &explanation, exception_type const &exception) : symbol(symbol) {
        console_digest_exception_wrapper(symbol, explanation, new actual_exception_wrapper<exception_type>(exception));
    }
    ~console_exception_hook();
};
void console_check_hooks();

// File operations
FILE *fopen_or_die(char const *filename, char const *mode, char const *error_message);
void skip_until_newline(FILE *handle);
std::string get_native_path(std::string source);

// Misc I/O: line eater
struct line_eater_type {};
extern line_eater_type line_eater;
std::istream &operator >> (std::istream &stream, line_eater_type const &);

// Bit
struct bit {
    char value;

    inline bit() : value(0) {}
    inline bit(bit const &that) : value(that.value) {}
    inline bit(bool that) : value(that ? 1 : 0) {}
    inline bit(int  that) : value(that ? 1 : 0) {}

    inline bit &operator = (bit const &that) { value = that.value;   return *this; }
    inline bit &operator = (bool that)       { value = that ? 1 : 0; return *this; }
    inline bit &operator = (int  that)       { value = that ? 1 : 0; return *this; }

    inline bit &operator |= (bit that) { value |= that.value; return *this; }
    inline bit &operator &= (bit that) { value &= that.value; return *this; }
    inline bit &operator ^= (bit that) { value ^= that.value; return *this; }

    inline operator bool  () const { return value != 0; }
    inline bool operator !() const { return value == 0; }
    inline bool operator  ~() const { return value == 0; }
};

// Macros

// fscanf_1 does the same fscanf with one target argument does
// additionally checks if this argument was read successfully
// and fails with a meaningful message if not.
#define fscanf_1(handle, fmt_string, param_1)                               \
    if (fscanf(handle, fmt_string, &param_1) != 1) {                        \
        die("Scanning parameter %s failed, file %s, line %d",               \
            #param_1, __FILE__, __LINE__                                    \
        );                                                                  \
    }
// "" to shutdown buggy highlighting in Nano

#endif
