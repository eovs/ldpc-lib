#define _CRT_SECURE_NO_WARNINGS // get rid of MS Visual Studio pesters about fscanf_s

#include "commons_portable.h"

// The independent part, here for convenience

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <map>

std::string format_to_string(char const *format, ...) {
    va_list arglist_1;
    va_start(arglist_1, format);
    char space[] = {0};
    int the_size = vsnprintf(space, 0, format, arglist_1);
    va_end(arglist_1);

    if (the_size == -1) {
        // This is a Visual Studio crap. They return -1 on truncation. What a mess :-E
        // We'll go by exponentiation; just O(log n) queries and we are done       :-E
        int xsize = (int)(strlen(format) + 1);
        do {
            char *placeholder = new char[xsize];
            va_list arglist_X;
            va_start(arglist_X, format);
            the_size = vsnprintf(placeholder, xsize, format, arglist_X);
            va_end(arglist_X);
            delete placeholder;
            xsize *= 2;
        } while (the_size == -1);
    }

    va_list arglist_2;
    va_start(arglist_2, format);
    char *real_chars = new char[the_size + 1];
    vsnprintf(real_chars, the_size + 1, format, arglist_2);
    va_end(arglist_2);
    std::string rv = real_chars;
    delete real_chars;
    return rv;
}

void skip_until_newline(FILE *handle) {
    char ch;
    while (fscanf(handle, "%c", &ch) == 1 && ch != 0x0d && ch != 0x0a);
    // Well, not consuming a second symbol from a two-character EOL
    // is safe because we read only numbers, and scanf on numbers
    // just consumes the white space.
}

int initial_random_seed = 1;	//2011;	//0;

line_eater_type line_eater;

std::istream &operator >> (std::istream &stream, line_eater_type const &) {
    char ch;
    while (stream.get(ch) && ch != '\n');
    if (!stream) {
        return stream;
    }
    if (stream.get(ch) && ch != '\r') {
        stream.unget();
    }
    return stream;
}

// Console handing: Platform-independent part

std::map< char, std::pair< std::string, std::string > > console_messages;
std::map< char, std::pair< std::string, exception_wrapper* > > console_exceptions;

exception_wrapper::~exception_wrapper() {}

void process_console_message(char symbol, std::string const &explanation, std::string const &message) {
    if (console_messages.count(symbol) != 0) {
        die("[console_message_hook] The symbol %c is already hooked: %s", symbol, console_messages[symbol].first.c_str());
    }
    if (console_exceptions.count(symbol) != 0) {
        die("[console_message_hook] The symbol %c is already hooked: %s", symbol, console_exceptions[symbol].first.c_str());
    }
    console_messages.insert(std::make_pair(symbol, std::make_pair(explanation, message)));
}

void console_digest_exception_wrapper(char symbol, std::string const &explanation, exception_wrapper *wrapper) {
    if (console_messages.count(symbol) != 0) {
        delete wrapper;
        die("[console_exception_hook] The symbol %c is already hooked: %s", symbol, console_messages[symbol].first.c_str());
    }
    if (console_exceptions.count(symbol) != 0) {
        delete wrapper;
        die("[console_exception_hook] The symbol %c is already hooked: %s", symbol, console_exceptions[symbol].first.c_str());
    }
    console_exceptions.insert(std::make_pair(symbol, std::make_pair(explanation, wrapper)));
}

console_message_hook::console_message_hook(char symbol, std::string const &explanation, std::string const &message) : symbol(symbol) {
    process_console_message(symbol, explanation, message);
}
console_message_hook::~console_message_hook() {
    console_messages.erase(symbol);
}
console_exception_hook::~console_exception_hook() {
    exception_wrapper *ptr = console_exceptions[symbol].second;
    console_exceptions.erase(symbol);
    delete ptr;
}

void console_check_hooks() {
    if (is_keyboard_hit()) {
        char ch = get_character();
        if (ch == '?') {
            printf("=========================================================\n");
            printf("  Current console hooks installed:\n");
            printf("  - '?': prints this message\n");
            for (std::map< char, std::pair< std::string, exception_wrapper* > >::iterator it = console_exceptions.begin(); it != console_exceptions.end(); ++it) {
                printf("  - '%c': %s\n", it->first, it->second.first.c_str());
            }
            for (std::map< char, std::pair< std::string, std::string > >::iterator it = console_messages.begin(); it != console_messages.end(); ++it) {
                printf("  - '%c': %s\n", it->first, it->second.first.c_str());
            }
            printf("=========================================================\n");
        } else if (console_exceptions.count(ch) != 0) {
            console_exceptions[ch].second->execute();
        } else if (console_messages.count(ch) != 0) {
            printf("%s\n", console_messages[ch].second.c_str());
        }
    }
}

// Compiler-dependent part: Random support

#if __cplusplus >= 201103L

    #define CPH_RANDOM_MESSAGE "[commons_portable.cpp] C++11 is detected, using <random> to implement next_random_* calls."

    #include <random>

    std::mt19937 generator(-1);

    bool random_initialized = false;

    void reset_random() {
        random_initialized = false;
    }

    void ensure_random_is_initialized() {
        if (!random_initialized) {
            if (initial_random_seed == 0) {
                std::random_device device;
                initial_random_seed = (int) (device());
                if (initial_random_seed == 0) {
                    initial_random_seed = 1;
                }
            }
            generator = std::mt19937(initial_random_seed);
            random_initialized = true;
        }
    }

    int next_random_int(int minInclusive, int maxExclusive) {
        ensure_random_is_initialized();
        std::uniform_int_distribution< int > dist(minInclusive, maxExclusive - 1);
        return dist(generator);
    }

    double next_random_01() {
        ensure_random_is_initialized();
        std::uniform_real_distribution< double > dist;
        return dist(generator);
    }

    double next_random_gaussian() {
        ensure_random_is_initialized();
        std::normal_distribution< double > dist;
        return dist(generator);
    }

    [[noreturn]]
    void die(char const *format, ...) {
        va_list arglist;
        va_start(arglist, format);
        vfprintf(stderr, format, arglist);
        va_end(arglist);
        fputs("\n", stderr);
        exit(1);
        throw 0;
    }

    FILE *fopen_or_die(char const *filename, char const *mode, char const *error_message) {
        FILE *result = fopen(filename, mode);
        if (result == NULL) {
            die(error_message);
        } else {
            return result;
        }
    }

#else

    #define CPH_RANDOM_MESSAGE "[commons_portable.cpp] C++11 is not detected, using rand() to implement next_random_* calls."

    #include <math.h>
    #include <time.h>

    bool random_initialized = false;

    void reset_random() {
        random_initialized = false;
    }

    void ensure_random_is_initialized() {
        if (!random_initialized) {
            if (initial_random_seed == 0) {
                initial_random_seed = (int) (time(NULL));
                if (initial_random_seed == 0) {
                    initial_random_seed = 1;
                }
            }
            srand(initial_random_seed);
            random_initialized = true;
        }
    }

    int next_random_int(int minInclusive, int maxExclusive) {
        ensure_random_is_initialized();
        return rand() % (maxExclusive - minInclusive) + minInclusive;
    }

    double next_random_01() {
        ensure_random_is_initialized();

        int bits_per_random = 0;
        for (int v = RAND_MAX; v > 0; v >>= 1, ++bits_per_random);

        // Collecting enough digits...
        int required_bits = 54;
        int collected_bits = 0;
        long long collected = 0;
        while (collected_bits < required_bits) {
            int value = rand();
            int required_diff = required_bits - collected_bits;
            if (required_diff >= bits_per_random) {
                collected = (collected << bits_per_random) + value;
                collected_bits += bits_per_random;
            } else {
                collected = (collected << required_diff) + (value >> (bits_per_random - required_diff));
                collected_bits += required_diff;
            }
        }
        if (collected < 0 || collected >= (1LL << required_bits)) {
            die("next_random_01 [C++98] does not work properly");
        }
        return (double) (collected) / (double) (1LL << required_bits);
    }

    double next_random_gaussian() {
        ensure_random_is_initialized();
        static bool has_next_next = false;
        static double next_next = 0;

        if (has_next_next) {
            has_next_next = false;
            return next_next;
        } else {
            double u, v, w;
            do {
                u = 2.0 * next_random_01() - 1.0;
                v = 2.0 * next_random_01() - 1.0;
                w = u * u + v * v;
            } while (w >= 1.0);
            w = sqrt(-2 * log(w) / w);
            next_next = v * w;
            has_next_next = true;
            return u * w;
        }
    }

    void die(char const *format, ...) {
        va_list arglist;
        va_start(arglist, format);
        vfprintf(stderr, format, arglist);
        va_end(arglist);
        fputs("\n", stderr);
        exit(1);
        throw 0;
    }

    FILE *fopen_or_die(char const *filename, char const *mode, char const *error_message) {
        FILE *result = fopen(filename, mode);
        if (result == NULL) {
            die(error_message);
            throw 0;
        } else {
            return result;
        }
    }

#endif


// Keyboard-hit: Platform-dependent part

#ifdef _WIN32

    #define CPH_KBHIT_MESSAGE "[commons_portable.cpp] Windows detected, using <conio.h> to implement keyboard-hits."

    #include <conio.h>

    bool is_keyboard_hit() {
        return _kbhit() != 0;
    }

    char get_character() {
        return _getch();
    }

    std::string get_native_path(std::string path) {
        for (std::string::size_type i = 0; i < path.size(); ++i) {
            if (path[i] == '/') {
                path[i] = '\\';
            }
        }
        return path;
    }

#elif __linux__

    #define CPH_KBHIT_MESSAGE "[commons_portable.cpp] Linux detected, using a termios/ioctl workaround to implement keyboard-hits."

    #include <stropts.h>
    #include <sys/ioctl.h>
    #include <sys/select.h>
    #include <termios.h>

    #define STDIN_HANDLE 0

    bool is_keyboard_hit() {
        static bool initialized = false;

        if (! initialized) {
            // Use termios to turn off line buffering
            termios term;
            tcgetattr(STDIN_HANDLE, &term);
            term.c_lflag &= ~ICANON;
            tcsetattr(STDIN_HANDLE, TCSANOW, &term);
            setbuf(stdin, NULL);
            initialized = true;
        }

        int bytesWaiting;
        ioctl(STDIN_HANDLE, FIONREAD, &bytesWaiting);
        return bytesWaiting;
    }

    char get_character() {
        return getchar();
    }

    std::string get_native_path(std::string path) {
        return path;
    }

#else
    #error "Cannot determine the platform (neither _WIN32 nor __linux__ were defined)"
#endif

void commons_portable_print_info() {
    printf("%s\n", CPH_RANDOM_MESSAGE);
    printf("%s\n", CPH_KBHIT_MESSAGE);
    if (initial_random_seed == 0) {
        printf("[commons_portable.cpp] initial_random_seed = 0, will init at random.\n");
    } else {
        printf("[commons_portable.cpp] initial_random_seed = %d, will use it as the seed.\n", initial_random_seed);
    }
}
