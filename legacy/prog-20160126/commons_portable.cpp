#include "commons_portable.h"

// The independent part, here for convenience

#include <stdio.h>
#include <stdlib.h>

FILE *fopen_or_die(char const *filename, char const *mode, char const *error_message) {
    FILE *result = fopen(filename, mode);
    if (result == NULL) {
        fputs(error_message, stderr);
        fputs("\n", stderr);
        exit(1);
    } else {
        return result;
    }
}

void skip_until_newline(FILE *handle) {
    char ch;
    while (fscanf(handle, "%c", &ch) == 1 && ch != 0x0d && ch != 0x0a);
    // Well, not consuming a second symbol from a two-character EOL
    // is safe because we read only numbers, and scanf on numbers
    // just consumes the white space.
}

int initial_random_seed = 0;

// Compiler-dependent part: Random support

#if __cplusplus >= 201103L

    #define CPH_RANDOM_MESSAGE "[commons_portable.cpp] C++11 is detected, using <random> to implement next_random_* calls."

    #include <random>

    std::mt19937 generator(-1);

    void ensure_random_is_initialized() {
        static bool random_initialized = false;
        if (!random_initialized) {
            if (initial_random_seed == 0) {
                std::random_device device;
                generator = std::mt19937(device());
            } else {
                generator = std::mt19937(initial_random_seed);
            }
            random_initialized = true;
        }
    }

    int next_random_int(int minInclusive, int maxExclusive) {
        std::uniform_int_distribution< int > dist(minInclusive, maxExclusive - 1);
        return dist(generator);
    }

#else

    #define CPH_RANDOM_MESSAGE "[commons_portable.cpp] C++11 is not detected, using rand() to implement next_random_* calls."

    #include <time.h>

    void ensure_random_is_initialized() {
        static bool random_initialized = false;
        if (!random_initialized) {
            if (initial_random_seed == 0) {
                srand((unsigned) (time(NULL)));
            } else {
                srand(initial_random_seed);
            }
            random_initialized = true;
        }
    }

    int next_random_int(int minInclusive, int maxExclusive) {
        ensure_random_is_initialized();
        return rand() % (maxExclusive - minInclusive) + minInclusive;
    }
#endif


// Keyboard-hit: Platform-dependent part

#ifdef _WIN32

    #define CPH_KBHIT_MESSAGE "[commons_portable.cpp] Windows detected, using <conio.h> to implement keyboard-hits."

    #include <conio.h>

    bool is_keyboard_hit() {
        return _kbhit();
    }

    char get_character() {
        return _getch();
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
