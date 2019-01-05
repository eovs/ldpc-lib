.PHONY: all tests s
.RECIPEPREFIX = >

all: main

tests: main
> ./main tests

clean:
> rm -f *.o main

CC=g++

SOURCES_CLEAN=data_structures.cpp equations.cpp commons_portable.cpp combinatorics.cpp code_generation.cpp base_matrix.cpp settings.cpp main_good_code_search.cpp main_unit_tests.cpp main_search_ggp.cpp main_simulation.cpp main.cpp
SOURCES_DIRTY=decoders.cpp bp_simulation.cpp direct_inverse_perm.cpp QAM_demodulator.cpp QAM_modulator.cpp

OBJECTS_CLEAN=$(SOURCES_CLEAN:.cpp=.o)
OBJECTS_DIRTY=$(SOURCES_DIRTY:.cpp=.o)

CFLAGS_CLEAN=-O3 -Wall -Wextra -ggdb
CFLAGS_DIRTY=$(CFLAGS_CLEAN) -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function -DSKIP_MEX

$(OBJECTS_CLEAN): CFLAGS := $(CFLAGS_CLEAN)
$(OBJECTS_DIRTY): CFLAGS := $(CFLAGS_DIRTY)

main: $(OBJECTS_CLEAN) $(OBJECTS_DIRTY)
> $(CC) $(CFLAGS_CLEAN) -o main $(OBJECTS_CLEAN) $(OBJECTS_DIRTY)

%.o: %.cpp
> $(CC) $(CFLAGS) -c $<


# SEARCH targets

single: main
> ./main search files/input32_16_single.jsonx   files/out_single.jsonx

double: main
> ./main search files/input32_16_double.jsonx   files/out_double.jsonx

qualcomm: main
> ./main search files/input32_16_qualcomm.jsonx files/out_qualcomm.jsonx

