#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <utility>

#include "commons_portable.h"
#include "main_modules.h"

using std::make_pair;
using std::map;
using std::string;

typedef int (*entry_type) (int, char *[]);
typedef map< string, entry_type >::const_iterator map_const_itr;

map< string, entry_type > entry_points;

void print_usage(char const *prog_name) {
    fprintf(stderr, "Usage: %s <entry-point> args...\n", prog_name);
    fprintf(stderr, "    where <entry-point> is one of the following:\n");
    for (map_const_itr it = entry_points.begin(); it != entry_points.end(); ++it) {
        fprintf(stderr, "        %s\n", it->first.c_str());
    }
}

int main(int argc, char *argv[]) {
    entry_points.insert(make_pair("search", &main_good_code_search));
    entry_points.insert(make_pair("tests", &main_unit_tests));
    entry_points.insert(make_pair("ggp", &main_search_ggp));
	entry_points.insert(make_pair("simulation", &main_simulation));

    if (argc <= 1) {
        fprintf(stderr, "Error: no entry point specified\n");
        print_usage(argv[0]);
        return 1;
    } else if (entry_points.find(argv[1]) == entry_points.end()) {
        fprintf(stderr, "Error: unknown entry point '%s'\n", argv[1]);
        print_usage(argv[0]);
        return 1;
    } else {
        return (entry_points[argv[1]])(argc - 1, argv + 1);
    }
}
