#include <utility>

#include "data_structures.h"

using std::vector;
using std::pair;
using std::make_pair;

/*************** vector_bag ***************/

typedef vector< pair< int, int > > vpii;

void merge_sum(vpii const &lv, vpii const &rv, vpii &merge, bool minus) {
    int li = 0, ri = 0;
    while (li < (int) lv.size() && ri < (int) rv.size()) {
        pair< int, int > const &lp = lv[li];
        pair< int, int > const &rp = rv[ri];
        if (lp.first < rp.first) {
            merge.push_back(lp);
            ++li;
        } else if (lp.first > rp.first) {
            merge.push_back(rp);
            if (minus) {
                pair< int, int > &last = merge.back();
                last.second = -last.second;
            }
            ++ri;
        } else {
            int sum = lp.second + (minus ? -rp.second : rp.second);
            if (sum != 0) {
                merge.push_back(make_pair(lp.first, sum));
            }
            ++li;
            ++ri;
        }
    }
    while (li < (int) lv.size()) {
        merge.push_back(lv[li++]);
    }
    while (ri < (int) rv.size()) {
        merge.push_back(rv[ri++]);
        if (minus) {
            pair< int, int > &last = merge.back();
            last.second = -last.second;
        }
    }
}

vector_bag::vector_bag() : contents() {}

vector_bag::vector_bag(int single_value) : contents(1, make_pair(single_value, 1)) {}

vector_bag::vector_bag(vector_bag const &that) : contents(that.contents) {}

vector_bag &vector_bag::operator = (vector_bag const &that) {
    contents = that.contents;
    return *this;
}

vector_bag &vector_bag::operator += (vector_bag const &that) {
    vector< pair< int, int > > merge;
    merge_sum(contents, that.contents, merge, false);
    contents.swap(merge);
    return *this;
}

vector_bag &vector_bag::operator -= (vector_bag const &that) {
    vector< pair< int, int > > merge;
    merge_sum(contents, that.contents, merge, true);
    contents.swap(merge);
    return *this;
}

vector_bag vector_bag::operator + (vector_bag const &that) const {
    vector_bag rv = *this;
    rv += that;
    return rv;
}

vector_bag vector_bag::operator - (vector_bag const &that) const {
    vector_bag rv = *this;
    rv -= that;
    return rv;
}

vector_bag vector_bag::operator - () const {
    vector_bag rv = *this;
    for (int i = 0, i_max = (int) rv.contents.size(); i < i_max; ++i) {
        rv.contents[i].second = -rv.contents[i].second;
    }
    return rv;
}

int vector_bag::compare_to(vector_bag const &that) const {
    vector< pair< int, int > > const &lv = contents;
    vector< pair< int, int > > const &rv = that.contents;
    int i = 0;
    while (i < (int) lv.size() && i < (int) rv.size()) {
        pair< int, int > const &lp = lv[i];
        pair< int, int > const &rp = rv[i];
        if (lp != rp) return lp < rp ? -1 : 1;
        ++i;
    }
    if (i != (int) lv.size()) return 1;
    if (i != (int) rv.size()) return -1;
    return 0;
}

void vector_bag::normalize() {
    if (contents.size() > 0 && contents[0].second < 0) {
        for (int i = 0, i_max = (int) contents.size(); i < i_max; ++i) {
            contents[i].second = -contents[i].second;
        }
    }
}

bool vector_bag::operator == (vector_bag const &that) const {
    return contents == that.contents;
}

bool vector_bag::operator != (vector_bag const &that) const {
    return contents != that.contents;
}

bool vector_bag::operator < (vector_bag const &that) const {
    return compare_to(that) < 0;
}

bool vector_bag::operator <= (vector_bag const &that) const {
    return compare_to(that) <= 0;
}

bool vector_bag::operator > (vector_bag const &that) const {
    return compare_to(that) > 0;
}

bool vector_bag::operator >= (vector_bag const &that) const {
    return compare_to(that) >= 0;
}

vector< pair< int, int > > const &vector_bag::data() const {
    return contents;
}
