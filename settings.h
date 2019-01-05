#ifndef __SETTINGS_H
#define __SETTINGS_H

#include <exception>
#include <stdexcept>
#include <istream>
#include <ostream>
#include <map>
#include <sstream>
#include <string>

#include "commons_portable.h"
#include "data_structures.h"

/*
 * The syntax rules are:
 * object = number                                                                    // type = STRING
 *        | string-in-quotes                                                          // type = STRING
 *        | 'array' '{' object* '}'                                                   // type = ARRAY
 *        | 'array' '@' string-in-quotes                                              // type = ARRAY // parses contents from a file
 *        | 'matrix' '(' number ',' number ')' '{' object* '}' // row 0, row 1, ...;  // type = MATRIX
 *        | 'sparse matrix' '(' number ',' number ')' '{' [number number object]* '}' // type = MATRIX
 *        | '{' [identifier '=' object]* '}'                                          // type = RECORD
 *        | '@' string-in-quotes // for an object in the file, filename is relative , resolves immediately
 */

class settings {
public:
    enum object_type { RECORD, STRING, ARRAY, MATRIX, OBJECT_REF, ARRAY_REF };
private:
    object_type my_type;
    std::string my_path;
    std::string my_path_prefix;

    std::string                       *data_string;
    std::vector< settings >           *data_vector;
    matrix< settings >                *data_matrix;
    std::map< std::string, settings > *data_record;

    void erase_contents();
    void to_stream_impl(std::ostream &stream, int indent, bool newline) const;
public:
    settings();
    settings(settings const &that);
    settings(std::istream &stream,
             std::string const &current_directory = ".",
             std::string const &path = "",
             std::string const &path_prefix = "");
    settings &operator = (settings const &that);
    ~settings();

    void swap(settings &that);

    static settings from_file(std::string const &filename);
    static settings from_stream(std::istream &stream);

    void to_stream(std::ostream &stream) const;
    void to_file(std::string const &filename, bool append = false) const;

    std::string get_full_path() const;
    std::string get_path() const;
    object_type get_type() const;

    bool can_select(std::string const &path) const;
    settings const &select(std::string const &path) const;
    settings &open(std::string const &path);

    template<typename result_t> void cast_to(result_t &target) const {
        if (my_type == STRING) {
            std::istringstream iss(*data_string);
            iss >> target;
        } else {
            die("Cannot convert a non-string object '%s' into what you ask!", get_full_path().c_str());
        }
    }

    void cast_to(std::string &target) const;
    void cast_to(settings &target) const;

    template<typename element_t> void cast_to(std::vector<element_t> &target) const {
        if (my_type == ARRAY) {
            std::vector<element_t> rv(data_vector->size());
            for (size_t i = 0; i < data_vector->size(); ++i) {
                data_vector->at(i).cast_to(rv[i]);
            }
            target = rv;
        } else {
            die("Cannot convert a non-array object '%s' into a vector!", get_full_path().c_str());
        }
    }

    template<typename element_t> void cast_to(matrix<element_t> &target) const {
        if (my_type == MATRIX) {
            matrix<element_t> rv(data_matrix->n_rows(), data_matrix->n_cols());
            for (int row = 0; row < data_matrix->n_rows(); ++row) {
                for (int col = 0; col < data_matrix->n_cols(); ++col) {
                    settings const &current = (*data_matrix)(row, col);
                    if (current.get_path().length() == 0) {
                        // Sparse matrices are the only case for this
                        rv(row, col) = element_t();
                    } else {
                        current.cast_to(rv(row, col));
                    }
                }
            }
            target = rv;
        } else {
            die("Cannot convert a non-matrix object '%s' into a matrix!", get_full_path().c_str());
        }
    }

    void set_object_ref(std::string const &path);
    void set_array_ref(std::string const &path);

    void set(settings const &value);

    template<typename value_t> void set(value_t const &value) {
        erase_contents();
        my_type = STRING;
        std::ostringstream oss;
        oss << value;
        data_string = new std::string(oss.str());
    }

    template<typename value_t> void set(std::vector<value_t> const &value) {
        erase_contents();
        my_type = ARRAY;
        data_vector = new std::vector<settings>(value.size());
        for (size_t i = 0, i_max = value.size(); i < i_max; ++i) {
            // TODO: my_path and my_path_prefix
            data_vector->at(i).set(value[i]);
        }
    }

    template<typename value_t> void set(matrix<value_t> const &value) {
        erase_contents();
        my_type = MATRIX;
        data_matrix = new matrix<settings>(value.n_rows(), value.n_cols());
        for (int r = 0, r_max = value.n_rows(); r < r_max; ++r) {
            for (int c = 0, c_max = value.n_cols(); c < c_max; ++c) {
                // TODO: my_path and my_path_prefix
                (*data_matrix)(r, c).set(value(r, c));
            }
        }
    }

    template<typename value_t> void append(value_t const &value) {
        if (my_type == ARRAY) {
            // TODO: my_path and my_path_prefix
            settings nv;
            nv.set(value);
            data_vector->push_back(nv);
        } else {
            die("Cannot append to a non-array object '%s'!", get_full_path().c_str());
        }
    }
};

struct interval_exception : std::invalid_argument {
    inline interval_exception(std::string const &msg) : std::invalid_argument(msg) {}
};

template<typename output_iterator>
void read_intervals(std::string const &data, output_iterator out) {
    int number_1 = 0, number_2 = 0;
    int state = 0;
        // 0 - before parsing a common ("first") number
        // 1 - inside first number
        // 2 - just finished parsing first number
        // 3 - hyphen found = before parsing second number
        // 4 - inside second number

    for (unsigned i = 0, i_max = (unsigned int)data.length(); i <= i_max; ++i) {
        char ch = i == i_max ? ' ' : data[i];
        const int SPC = 0, DIG = 1, HYP = 2, ETC = 3;
        int category = ch <= ' ' ? SPC : ('0' <= ch && ch <= '9' ? DIG : (ch == '-' ? HYP : ETC));
        if (category == 3) {
            throw interval_exception(format_to_string("String '%s', index %u: unexpected symbol '%c'", data.c_str(), i + 1, ch));
        }
        switch (state) {
            case 0: switch (category) {
                case SPC:                      state = 0; break;
                case DIG: number_1 = ch - '0'; state = 1; break;
                case HYP: throw interval_exception(format_to_string("String '%s', index %u: unexpected '-'", data.c_str(), i + 1));
            } break;
            case 1: switch (category) {
                case SPC: *out++ = number_1;                   state = 2; break;
                case DIG: number_1 = 10 * number_1 + ch - '0'; state = 1; break;
                case HYP: *out++ = number_1;                   state = 3; break;
            } break;
            case 2: switch (category) {
                case SPC:                      state = 2; break;
                case DIG: number_1 = ch - '0'; state = 1; break;
                case HYP:                      state = 3; break;
            } break;
            case 3: switch (category) {
                case SPC:                      state = 3; break;
                case DIG: number_2 = ch - '0'; state = 4; break;
                case HYP: throw interval_exception(format_to_string("String '%s', index %u: unexpected '-'", data.c_str(), i + 1));
            } break;
            case 4: switch (category) {
                case SPC: {
                    if (number_2 < number_1) {
                        throw interval_exception(format_to_string("String '%s': an empty interval %d-%d requested", data.c_str(), number_1, number_2));
                    }
                    for (int j = number_1 + 1; j <= number_2; *out++ = j++);
                    state = 0;
                    break;
                }
                case DIG: number_2 = 10 * number_2 + ch - '0';                       state = 4; break;
                case HYP: throw interval_exception(format_to_string("String '%s', index %u: unexpected '-'", data.c_str(), i + 1));
            } break;
        }
    }
    if (state != 0 && state != 2) {
        throw interval_exception(format_to_string("String '%s': the last interval is not finished", data.c_str()));
    }
}

#endif
