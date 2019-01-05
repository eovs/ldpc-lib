#include "settings.h"

#include <ios>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <utility>

using std::map;
using std::ios_base;
using std::istream;
using std::ifstream;
using std::istringstream;
using std::make_pair;
using std::pair;
using std::string;
using std::ostream;
using std::ofstream;
using std::ostringstream;
using std::vector;

settings settings::from_stream(istream &stream) {
    return settings(stream);
}

settings settings::from_file(string const &filename) {
    string native = get_native_path(filename);
    ifstream file(native.c_str());
    if (!file) {
        die("File not found: %s", native.c_str());
    }
    string dir;
    size_t last_slash = filename.rfind('/');
    if (last_slash == string::npos) {
        last_slash = filename.rfind('\\');
    }
    if (last_slash == string::npos) {
        dir = ".";
    } else {
        dir = filename.substr(0, last_slash);
    }
    return settings(file, dir, "", filename);
}

void settings::erase_contents() {
    if (data_string) delete data_string;
    if (data_vector) delete data_vector;
    if (data_matrix) delete data_matrix;
    if (data_record) delete data_record;

    data_string = NULL;
    data_vector = NULL;
    data_matrix = NULL;
    data_record = NULL;
}

settings::~settings() {
    erase_contents();
}

void settings::swap(settings &that) {
    std::swap(my_type, that.my_type);
    std::swap(my_path, that.my_path);
    std::swap(my_path_prefix, that.my_path_prefix);
    std::swap(data_string, that.data_string);
    std::swap(data_vector, that.data_vector);
    std::swap(data_matrix, that.data_matrix);
    std::swap(data_record, that.data_record);
}

settings &settings::operator = (settings const &that) {
    my_type = that.my_type;
    my_path = that.my_path;
    my_path_prefix = that.my_path_prefix;
    data_string = that.data_string ? new string(*that.data_string) : 0;
    data_vector = that.data_vector ? new vector<settings>(*that.data_vector) : 0;
    data_matrix = that.data_matrix ? new matrix<settings>(*that.data_matrix) : 0;
    data_record = that.data_record ? new map<string, settings>(*that.data_record) : 0;
    return *this;
}

// The default constructor. This will create an empty object, i.e. "{}"
settings::settings()
    : my_type(settings::RECORD)
    , my_path("")
    , my_path_prefix("")
    , data_string(0)
    , data_vector(0)
    , data_matrix(0)
    , data_record(new map<string, settings>())
{}

// The copy constructor. Nothing special except data values are cloned.
settings::settings(settings const &that)
    : my_type(that.my_type)
    , my_path(that.my_path)
    , my_path_prefix(that.my_path_prefix)
    , data_string(that.data_string ? new string(*that.data_string) : 0)
    , data_vector(that.data_vector ? new vector<settings>(*that.data_vector) : 0)
    , data_matrix(that.data_matrix ? new matrix<settings>(*that.data_matrix) : 0)
    , data_record(that.data_record ? new map<string, settings>(*that.data_record) : 0)
{}

// Various private stuff around
// --------------------------------------------------------------------
void skip_ws(istream &stream) {
    while (true) {
        while (!stream.eof() && stream.peek() <= ' ') {
            stream.get();
        }
        if (!stream.eof() && stream.peek() == '/') {
            while (!stream.eof() && stream.get() != '\n');
        } else {
            break;
        }
    }
}

void expect(istream &stream, string const &path, string const &path_prefix, string const &text) {
    string accumulator;
    for (int i = 0, i_max = (int) (text.length()); i < i_max; ++i) {
        if (stream.eof()) {
            die("Unexpected end of file where '%s' expected (current context: '%s%s')", text.c_str(), path_prefix.c_str(), path.c_str());
        }
        char ch = (char) stream.get();
        accumulator.push_back(ch);
        if (ch != text[i]) {
            char const *dots = text.length() == 1 ? "" : "...";
            die("Expected '%s' found '%s%s' (current context: '%s%s')", text.c_str(), accumulator.c_str(), dots, path_prefix.c_str(), path.c_str());
        }
    }
}

string digest_quoted_string(istream &stream, string const &path, string const &path_prefix) {
    string buf;
    expect(stream, path, path_prefix, "\"");
    while (!stream.eof() && stream.peek() != '"') {
        buf += (char) stream.get();
    }
    expect(stream, path, path_prefix, "\"");
    return buf;
}

pair<string, string> resolve_path(string const &relative_path, string const &current_directory) {
    string new_file = current_directory + "/" + relative_path;
    string new_dir = new_file;
    size_t last_slash = new_dir.rfind('/');
    new_dir.resize(last_slash);
    string native = get_native_path(new_file);
    return make_pair(native, new_dir);
}

bool is_identifier_part(char v) {
    return v == '_' || (v >= 'a' && v <= 'z') || (v >= 'A' && v <= 'Z') || (v >= '0' && v <= '9') || v == '-';
}

string digest_identifier(istream &stream) {
    string buf;
    while (!stream.eof() && is_identifier_part(stream.peek())) {
        buf += (char) stream.get();
    }
    return buf;
}

bool is_number_part(char v) {
    return v == '.' || (v >= '0' && v <= '9') || v == '-' || v == 'e' || v == 'E';
}

string digest_number(istream &stream) {
    string buf;
    while (!stream.eof() && is_number_part(stream.peek())) {
        buf += (char) stream.get();
    }
    return buf;
}
// --------------------------------------------------------------------

// The constructor which actually reads and parses the configuration.
settings::settings(istream &stream, string const &current_directory, string const &path, string const &path_prefix)
    : my_path(path)
    , my_path_prefix(path_prefix)
    , data_string(0)
    , data_vector(0)
    , data_matrix(0)
    , data_record(0)
{
    if (!stream) {
        die("Stream is invalid. Current directory = '%s', path = '%s', path_prefix = '%s'",
            current_directory.c_str(), path.c_str(), path_prefix.c_str());
    }
    skip_ws(stream);
    if (stream.eof()) {
        // behave as if the default constructor is called
        my_type = RECORD;
        data_record = new map<string, settings>();
    } else {
        switch (stream.peek()) {
            case '"': { // string in quotes
                my_type = STRING;
                data_string = new string(digest_quoted_string(stream, path, path_prefix));
                break;
            }
            case 'a': { // array
                my_type = ARRAY;
                data_vector = new vector<settings>();

                expect(stream, path, path_prefix, "array");
                skip_ws(stream);
                if (!stream.eof() && stream.peek() == '{') {
                    // explicit specified contents
                    expect(stream, path, path_prefix, "{");
                    skip_ws(stream);
                    while (!stream.eof() && stream.peek() != '}') {
                        ostringstream oss;
                        oss << path << "/" << data_vector->size();
                        data_vector->push_back(settings(stream, current_directory, oss.str(), path_prefix));
                        skip_ws(stream);
                    }
                    expect(stream, path, path_prefix, "}");
                } else {
                    // contents to be read from a file
                    expect(stream, path, path_prefix, "@");
                    skip_ws(stream);
                    pair<string, string> paths = resolve_path(digest_quoted_string(stream, path, path_prefix), current_directory);
                    ifstream file(paths.first.c_str());
                    skip_ws(file);
                    while (!file.eof()) {
                        ostringstream oss;
                        oss << path << "/" << data_vector->size();
                        data_vector->push_back(settings(file, paths.second, oss.str(), path_prefix));
                        skip_ws(file);
                    }
                }
                break;
            }
            case 'm': { // matrix
                my_type = MATRIX;
                expect(stream, path, path_prefix, "matrix");
                skip_ws(stream);
                expect(stream, path, path_prefix, "(");
                int n_rows, n_cols;
                settings(stream, current_directory, path + "/#rows", path_prefix).cast_to(n_rows);
                settings(stream, current_directory, path + "/#cols", path_prefix).cast_to(n_cols);
                skip_ws(stream);
                expect(stream, path, path_prefix, ")");

                data_matrix = new matrix<settings>(n_rows, n_cols);

                skip_ws(stream);
                expect(stream, path, path_prefix, "{");
                for (int row = 0; row < n_rows; ++row) {
                    for (int col = 0; col < n_cols; ++col) {
                        ostringstream oss;
                        oss << "/(" << row << "," << col << ")";
                        (*data_matrix)(row, col) = settings(stream, current_directory, path + oss.str(), path_prefix);
                    }
                }
                skip_ws(stream);
                expect(stream, path, path_prefix, "}");
                break;
            }
            case 's': { // sparse matrix
                my_type = MATRIX;
                expect(stream, path, path_prefix, "sparse");
                skip_ws(stream);
                expect(stream, path, path_prefix, "matrix");
                skip_ws(stream);
                expect(stream, path, path_prefix, "(");
                int n_rows, n_cols;
                settings(stream, current_directory, path + "/#rows", path_prefix).cast_to(n_rows);
                settings(stream, current_directory, path + "/#cols", path_prefix).cast_to(n_cols);
                skip_ws(stream);
                expect(stream, path, path_prefix, ")");

                data_matrix = new matrix<settings>(n_rows, n_cols);
                skip_ws(stream);
                expect(stream, path, path_prefix, "{");
                skip_ws(stream);
                int count = 0;
                while (!stream.eof() && stream.peek() != '}') {
                    int row, col;
                    ostringstream xr, xc, xv;
                    xr << "/row#" << count;
                    xc << "/col#" << count;
                    xv << "/val#" << count;
                    string pxr = path + xr.str();
                    settings(stream, current_directory, pxr, path_prefix).cast_to(row);
                    if (row < 0 || row >= n_rows) {
                        die("Row index out of bounds: %d, number of rows: %d (current context: '%s%s')", row, n_rows, path_prefix.c_str(), pxr.c_str());
                    }
                    string pxc = path + xc.str();
                    settings(stream, current_directory, pxc, path_prefix).cast_to(col);
                    if (col < 0 || col >= n_cols) {
                        die("Column index out of bounds: %d, number of columns: %d (current context: '%s%s')", col, n_cols, path_prefix.c_str(), pxc.c_str());
                    }
                    string pxv = path + xv.str();
                    (*data_matrix)(row, col) = settings(stream, current_directory, pxv, path_prefix);
                    skip_ws(stream);
                    ++count;
                }
                expect(stream, path, path_prefix, "}");
                break;
            }
            case '{': { // record
                my_type = RECORD;
                data_record = new map<string, settings>();
                expect(stream, path, path_prefix, "{");
                skip_ws(stream);
                while (!stream.eof() && stream.peek() != '}') {
                    string id = digest_identifier(stream);
                    skip_ws(stream);
                    expect(stream, path, path_prefix, "=");
                    settings field(stream, current_directory, path + "/" + id, path_prefix);
                    skip_ws(stream);
                    data_record->insert(make_pair(id, field));
                }
                expect(stream, path, path_prefix, "}");
                break;
            }
            case '@': { // file reference
                expect(stream, path, path_prefix, "@");
                pair<string, string> paths = resolve_path(digest_quoted_string(stream, path, path_prefix), current_directory);
                ifstream file(paths.first.c_str());
                settings referenced(file, paths.second, path, path_prefix);
                swap(referenced);
                break;
            }
            default: { // number
                my_type = STRING;
                data_string = new string(digest_number(stream));
                break;
            }
        }
    }
}

std::string           settings::get_full_path() const { return my_path_prefix + my_path; }
std::string           settings::get_path()      const { return my_path;                  }
settings::object_type settings::get_type()      const { return my_type;                  }

bool settings::can_select(string const &path) const {
    if (path.length() == 0) {
        return true;
    }
    if (my_type != RECORD) {
        return false;
    }
    size_t slashIndex = path.find('/'); //npos is processed correctly
    string key = path.substr(0, slashIndex);
    string etc = slashIndex == string::npos ? "" : path.substr(slashIndex + 1);
    map<string, settings>::iterator itr = data_record->find(key);
    if (itr != data_record->end()) {
        return itr->second.can_select(etc);
    } else {
        map<string, settings>::iterator def = data_record->find("defaults");
        if (def != data_record->end()) {
            return def->second.can_select(path);
        } else {
            return false;
        }
    }
}

settings const &settings::select(string const &path) const {
    if (path.length() == 0) {
        // 0. Exact hit. Return itself.
        return *this;
    }
    if (my_type != RECORD) {
        die("Object '%s' is not a record, it is impossible to select from it!", get_full_path().c_str());
    }
    size_t slashIndex = path.find('/'); //npos is processed correctly
    string key = path.substr(0, slashIndex);
    string etc = slashIndex == string::npos ? "" : path.substr(slashIndex + 1);
    map<string, settings>::iterator itr = data_record->find(key);
    if (itr != data_record->end()) {
        // 1. The key is found, continue selecting.
        return itr->second.select(etc);
    } else {
        // 2. The key is not found. Looking for 'defaults'
        map<string, settings>::iterator def = data_record->find("defaults");
        if (def != data_record->end()) {
            // 2.1. The 'defaults' object is found.
            if (def->second.can_select(path)) {
                // 2.1.1. Forwarding...
                return def->second.select(path);
            } else {
                // 2.1.2. Defaults does not contain the field. Complaining the right way...
                die("Object '%s' does not contain a field '%s' ('defaults' has also been checked)!", get_full_path().c_str(), key.c_str());
                // never throws
                throw 239;
            }
        } else {
            // 2.2. No 'defaults' object
            die("Object '%s' does not contain a field '%s'!", get_full_path().c_str(), key.c_str());
            // never throws
            throw 239;
        }
    }
}

// Opening. Creates missing records and overwrites defaults
settings &settings::open(string const &path) {
    if (path.length() == 0) {
        // 0. Exact hit. Return itself.
        return *this;
    }
    if (my_type != RECORD) {
        die("Object '%s' is not a record, it is impossible to open it!", get_full_path().c_str());
    }
    size_t slashIndex = path.find('/'); //npos is processed correctly
    string key = path.substr(0, slashIndex);
    string etc = slashIndex == string::npos ? "" : path.substr(slashIndex + 1);
    map<string, settings>::iterator itr = data_record->find(key);
    if (itr != data_record->end()) {
        // 1. The key is found, continue opening.
        return itr->second.open(etc);
    } else {
        // 2. The key is not found. Creating it.
        settings &new_setting = (*data_record)[key];
        new_setting.my_path_prefix = my_path_prefix;
        new_setting.my_path = my_path + "/" + key;
        return new_setting.open(etc);
    }
}

string indentor(int howMuch) {
    return string(howMuch * 2, ' ');
}

void settings::to_stream_impl(std::ostream &stream, int indent, bool newline) const {
    char last = newline ? '\n' : ' ';
    switch (my_type) {
        case RECORD: {
            stream << "{\n";
            map<string, settings>::const_iterator itr = data_record->begin();
            while (itr != data_record->end()) {
                stream << indentor(indent + 1) << itr->first << " = ";
                itr->second.to_stream_impl(stream, indent + 1, true);
                ++itr;
            }
            stream << indentor(indent) << "}" << last;
            break;
        }
        case STRING: {
            bool is_number = true;
            for (size_t i = 0, i_max = data_string->length(); i < i_max; ++i) {
                is_number &= is_number_part((*data_string)[i]);
            }
            if (is_number) {
                stream << *data_string << last;
            } else {
                stream << '"' << *data_string << '"' << last;
            }
            break;
        }
        case ARRAY: {
            stream << "array {";
            bool has_nonstring = false;
            for (size_t i = 0, i_max = data_vector->size(); i < i_max; ++i) {
                has_nonstring |= data_vector->at(i).my_type != STRING;
            }
            stream << (has_nonstring ? "\n" : " ");
            for (size_t i = 0, i_max = data_vector->size(); i < i_max; ++i) {
                if (has_nonstring) {
                    stream << indentor(indent + 1);
                }
                data_vector->at(i).to_stream_impl(stream, indent + 1, has_nonstring);
            }
            if (has_nonstring) {
                stream << indentor(indent);
            }
            stream << "}" << last;
            break;
        }
        case MATRIX: {
            stream << "matrix (" << data_matrix->n_rows() << " " << data_matrix->n_cols() << ") {\n";

            unsigned field_length = 0;
            for (int r = 0, r_max = data_matrix->n_rows(); r < r_max; ++r) {
                for (int c = 0, c_max = data_matrix->n_cols(); c < c_max; ++c) {
                    ostringstream oss;
                    (*data_matrix)(r, c).to_stream_impl(oss, indent + 1, false);
                    // no whitespace here
                    field_length = std::max(field_length, (unsigned) oss.str().length() - 1);
                }
            }
            for (int r = 0, r_max = data_matrix->n_rows(); r < r_max; ++r) {
                stream << indentor(indent + 1);
                for (int c = 0, c_max = data_matrix->n_cols(); c < c_max; ++c) {
                    stream << std::setw(field_length);
                    (*data_matrix)(r, c).to_stream_impl(stream, indent + 1, c + 1 == c_max);
                }
            }
            stream << indentor(indent) << "}" << last;
            break;
        }
        case OBJECT_REF: {
            stream << "@\"" << *data_string << '"' << last;
            break;
        }
        case ARRAY_REF: {
            stream << "array @\"" << *data_string << '"' << last;
            break;
        }
    }
}


void settings::to_stream(ostream &stream) const {
    to_stream_impl(stream, 0, true);
}

void settings::to_file(string const &filename, bool append) const {
    string native = get_native_path(filename);
    ofstream out(native.c_str(), append ? ios_base::app : ios_base::out);
    to_stream(out);
}

void settings::set(settings const &value) {
    switch (value.my_type) {
        case RECORD: {
            erase_contents();
            my_type = RECORD;
            data_record = new std::map<std::string, settings>();
            std::map<std::string, settings>::const_iterator itr = value.data_record->begin();
            while (itr != value.data_record->end()) {
                open(itr->first).set(itr->second);
                ++itr;
            }
            break;
        }
        case STRING: {
            set(*value.data_string);
            break;
        }
        case ARRAY: {
            set(*value.data_vector);
            break;
        }
        case MATRIX: {
            set(*value.data_matrix);
            break;
        }
        case OBJECT_REF: {
            set_object_ref(*value.data_string);
            break;
        }
        case ARRAY_REF: {
            set_array_ref(*value.data_string);
            break;
        }
    }
}

void settings::cast_to(string &target) const {
    if (my_type == STRING) {
        target = *data_string;
    } else {
        die("Cannot convert a non-string object '%s' into a string!", get_full_path().c_str());
    }
}

void settings::cast_to(settings &target) const {
    target = *this;
}

void settings::set_object_ref(string const &path) {
    erase_contents();
    data_string = new string(path);
    my_type = OBJECT_REF;
}

void settings::set_array_ref(string const &path) {
    erase_contents();
    data_string = new string(path);
    my_type = ARRAY_REF;
}

