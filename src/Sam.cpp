#include "Sam.hpp"

#include <stdexcept>
#include <iostream>
#include <istream>
#include <limits>

using std::numeric_limits;
using std::runtime_error;
using std::streamsize;
using std::ifstream;
using std::cerr;

namespace gfase {


SamElement::SamElement(string& read_name, string& ref_name, int64_t line, int16_t flag, int8_t mapq) :
        read_name(read_name),
        ref_name(ref_name),
        line(line),
        flag(flag),
        mapq(mapq)
{}


//SamElement::SamElement(SamElement&& other) noexcept:
//        read_name(std::move(other.read_name)),
//        ref_name(std::move(other.ref_name)),
//        line(other.line),
//        flag(other.flag),
//        mapq(other.mapq)
//{}


//SamElement::SamElement(const SamElement& other) noexcept:
//        read_name(other.read_name),
//        ref_name(other.ref_name),
//        line(other.line),
//        flag(other.flag),
//        mapq(other.mapq)
//{}


SamElement::SamElement() :
        read_name(),
        ref_name(),
        line(-1),
        flag(-1),
        mapq(-1)
{}


//SamElement& SamElement::operator=(SamElement&& other) noexcept{
//    if (this != other){
//
//    }
//}


bool SamElement::is_first_mate() const {
    return (int16_t(flag) >> 6) & int16_t(1);
}


bool SamElement::is_second_mate() const {
    return (int16_t(flag) >> 7) & int16_t(1);
}


bool SamElement::is_not_primary() const {
    return (int16_t(flag) >> 8) & int16_t(1);
}


bool SamElement::is_supplementary() const {
    return (int16_t(flag) >> 11) & int16_t(1);
}


// Hash/compare using line number to guarantee unique
bool sam_comparator(const SamElement& a, const SamElement& b) {
    return a.line > b.line;
}


void for_element_in_sam_file(path sam_path, const function<void(SamElement& e)>& f){
    ifstream file(sam_path);

    if (not file.is_open() or not file.good()){
        throw runtime_error("ERROR: could not read input file: " + sam_path.string());
    }

    char c;
    size_t n_delimiters = 0;
    int64_t n_lines = 0;

    string read_name;
    string flag_token;
    string ref_name;
    string mapq_token;

    char header_delimiter = '@';

    // If this is a header line, skip to the next line and increment n_lines until we are no longer on a header
    while (file.peek() == header_delimiter){
        cerr << string(1, char(file.peek())) << '\n';
        file.ignore(numeric_limits<streamsize>::max(), '\n');
        n_lines++;
    }

    while (file.get(c)){
        if (c == '\n'){
            n_delimiters = 0;

            SamElement e(read_name, ref_name, n_lines, int16_t(stoi(flag_token)), int8_t(stoi(mapq_token)));

            f(e);

            read_name.clear();
            flag_token.clear();
            ref_name.clear();
            mapq_token.clear();

            n_lines++;
        }
        else if (c == '\t'){
            n_delimiters++;
        }
        else{
            if (n_delimiters == 0){
                read_name += c;
            }
            else if (n_delimiters == 1){
                flag_token += c;
            }
            else if (n_delimiters == 2){
                ref_name += c;
            }
            else if (n_delimiters == 4){
                mapq_token += c;
            }
        }
    }
}

}


//bool operator<(const gfase::SamElement& a, const gfase::SamElement& b) {
//    return a.line < b.line;
//}
//
//
//bool operator>(const gfase::SamElement& a, const gfase::SamElement& b) {
//    return a.line > b.line;
//}
//
//
//bool operator==(const gfase::SamElement& a, const gfase::SamElement& b) {
//    return a.line == b.line;
//}


ostream& operator<<(ostream& o, const gfase::SamElement& a){
    o << a.read_name << '\n';

    o << '\t' << a.ref_name << '\n';
    o << '\t' << "line: " << a.line << '\n';
    o << '\t' << "flags: " << bitset<sizeof(a.flag)*8>(a.flag) << " = " << a.flag << '\n';
    o << '\t' << "pair: " << a.is_first_mate() << ' ' << a.is_second_mate() << '\n';
    o << '\t' << "secondary: " << a.is_not_primary() << '\n';
    o << '\t' << "supplementary: " << a.is_supplementary() << '\n';

    return o;
}

