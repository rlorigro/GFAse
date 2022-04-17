#include "SAMElement.hpp"


namespace gfase {

SAMElement::SAMElement(string& read_name, string& ref_name, int32_t line, int16_t flag, int8_t mapq) :
        read_name(read_name),
        ref_name(ref_name),
        line(line),
        flag(flag),
        mapq(mapq) {}


//SAMElement::SAMElement(SAMElement&& other) noexcept:
//        read_name(std::move(other.read_name)),
//        ref_name(std::move(other.ref_name)),
//        line(other.line),
//        flag(other.flag),
//        mapq(other.mapq)
//{}


//SAMElement::SAMElement(const SAMElement& other) noexcept:
//        read_name(other.read_name),
//        ref_name(other.ref_name),
//        line(other.line),
//        flag(other.flag),
//        mapq(other.mapq)
//{}


SAMElement::SAMElement() :
        read_name(),
        ref_name(),
        line(-1),
        flag(-1),
        mapq(-1)
{}


//SAMElement& SAMElement::operator=(SAMElement&& other) noexcept{
//    if (this != other){
//
//    }
//}


bool SAMElement::is_first_mate() const {
    return (int16_t(flag) >> 6) & int16_t(1);
}


bool SAMElement::is_second_mate() const {
    return (int16_t(flag) >> 7) & int16_t(1);
}


bool SAMElement::is_not_primary() const {
    return (int16_t(flag) >> 8) & int16_t(1);
}


bool SAMElement::is_supplementary() const {
    return (int16_t(flag) >> 11) & int16_t(1);
}


// Hash/compare using line number to guarantee unique
bool sam_comparator(const SAMElement& a, const SAMElement& b) {
    return a.line > b.line;
}

}


//bool operator<(const gfase::SAMElement& a, const gfase::SAMElement& b) {
//    return a.line < b.line;
//}
//
//
//bool operator>(const gfase::SAMElement& a, const gfase::SAMElement& b) {
//    return a.line > b.line;
//}
//
//
//bool operator==(const gfase::SAMElement& a, const gfase::SAMElement& b) {
//    return a.line == b.line;
//}


ostream& operator<<(ostream& o, const gfase::SAMElement& a){
    o << a.read_name << '\n';

    o << '\t' << a.ref_name << '\n';
    o << '\t' << "line: " << a.line << '\n';
    o << '\t' << "flags: " << bitset<sizeof(a.flag)*8>(a.flag) << " = " << a.flag << '\n';
    o << '\t' << "pair: " << a.is_first_mate() << ' ' << a.is_second_mate() << '\n';
    o << '\t' << "secondary: " << a.is_not_primary() << '\n';
    o << '\t' << "supplementary: " << a.is_supplementary() << '\n';

    return o;
}

