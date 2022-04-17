#ifndef GFASE_SAMELEMENT_HPP
#define GFASE_SAMELEMENT_HPP

#include <iostream>
#include <bitset>
#include <string>

using std::ostream;
using std::bitset;
using std::string;

namespace gfase {


class SAMElement {
public:
    string read_name;
    string ref_name;
    int32_t line;
    int16_t flag;
    int8_t mapq;

    SAMElement();
    SAMElement(string& read_name, string& ref_name, int32_t line, int16_t flag, int8_t mapq);
//    SAMElement(SAMElement&& other) noexcept;
//    SAMElement(const SAMElement& other) noexcept;
//    SAMElement& operator=(SAMElement&& other) noexcept;
    bool is_first_mate() const;
    bool is_second_mate() const;
    bool is_not_primary() const;
    bool is_supplementary() const;
};

}

//bool operator<(const gfase::SAMElement& a, const gfase::SAMElement& b);
//
//bool operator>(const gfase::SAMElement& a, const gfase::SAMElement& b);
//
//bool operator==(const gfase::SAMElement& a, const gfase::SAMElement& b);

ostream& operator<<(ostream& o, const gfase::SAMElement& a);


namespace std {
template<>
struct less<gfase::SAMElement>
{
    bool operator()(const gfase::SAMElement& a, const gfase::SAMElement& b) const
    {
        return a.line < b.line;
    }
};
}


#endif //GFASE_SAMELEMENT_HPP
