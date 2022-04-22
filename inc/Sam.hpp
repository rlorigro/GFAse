#ifndef GFASE_SAM_HPP
#define GFASE_SAM_HPP

#include "Filesystem.hpp"

using ghc::filesystem::path;

#include <functional>
#include <iostream>
#include <bitset>
#include <string>

using std::function;
using std::ostream;
using std::bitset;
using std::string;


namespace gfase {


class SamElement {
public:
    string read_name;
    string ref_name;
    int64_t line;
    int16_t flag;
    int8_t mapq;

    SamElement();
    SamElement(string& read_name, string& ref_name, int64_t line, int16_t flag, int8_t mapq);
//    SamElement(SamElement&& other) noexcept;
//    SamElement(const SamElement& other) noexcept;
//    SamElement& operator=(SamElement&& other) noexcept;
    bool is_first_mate() const;
    bool is_second_mate() const;
    bool is_not_primary() const;
    bool is_supplementary() const;
};


void for_element_in_sam_file(path sam_path, const function<void(SamElement& e)>& f);


}

//bool operator<(const gfase::SamElement& a, const gfase::SamElement& b);
//
//bool operator>(const gfase::SamElement& a, const gfase::SamElement& b);
//
//bool operator==(const gfase::SamElement& a, const gfase::SamElement& b);

ostream& operator<<(ostream& o, const gfase::SamElement& a);


namespace std {
template<>
struct less<gfase::SamElement>
{
    bool operator()(const gfase::SamElement& a, const gfase::SamElement& b) const
    {
        return a.line < b.line;
    }
};
}


#endif //GFASE_SAM_HPP
