#ifndef GFASE_SAMELEMENT_HPP
#define GFASE_SAMELEMENT_HPP

#include <string>

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
    SAMElement(SAMElement&& other) noexcept;
    SAMElement(const SAMElement& other);
    bool is_first_mate() const;
    bool is_second_mate() const;
};

bool sam_comparator(const SAMElement& a, const SAMElement& b);

}

bool operator<(const gfase::SAMElement& a, const gfase::SAMElement& b);

bool operator>(const gfase::SAMElement& a, const gfase::SAMElement& b);


#endif //GFASE_SAMELEMENT_HPP
