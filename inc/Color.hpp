#ifndef OVERLAP_ANALYSIS_COLOR_HPP
#define OVERLAP_ANALYSIS_COLOR_HPP


#include <vector>
#include <array>
#include <string>

using std::vector;
using std::array;
using std::string;


template <typename T> inline std::string int_to_hex(T val, size_t width=sizeof(T)*2);


string rgb_to_hex(double r, double g, double b);


class ColorMap {
    /// Methods ///
    virtual array<double,3> get_rgb(double x)=0;
};


class MplRainbow: public ColorMap {
public:
    /// Attributes ///
    static const vector <array <double,3> > rgbs;

    /// Methods ///
    MplRainbow()=default;
    array<double,3> get_rgb(double x);
    string get_svg_color(double x);
};


class MplGnuplot: public ColorMap {
public:
    /// Attributes ///
    static const vector <array <double,3> > rgbs;

    /// Methods ///
    MplGnuplot()=default;
    array<double,3> get_rgb(double x);
    string get_svg_color(double x);
};


class Viridis: public ColorMap {
public:
    /// Attributes ///
    static const vector <array <double,3> > rgbs;

    /// Methods ///
    Viridis()=default;
    array<double,3> get_rgb(double x);
    string get_svg_color(double x);
};


class Seismic: public ColorMap {
public:
    /// Attributes ///
    static const vector <array <double,3> > rgbs;

    /// Methods ///
    Seismic()=default;
    array<double,3> get_rgb(double x);
    string get_svg_color(double x);
};


class Coolwarm: public ColorMap {
public:
    /// Attributes ///
    static const vector <array <double,3> > rgbs;

    /// Methods ///
    Coolwarm()=default;
    array<double,3> get_rgb(double x);
    string get_svg_color(double x);
};




#endif //OVERLAP_ANALYSIS_COLOR_HPP
