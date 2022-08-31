#include "SvgPlot.hpp"
#include <cmath>

using std::min;


SvgPlot::SvgPlot(
        path output_path,
        size_t width,
        size_t height,
        size_t x_min,
        size_t x_max,
        size_t y_min,
        size_t y_max,
        bool axes):
        image_output_path(output_path),
        file(output_path),
        width(width),
        height(height),
        x_min(x_min),
        x_max(x_max),
        y_min(y_min),
        y_max(y_max),
        axes(axes)
{
    if (not file.is_open() and file.good()){
        throw runtime_error("ERROR: could not write to file: " + image_output_path.string());
    }

    if (axes){
        padding = (double((x_max-x_min)+(y_max-y_min))/2)*0.1;
        width *= 1.1;
        height *= 1.2;
    }
    else{
        padding = 0;
    }

    file << "<svg width='" << width + padding << "' height='" << height + padding << "' xmlns='http://www.w3.org/2000/svg' "
         << "viewBox='" << x_min << ' ' << y_min << ' ' << x_max + padding << ' ' << y_max + 2*padding << "'>" << '\n';

    // Transformation group means that the svg elements retain their given coordinates in XML, but are shifted
    file << "<g transform='translate(" << padding << ',' << padding << ")' >" << '\n';

    if (axes) {
        add_line(x_min, y_min, x_min, y_max, padding/200, "black");
        add_line(x_min, y_min, x_max, y_min, padding/200, "black");
        add_line(x_max, y_min, x_max, y_max, padding/200, "black");
        add_line(x_min, y_max, x_max, y_max, padding/200, "black");

        double tick_scale = round(2*(double(x_max - x_min)/10))/2;

        // Transformation group for x ticks
        file << "<g transform='translate(" << 0 << ',' << 0 << ")' >" << '\n';

        for (size_t i=0; (i*tick_scale)<=x_max; i++){
            add_line(
                    double(i*tick_scale),
                    double(y_max),
                    double(i*tick_scale),
                    double(y_max+padding/5),
                    padding/200,
                    "black");
        }

        file << "</g>" << '\n';

        // Transformation group for y ticks
        file << "<g transform='translate(" << -padding/5 << ',' << 0 << ")' >" << '\n';

        for (size_t i=0; (i*tick_scale)<=y_max; i++){
            add_line(
                    double(x_min),
                    double(i*tick_scale),
                    double(x_min+padding/5),
                    double(i*tick_scale),
                    padding/200,
                    "black");
        }

        file << "</g>" << '\n';

    }
}


SvgPlot::~SvgPlot(){
    file << "</g>" << '\n';
    file << "</svg>" << '\n';
    file.close();
}


void SvgPlot::check_file(){
    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: cannot plot with closed svg file");
    }
}
