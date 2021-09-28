#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "graph_utility.hpp"
#include "Filesystem.hpp"
#include "GfaReader.hpp"
#include "CLI11.hpp"

#include "bdsg/hash_graph.hpp"

#include <string>
#include <array>

using gfase::IncrementalIdMap;
using gfase::for_each_connected_component;
using gfase::split_connected_components;
using gfase::handle_graph_to_gfa;
using gfase::print_graph_paths;
using gfase::plot_graph;
using ghc::filesystem::path;

using bdsg::HashGraph;
using bdsg::MutablePathMutableHandleGraph;
using bdsg::MutablePathDeletableHandleGraph;
using handlegraph::path_handle_t;
using handlegraph::step_handle_t;
using handlegraph::handle_t;

using std::string;
using std::cout;
using std::cerr;

///
/// Aqua 	#00FFFF 	0,255,255
///  Aquamarine 	#7FFFD4 	127,255,212
///  Beige 	#F5F5DC 	245,245,220
///  Black 	#000000 	0,0,0
///  Blue 	#0000FF 	0,0,255
///  BlueViolet 	#8A2BE2 	138,43,226
///  Brown 	#A52A2A 	165,42,42
///  CadetBlue 	#5F9EA0 	95,158,160
///  Chartreuse 	#7FFF00 	127,255,0
///  Coral 	#FF7F50 	255,127,80
///  CornflowerBlue 	#6495ED 	100,149,237
///  Crimson 	#DC143C 	220,20,60
///  Cyan 	#00FFFF 	0,255,255
///  DarkBlue 	#00008B 	0,0,139
///  DarkCyan 	#008B8B 	0,139,139
///  DarkGoldenrod 	#B8860B 	184,134,11
///  DarkGray 	#A9A9A9 	169,169,169
///  DarkGreen 	#006400 	0,100,0
///  DarkKhaki 	#BDB76B 	189,183,107
///  DarkMagenta 	#8B008B 	139,0,139
///  DarkOliveGreen 	#556B2F 	85,107,47
///  DarkOrange 	#FF8C00 	255,140,0
///  DarkOrchid 	#9932CC 	153,50,204
///  DarkRed 	#8B0000 	139,0,0
///  DarkSalmon 	#E9967A 	233,150,122
///  DarkSeaGreen 	#8FBC8F 	143,188,143
///  DarkSlateBlue 	#483D8B 	72,61,139
///  DarkSlateGray 	#2F4F4F 	47,79,79
///  DarkTurquoise 	#00CED1 	0,206,209
///  DarkViolet 	#9400D3 	148,0,211
///  DeepPink 	#FF1493 	255,20,147
///  DeepSkyBlue 	#00BFFF 	0,191,255
///  DodgerBlue 	#1E90FF 	30,144,255
///  FireBrick 	#B22222 	178,34,34
///  ForestGreen 	#228B22 	34,139,34
///  Fuchsia 	#FF00FF 	255,0,255
///  Gold 	#FFD700 	255,215,0
///  Goldenrod 	#DAA520 	218,165,32
///  Gray 	#808080 	128,128,128
///  Green 	#008000 	0,128,0
///  GreenYellow 	#ADFF2F 	173,255,47
///  HotPink 	#FF69B4 	255,105,180
///  IndianRed 	#CD5C5C 	205,92,92
///  Indigo 	#4B0082 	75,0,130
///  Khaki 	#F0E68C 	240,230,140
///  Lavender 	#E6E6FA 	230,230,250
///  LavenderBlush 	#FFF0F5 	255,240,245
///  LawnGreen 	#7CFC00 	124,252,0
///  LightBlue 	#ADD8E6 	173,216,230
///  LightCoral 	#F08080 	240,128,128
///  LightGreen 	#90EE90 	144,238,144
///  LightPink 	#FFB6C1 	255,182,193
///  LightSalmon 	#FFA07A 	255,160,122
///  LightSeaGreen 	#20B2AA 	32,178,170
///  LightSkyBlue 	#87CEFA 	135,206,250
///  Lime 	#00FF00 	0,255,0
///  LimeGreen 	#32CD32 	50,205,50
///  Magenta 	#FF00FF 	255,0,255
///  Maroon 	#800000 	128,0,0
///  MediumAquamarine 	#66CDAA 	102,205,170
///  MediumBlue 	#0000CD 	0,0,205
///  MediumOrchid 	#BA55D3 	186,85,211
///  MediumPurple 	#9370DB 	147,112,219
///  MediumSeaGreen 	#3CB371 	60,179,113
///  MediumSlateBlue 	#7B68EE 	123,104,238
///  MediumSpringGreen 	#00FA9A 	0,250,154
///  MediumTurquoise 	#48D1CC 	72,209,204
///  MediumVioletRed 	#C71585 	199,21,133
///  MidnightBlue 	#191970 	25,25,112
///  Navy 	#000080 	0,0,128
///  Olive 	#808000 	128,128,0
///  OliveDrab 	#6B8E23 	107,142,35
///  Orange 	#FFA500 	255,165,0
///  OrangeRed 	#FF4500 	255,69,0
///  Orchid 	#DA70D6 	218,112,214
///  PaleGreen 	#98FB98 	152,251,152
///  PaleTurquoise 	#AFEEEE 	175,238,238
///  PaleVioletRed 	#DB7093 	219,112,147
///  PeachPuff 	#FFDAB9 	255,218,185
///  Pink 	#FFC0CB 	255,192,203
///  Plum 	#DDA0DD 	221,160,221
///  PowderBlue 	#B0E0E6 	176,224,230
///  Purple 	#800080 	128,0,128
///  Red 	#FF0000 	255,0,0
///  RosyBrown 	#BC8F8F 	188,143,143
///  RoyalBlue 	#4169E1 	65,105,225
///  Salmon 	#FA8072 	250,128,114
///  SandyBrown 	#F4A460 	244,164,96
///  SeaGreen 	#2E8B57 	46,139,87
///  SkyBlue 	#87CEEB 	135,206,235
///  SlateBlue 	#6A5ACD 	106,90,205
///  SpringGreen 	#00FF7F 	0,255,127
///  SteelBlue 	#4682B4 	70,130,180
///  Teal 	#008080 	0,128,128
///  Thistle 	#D8BFD8 	216,191,216
///  Tomato 	#FF6347 	255,99,71
///  Turquoise 	#40E0D0 	64,224,208
///  Violet 	#EE82EE 	238,130,238
///  Yellow 	#FFFF00 	255,255,0
///  YellowGreen 	#9ACD32 	154,205,50
///


class Colors{
public:
    static const vector <string> hex;
    static const vector <char> hex_index_to_char;
    static const unordered_map<char,uint8_t> char_to_hex_index;
    static uint8_t blend_value(uint8_t a, uint8_t b, float t);
    static array<uint8_t,3> blend_rgb(array<uint8_t,3> a, array<uint8_t,3> b);
    static string rgb_to_hex(array<uint8_t,3> rgb);
    static array<uint8_t,3> hex_to_rgb(string hex);
    static string dec_to_hex(uint8_t dec);
};


const vector<string> Colors::hex = {
        "#00FFFF",      // Aqua
        "#7FFFD4",      // Aquamarine
        "#F5F5DC",      // Beige
        "#000000",      // Black
        "#0000FF",      // Blue
        "#8A2BE2",      // BlueViolet
        "#A52A2A",      // Brown
        "#5F9EA0",      // CadetBlue
        "#7FFF00",      // Chartreuse
        "#FF7F50",      // Coral
        "#6495ED",      // CornflowerBlue
        "#DC143C",      // Crimson
        "#00FFFF",      // Cyan
        "#00008B",      // DarkBlue
        "#008B8B",      // DarkCyan
        "#B8860B",      // DarkGoldenrod
        "#A9A9A9",      // DarkGray
        "#006400",      // DarkGreen
        "#BDB76B",      // DarkKhaki
        "#8B008B",      // DarkMagenta
        "#556B2F",      // DarkOliveGreen
        "#FF8C00",      // DarkOrange
        "#9932CC",      // DarkOrchid
        "#8B0000",      // DarkRed
        "#E9967A",      // DarkSalmon
        "#8FBC8F",      // DarkSeaGreen
        "#483D8B",      // DarkSlateBlue
        "#2F4F4F",      // DarkSlateGray
        "#00CED1",      // DarkTurquoise
        "#9400D3",      // DarkViolet
        "#FF1493",      // DeepPink
        "#00BFFF",      // DeepSkyBlue
        "#1E90FF",      // DodgerBlue
        "#B22222",      // FireBrick
        "#228B22",      // ForestGreen
        "#FF00FF",      // Fuchsia
        "#FFD700",      // Gold
        "#DAA520",      // Goldenrod
        "#808080",      // Gray
        "#008000",      // Green
        "#ADFF2F",      // GreenYellow
        "#FF69B4",      // HotPink
        "#CD5C5C",      // IndianRed
        "#4B0082",      // Indigo
        "#F0E68C",      // Khaki
        "#E6E6FA",      // Lavender
        "#FFF0F5",      // LavenderBlush
        "#7CFC00",      // LawnGreen
        "#ADD8E6",      // LightBlue
        "#F08080",      // LightCoral
        "#90EE90",      // LightGreen
        "#FFB6C1",      // LightPink
        "#FFA07A",      // LightSalmon
        "#20B2AA",      // LightSeaGreen
        "#87CEFA",      // LightSkyBlue
        "#00FF00",      // Lime
        "#32CD32",      // LimeGreen
        "#FF00FF",      // Magenta
        "#800000",      // Maroon
        "#66CDAA",      // MediumAquamarine
        "#0000CD",      // MediumBlue
        "#BA55D3",      // MediumOrchid
        "#9370DB",      // MediumPurple
        "#3CB371",      // MediumSeaGreen
        "#7B68EE",      // MediumSlateBlue
        "#00FA9A",      // MediumSpringGreen
        "#48D1CC",      // MediumTurquoise
        "#C71585",      // MediumVioletRed
        "#191970",      // MidnightBlue
        "#000080",      // Navy
        "#808000",      // Olive
        "#6B8E23",      // OliveDrab
        "#FFA500",      // Orange
        "#FF4500",      // OrangeRed
        "#DA70D6",      // Orchid
        "#98FB98",      // PaleGreen
        "#AFEEEE",      // PaleTurquoise
        "#DB7093",      // PaleVioletRed
        "#FFDAB9",      // PeachPuff
        "#FFC0CB",      // Pink
        "#DDA0DD",      // Plum
        "#B0E0E6",      // PowderBlue
        "#800080",      // Purple
        "#FF0000",      // Red
        "#BC8F8F",      // RosyBrown
        "#4169E1",      // RoyalBlue
        "#FA8072",      // Salmon
        "#F4A460",      // SandyBrown
        "#2E8B57",      // SeaGreen
        "#87CEEB",      // SkyBlue
        "#6A5ACD",      // SlateBlue
        "#00FF7F",      // SpringGreen
        "#4682B4",      // SteelBlue
        "#008080",      // Teal
        "#D8BFD8",      // Thistle
        "#FF6347",      // Tomato
        "#40E0D0",      // Turquoise
        "#EE82EE",      // Violet
        "#FFFF00",      // Yellow
        "#9ACD32"       // YellowGreen
};


const vector<char> Colors::hex_index_to_char = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};
const unordered_map<char,uint8_t> Colors::char_to_hex_index = {
        {'0',0},
        {'1',1},
        {'2',2},
        {'3',3},
        {'4',4},
        {'5',5},
        {'6',6},
        {'7',7},
        {'8',8},
        {'9',9},
        {'A',10},
        {'B',11},
        {'C',12},
        {'D',13},
        {'E',14},
        {'F',15}
};


uint8_t Colors::blend_value(uint8_t a, uint8_t b, float t){
    double a_normal = double(a)/255;
    double b_normal = double(b)/255;
    double c_normal = sqrt((1-t) * pow(a_normal,2.2) + t * pow(b_normal,2.2));

    return uint8_t(round(255 * c_normal));
}


array<uint8_t,3> Colors::blend_rgb(array<uint8_t,3> a, array<uint8_t,3> b){
    array<uint8_t,3> blended_rgb;

    for (size_t i=0; i<blended_rgb.size(); i++){
        blended_rgb[i] = blend_value(a[i], b[i], 0.5);
    }

    return blended_rgb;
}


string Colors::dec_to_hex(uint8_t dec){
    uint16_t v1 = dec & 0b1111;
    uint16_t v2 = (dec >> 4) & 0b1111;

    return string(1,Colors::hex_index_to_char.at(v2)) + string(1,Colors::hex_index_to_char.at(v1));
}


string Colors::rgb_to_hex(array<uint8_t,3> rgb){
    return "#" + dec_to_hex(rgb[0]) + dec_to_hex(rgb[1]) + dec_to_hex(rgb[2]);
}


array<uint8_t,3> Colors::hex_to_rgb(string hex){
    hex = hex.substr(1,hex.size()-1);
    array<uint8_t,3> rgb = {0,0,0};

    for (size_t i=0; i<3; i++){
        auto h1 = hex[2*i];
        auto h2 = hex[2*i+1];

        uint8_t v = 0;
        v |= Colors::char_to_hex_index.at(h1);
        v <<= 4;
        v &= 0b11110000;
        v |= Colors::char_to_hex_index.at(h2);

        rgb[i] = v;
    }

    return rgb;
}


void create_color_table(path gfa_path){
    HashGraph graph;
    IncrementalIdMap<string> id_map;

    gfa_to_handle_graph(graph, id_map, gfa_path, false);

    path out_path = gfa_path;
    out_path.replace_extension(".path_colors.csv");

    ofstream file(out_path);
    file << "Name" << ',' << "Color" << ',' << "Paths" << '\n';

    graph.for_each_handle([&](const handle_t& h){
        auto node_name = id_map.get_name(graph.get_id(h));

        unordered_set<string> paths_on_handle;

        graph.for_each_step_on_handle(h, [&](const step_handle_t& s){
            auto p = graph.get_path_handle_of_step(s);
            auto path_name = graph.get_path_name(p);

            paths_on_handle.emplace(path_name);
        });

        string hex;

        if (paths_on_handle.size() == 1){
            string path_name = *paths_on_handle.begin();

            auto hash_code = std::hash<std::string>()(path_name);
            hex = Colors::hex.at(hash_code % Colors::hex.size());
            auto rgb = Colors::hex_to_rgb(hex);
            auto hex2 = Colors::rgb_to_hex(rgb);

            cerr << hex << " " << to_string(rgb[0]) << ',' << to_string(rgb[1]) << ',' << to_string(rgb[2]) << " " << hex2 << '\n';
        }
        else if (paths_on_handle.size() == 2){
            auto iter = paths_on_handle.begin();

            string path_name_a = *iter;
            iter++;
            string path_name_b = *iter;

            auto hash_code_a = std::hash<std::string>()(path_name_a);
            auto hash_code_b = std::hash<std::string>()(path_name_b);

            auto hex_a = Colors::hex.at(hash_code_a % Colors::hex.size());
            auto hex_b = Colors::hex.at(hash_code_b % Colors::hex.size());

            auto rgb_a = Colors::hex_to_rgb(hex_a);
            auto rgb_b = Colors::hex_to_rgb(hex_b);

            auto blended_rgb = Colors::blend_rgb(rgb_a, rgb_b);

            hex = Colors::rgb_to_hex(blended_rgb);
        }
        else if (paths_on_handle.size() > 2){
            hex = "#FFFFFF";
        }
        else{
            // Near-black for no path
            hex = "#242424";
        }

        string paths;
        for (auto& path_name: paths_on_handle){
            paths += path_name + "_";
        }
        paths = paths.substr(0, paths.size()-1);

        file << node_name << ',' << hex << ',' << paths << '\n';
    });
}


int main (int argc, char* argv[]){
    path gfa_path;

    CLI::App app{"App description"};

    app.add_option(
            "-i,--input_gfa",
            gfa_path,
            "Path to GFA containing phased non-overlapping segments")
            ->required();

    CLI11_PARSE(app, argc, argv);

    create_color_table(gfa_path);

    return 0;
}