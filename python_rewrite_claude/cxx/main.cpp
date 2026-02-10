#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "field_data.h"
#include "trace.h"

namespace {

void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " --apar-file FILE --grid-file FILE [options]\n"
              << "Options:\n"
              << "  --direction DIR   Tracing direction: 1 or -1 (default: 1)\n"
              << "  --nlines N        Number of lines to trace (default: 256)\n"
              << "  --nturns N        Number of turns (default: 25)\n"
              << "  --lines LIST      Comma-separated list, e.g. \"50,75,100\"\n"
              << "  --help            Show this message\n";
}

std::vector<int> parse_lines_csv(const std::string& text) {
    std::vector<int> lines;
    size_t pos = 0;
    while (pos < text.size()) {
        size_t comma = text.find(',', pos);
        if (comma == std::string::npos) {
            comma = text.size();
        }
        std::string token = text.substr(pos, comma - pos);
        if (!token.empty()) {
            lines.push_back(std::atoi(token.c_str()));
        }
        pos = comma + 1;
    }
    return lines;
}

}  // namespace

int main(int argc, char** argv) {
    std::string apar_file;
    std::string grid_file;
    int direction = 1;
    int nlines = 256;
    int nturns = 25;
    std::vector<int> lines;
    bool lines_given = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "--apar-file" && i + 1 < argc) {
            apar_file = argv[++i];
        } else if (arg == "--grid-file" && i + 1 < argc) {
            grid_file = argv[++i];
        } else if (arg == "--direction" && i + 1 < argc) {
            direction = std::atoi(argv[++i]);
        } else if (arg == "--nlines" && i + 1 < argc) {
            nlines = std::atoi(argv[++i]);
        } else if (arg == "--nturns" && i + 1 < argc) {
            nturns = std::atoi(argv[++i]);
        } else if (arg == "--lines" && i + 1 < argc) {
            lines_given = true;
            lines = parse_lines_csv(argv[++i]);
        } else {
            std::cerr << "Unknown or incomplete argument: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    if (apar_file.empty() || grid_file.empty()) {
        std::cerr << "Error: --apar-file and --grid-file are required.\n";
        print_usage(argv[0]);
        return 1;
    }
    if (direction != 1 && direction != -1) {
        std::cerr << "Error: --direction must be 1 or -1.\n";
        return 1;
    }

    TraceOptions options;
    options.direction = direction;
    options.nlines = nlines;
    options.nturns = nturns;
    if (lines_given) {
        options.lines = lines;
    }

    try {
        FieldData field;
        field.initialize(grid_file, apar_file);
        trace_field_lines(field, options);
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
