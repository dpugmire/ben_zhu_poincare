#include "MatlabComparator.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

std::vector<std::vector<double>> MatlabComparator::loadDataRows(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Failed to open file: " + path);
    }

    std::vector<std::vector<double>> rows;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        std::istringstream iss(line);
        std::vector<double> values;
        double value = 0.0;
        while (iss >> value) {
            values.push_back(value);
        }

        // Skip header/non-numeric rows.
        if (values.empty()) {
            continue;
        }

        rows.push_back(values);
    }

    return rows;
}

CompareSummary MatlabComparator::compareFiles(const std::string& generatedPath,
                                              const std::string& matlabPath,
                                              const std::string& label,
                                              double tolerance) const {
    const auto a = loadDataRows(generatedPath);
    const auto b = loadDataRows(matlabPath);

    CompareSummary out;
    out.label = label;
    out.aRows = a.size();
    out.bRows = b.size();

    const size_t rows = std::min(a.size(), b.size());
    out.comparedRows = rows;

    double maxAbs = 0.0;
    double sumSq = 0.0;
    size_t count = 0;

    for (size_t i = 0; i < rows; ++i) {
        const auto& ra = a[i];
        const auto& rb = b[i];

        if (ra.size() < 5 || rb.size() < 5) {
            continue;
        }

        // Compare x,y,z columns only (MATLAB-compatible numeric payload).
        for (size_t c = 2; c <= 4; ++c) {
            const double diff = ra[c] - rb[c];
            const double adiff = std::fabs(diff);
            maxAbs = std::max(maxAbs, adiff);
            sumSq += diff * diff;
            ++count;
        }
    }

    out.maxAbsError = maxAbs;
    out.l2Error = (count > 0) ? std::sqrt(sumSq / static_cast<double>(count)) : 0.0;

    out.pass = (out.aRows == out.bRows) && (out.maxAbsError < tolerance);
    return out;
}
