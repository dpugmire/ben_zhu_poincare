#ifndef CODEX_CXX2_MATLABCOMPARATOR_H
#define CODEX_CXX2_MATLABCOMPARATOR_H

#include <string>
#include <vector>

#include "Types.h"

class MatlabComparator {
public:
    CompareSummary compareFiles(const std::string& generatedPath,
                                const std::string& matlabPath,
                                const std::string& label,
                                double tolerance) const;

private:
    static std::vector<std::vector<double>> loadDataRows(const std::string& path);
};

#endif
