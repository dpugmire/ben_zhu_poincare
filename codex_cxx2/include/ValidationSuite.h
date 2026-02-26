#ifndef CODEX_CXX2_VALIDATIONSUITE_H
#define CODEX_CXX2_VALIDATIONSUITE_H

#include <string>
#include <vector>

#include "Types.h"

struct ValidationConfig {
    std::string referenceDir;
    std::string outputDir;
    std::string aparSinglePath;
    std::string aparCircPath;
    std::string divertorFilter = "all";  // all|single|circ
    std::vector<int> linesFilter;          // empty => all discovered
    TraceOptions traceOptions;
    double tolerance = 1.0e-8;
};

class ValidationSuite {
public:
    std::vector<ValidationResult> run(const ValidationConfig& config) const;

private:
    static std::vector<ValidationCase> discoverCases(const ValidationConfig& config);
};

#endif
