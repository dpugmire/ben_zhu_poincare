#ifndef CODEX_CXX2_POINCAREOUTPUT_H
#define CODEX_CXX2_POINCAREOUTPUT_H

#include <cstddef>
#include <string>

#include "Types.h"

class PoincareOutput {
public:
    void writeLineOutputs(const LineTraceResult& line,
                          const std::string& outputDir,
                          const std::string& divertorTag) const;

    void writeLineOutputs(const PackedLineTraceBatch& batch,
                          std::size_t lineIndex,
                          const std::string& outputDir,
                          const std::string& divertorTag) const;

    void writeCombinedOutputs(const std::vector<LineTraceResult>& lines,
                              const std::string& outputDir) const;

    void writeCombinedOutputs(const PackedLineTraceBatch& batch,
                              const std::string& outputDir) const;
};

#endif
