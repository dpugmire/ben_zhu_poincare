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

    void writeLineOutputsFlat(double iline,
                              std::size_t seedIndex,
                              int maxTrajPerSeed,
                              int maxPuncPerSeed,
                              int trajCount,
                              int punctureCount,
                              const std::vector<Point3D>& trajectories,
                              const std::vector<PuncturePoint>& punctures,
                              const std::vector<std::uint8_t>* punctureValid,
                              const std::string& outputDir,
                              const std::string& divertorTag) const;

    void writeCombinedOutputs(const std::vector<LineTraceResult>& lines,
                              const std::string& outputDir) const;

    void writeCombinedOutputsFlat(const std::vector<double>& ilinePerSeed,
                                  const std::vector<int>& trajCountPerSeed,
                                  const std::vector<int>& punctureCountPerSeed,
                                  int maxTrajPerSeed,
                                  int maxPuncPerSeed,
                                  const std::vector<Point3D>& trajectories,
                                  const std::vector<PuncturePoint>& punctures,
                                  const std::vector<std::uint8_t>* punctureValid,
                                  const std::string& outputDir) const;
};

#endif
