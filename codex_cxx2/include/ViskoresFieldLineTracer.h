#ifndef CODEX_CXX2_VISKORESFIELDLINETRACER_H
#define CODEX_CXX2_VISKORESFIELDLINETRACER_H

#include <vector>

#include "AparData.h"
#include "Types.h"

#if defined(CODEX_USE_VISKORES)

class ViskoresFieldLineTracer
{
public:
  ViskoresFieldLineTracer(const AparData& data, const TraceOptions& options);

  int maxStatesPerSeed() const
  {
    return maxStatesPerSeed_;
  }

  int maxTrajPerSeed() const
  {
    return maxTrajPerSeed_;
  }

  int maxPuncPerSeed() const
  {
    return maxPuncPerSeed_;
  }

  void traceLines(const std::vector<Point3D>& seeds,
                  const std::vector<CodeXId>& globalSeedIndices,
                  const TraceOutputViews& outputs,
                  std::vector<double>& ilinePerSeed,
                  std::vector<int>& endRegionPerSeed,
                  std::vector<double>& connectionLengthPerSeed,
                  std::vector<int>& stateCountPerSeed,
                  std::vector<int>& trajCountPerSeed,
                  std::vector<int>& punctureCountPerSeed,
                  std::vector<TraceStatus>& traceStatuses,
                  double* deviceInvokeSeconds = nullptr,
                  double* hostPostprocessSeconds = nullptr) const;

private:
  const AparData& data_;
  TraceOptions options_;
  int maxStatesPerSeed_ = 1;
  int maxTrajPerSeed_ = 1;
  int maxPuncPerSeed_ = 1;
};

#endif

#endif
