#ifndef CODEX_CXX2_FIELDLINEINTEGRATOR_H
#define CODEX_CXX2_FIELDLINEINTEGRATOR_H

#include "AparFieldModel.h"
#include "Types.h"

class FieldLineIntegrator
{
public:
  FieldLineIntegrator(const AparFieldModel& model, const TraceOptions& options);

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

  TraceStatus traceLine(const Point3D& seedInd,
                        std::size_t seedIndex,
                        const TraceOutputViews& outputs,
                        int& stateCount,
                        int& trajCount,
                        int& punctureCount,
                        int& endRegion,
                        double& connectionLength,
                        double& iline) const;

  TraceStatus traceLine(const Point3D& seedInd, LineTraceResult& out) const;

private:
  const AparFieldModel& model_;
  TraceOptions options_;
  int maxStatesPerSeed_ = 1;
  int maxTrajPerSeed_ = 1;
  int maxPuncPerSeed_ = 1;

  void rk4Step(const XZPoint& start,
               int yStart,
               int region,
               int direction,
               XZPoint& end) const;
};

#endif
