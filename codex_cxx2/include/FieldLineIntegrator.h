#ifndef CODEX_CXX2_FIELDLINEINTEGRATOR_H
#define CODEX_CXX2_FIELDLINEINTEGRATOR_H

#include "AparFieldModel.h"
#include "Types.h"

class FieldLineIntegrator
{
public:
  explicit FieldLineIntegrator(const AparFieldModel& model);

  void traceLine(const Point3D& seedInd,
                 std::size_t seedOffset,
                 const TraceOptions& options,
                 int maxStatesPerSeed,
                 int maxTrajPerSeed,
                 int maxPuncPerSeed,
                 std::vector<TrajectoryState>& states,
                 std::vector<Point3D>& trajectories,
                 std::vector<PuncturePoint>& punctures,
                 std::vector<int>& stateCountPerSeed,
                 std::vector<int>& trajCountPerSeed,
                 std::vector<int>& punctureCountPerSeed,
                 std::vector<int>& endRegionPerSeed,
                 std::vector<double>& connectionLengthPerSeed,
                 std::vector<double>& ilinePerSeed) const;

  void traceLine(const Point3D& seedInd, const TraceOptions& options, LineTraceResult& out) const;

private:
  const AparFieldModel& model_;

  void rk4Step(const XZPoint& start,
               int yStart,
               int region,
               int direction,
               XZPoint& end) const;
};

#endif
