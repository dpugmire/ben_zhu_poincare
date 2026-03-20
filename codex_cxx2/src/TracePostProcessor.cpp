#include "TracePostProcessor.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace
{

struct CrossingEval
{
  Point2D ind;
  double zvalue = 0.0;
  Point3D xyz;
};

CrossingEval evaluateCrossing(const AparFieldModel::ExecutionAccessor& model,
                              const TrajectoryState* states,
                              std::size_t stateBase,
                              int tc0,
                              int tc1,
                              int direction,
                              double alpha)
{
  const AparData& d = model.data();

  alpha = std::max(0.0, std::min(1.0, alpha));
  const double beta = 1.0 - alpha;

  const TrajectoryState& s0 = states[stateBase + static_cast<std::size_t>(tc0)];
  const TrajectoryState& s1 = states[stateBase + static_cast<std::size_t>(tc1)];

  CrossingEval out;
  out.ind.x = beta * s0.ind.x + alpha * s1.ind.x;
  out.ind.y = beta * s0.ind.y + alpha * s1.ind.y;

  out.zvalue = beta * s0.rawZ + alpha * s1.rawZ;
  if (std::fabs(s0.rawZ - s1.rawZ) > 1.0)
  {
    const double z0 = d.wrapZ(s0.rawZ);
    const double z1 = d.wrapZ(s1.rawZ);
    out.zvalue = beta * z0 + alpha * z1;
  }

  if (static_cast<int>(std::round(s0.ind.y)) == d.nypf2 && direction == 1 && out.ind.x < static_cast<double>(d.ixsep) + 0.5)
  {
    out.ind.y = beta * s0.ind.y + alpha * static_cast<double>(d.nypf2 + 1);
  }
  else if (static_cast<int>(std::round(s0.ind.y)) == (d.nypf1 + 1) && direction == -1 &&
           out.ind.x < static_cast<double>(d.ixsep) + 0.5)
  {
    out.ind.y = beta * static_cast<double>(d.nypf2 + 1) + alpha * s1.ind.y;
    const double shift = model.interp1(d.xiarray, d.shiftAngle, out.ind.x);
    out.zvalue = d.wrapZ(out.zvalue - shift);
  }
  else if (tc0 > 0)
  {
    const int yPrev = static_cast<int>(std::round(states[stateBase + static_cast<std::size_t>(tc0 - 1)].ind.y));
    if (yPrev == d.nypf2 || yPrev == (d.nypf1 + 1))
    {
      const double z0 = model.interp1(d.ziarray, d.zarray, s0.ind.z);
      const double z1 = model.interp1(d.ziarray, d.zarray, s1.ind.z);
      out.zvalue = beta * z0 + alpha * z1;
    }
  }

  out.zvalue = d.wrapZ(out.zvalue);
  out.xyz = model.reconstructPunctureXYZ(out.ind, out.zvalue);
  return out;
}

bool hasSignChange(double a, double b)
{
  return (a <= 0.0 && b >= 0.0) || (a >= 0.0 && b <= 0.0);
}

void tryDetectPunctureOnLastSegment(const AparFieldModel::ExecutionAccessor& model,
                                    const TrajectoryState* states,
                                    const Point3D* trajectories,
                                    std::size_t stateBase,
                                    std::size_t trajBase,
                                    int stateCount,
                                    int direction,
                                    int maxPuncturesForSeed,
                                    PuncturePoint* punctures,
                                    std::uint8_t* punctureValid,
                                    std::size_t punctureBase,
                                    int& punctureCount,
                                    double& lastFitRoot)
{
  if (punctureCount >= maxPuncturesForSeed)
  {
    return;
  }
  if (stateCount < 2)
  {
    return;
  }

  const int tc1 = stateCount - 1;
  const int tc0 = tc1 - 1;

  const double xPrev = trajectories[trajBase + static_cast<std::size_t>(tc0)].x;
  const double xCurr = trajectories[trajBase + static_cast<std::size_t>(tc1)].x;
  if (!hasSignChange(xPrev, xCurr))
  {
    return;
  }

  double alpha = 0.5;
  const double denom = xCurr - xPrev;
  if (std::fabs(denom) > 1.0e-20)
  {
    alpha = -xPrev / denom;
  }
  alpha = std::max(0.0, std::min(1.0, alpha));

  CrossingEval crossing = evaluateCrossing(model, states, stateBase, tc0, tc1, direction, alpha);
  double bestAlpha = alpha;
  double bestAbsX = std::fabs(crossing.xyz.x);
  CrossingEval bestCrossing = crossing;

  constexpr double xTol = 1.0e-4;
  constexpr double alphaTol = 1.0e-6;
  constexpr int maxIter = 10;

  if (bestAbsX > xTol)
  {
    double aLeft = 0.0;
    double aRight = 1.0;
    CrossingEval left = evaluateCrossing(model, states, stateBase, tc0, tc1, direction, aLeft);
    CrossingEval right = evaluateCrossing(model, states, stateBase, tc0, tc1, direction, aRight);
    double fLeft = left.xyz.x;
    double fRight = right.xyz.x;

    if (std::fabs(fLeft) < bestAbsX)
    {
      bestAbsX = std::fabs(fLeft);
      bestAlpha = aLeft;
      bestCrossing = left;
    }
    if (std::fabs(fRight) < bestAbsX)
    {
      bestAbsX = std::fabs(fRight);
      bestAlpha = aRight;
      bestCrossing = right;
    }

    if (hasSignChange(fLeft, fRight))
    {
      for (int iter = 0; iter < maxIter; ++iter)
      {
        if ((aRight - aLeft) <= alphaTol)
        {
          break;
        }

        double aNext = 0.5 * (aLeft + aRight);
        const double secantDenom = fRight - fLeft;
        if (std::fabs(secantDenom) > 1.0e-20)
        {
          const double secant = aLeft - fLeft * (aRight - aLeft) / secantDenom;
          if (secant > aLeft + alphaTol && secant < aRight - alphaTol)
          {
            aNext = secant;
          }
        }

        CrossingEval next = evaluateCrossing(model, states, stateBase, tc0, tc1, direction, aNext);
        const double fNext = next.xyz.x;
        const double absNext = std::fabs(fNext);
        if (absNext < bestAbsX)
        {
          bestAbsX = absNext;
          bestAlpha = aNext;
          bestCrossing = next;
        }
        if (absNext <= xTol)
        {
          break;
        }

        if (hasSignChange(fLeft, fNext))
        {
          aRight = aNext;
          fRight = fNext;
        }
        else
        {
          aLeft = aNext;
          fLeft = fNext;
        }
      }
    }
  }

  alpha = bestAlpha;
  crossing = bestCrossing;

  const double fitRoot = static_cast<double>(tc0 + 1) + alpha;
  constexpr double dedupEps = 1.0e-5;
  if (std::fabs(fitRoot - lastFitRoot) < dedupEps)
  {
    return;
  }

  if (crossing.xyz.y <= 0.0)
  {
    return;
  }

  PuncturePoint puncture;
  puncture.xyz = crossing.xyz;
  puncture.thetaPsi.x = model.thetaFromY(crossing.ind.y);
  puncture.thetaPsi.y = model.psiFromX(crossing.ind.x);

  int step = static_cast<int>(std::floor(fitRoot));
  if (step < 1)
  {
    step = 1;
  }
  if (step >= stateCount)
  {
    step = stateCount - 1;
  }
  puncture.step = step;

  const std::size_t punctureIndex = punctureBase + static_cast<std::size_t>(punctureCount);
  punctures[punctureIndex] = puncture;
  if (punctureValid != nullptr)
  {
    punctureValid[punctureIndex] = static_cast<std::uint8_t>(1);
  }
  ++punctureCount;
  lastFitRoot = fitRoot;
}

} // namespace

namespace TracePostProcessor
{

void rebuildSeedOutputs(const AparFieldModel& model,
                        const TraceOptions& options,
                        std::size_t seedIndex,
                        int maxStatesPerSeed,
                        int maxTrajPerSeed,
                        int maxPuncPerSeed,
                        const TraceOutputViews& outputs,
                        int stateCount,
                        int& trajCount,
                        int& punctureCount,
                        double& connectionLength)
{
  trajCount = 0;
  punctureCount = 0;
  connectionLength = 0.0;

  if (outputs.states == nullptr || outputs.trajectories == nullptr || outputs.punctures == nullptr)
  {
    return;
  }
  if (maxStatesPerSeed <= 0 || maxTrajPerSeed <= 0 || maxPuncPerSeed <= 0)
  {
    return;
  }

  const std::size_t stateBase = seedIndex * static_cast<std::size_t>(maxStatesPerSeed);
  const std::size_t trajBase = seedIndex * static_cast<std::size_t>(maxTrajPerSeed);
  const std::size_t punctureBase = seedIndex * static_cast<std::size_t>(maxPuncPerSeed);

  const int clampedStateCount = std::max(0, std::min(stateCount, maxStatesPerSeed));
  const int clampedTrajCount = std::max(0, std::min(clampedStateCount, maxTrajPerSeed));
  const int maxPunctureCount = std::max(0, std::min(options.npMax, maxPuncPerSeed));

  if (stateBase + static_cast<std::size_t>(clampedStateCount) > outputs.statesSize ||
      trajBase + static_cast<std::size_t>(clampedTrajCount) > outputs.trajectoriesSize ||
      punctureBase + static_cast<std::size_t>(maxPuncPerSeed) > outputs.puncturesSize)
  {
    return;
  }

  if (outputs.punctureValid != nullptr)
  {
    if (punctureBase + static_cast<std::size_t>(maxPuncPerSeed) > outputs.punctureValidSize)
    {
      return;
    }
    std::fill(outputs.punctureValid + punctureBase,
              outputs.punctureValid + punctureBase + static_cast<std::size_t>(maxPuncPerSeed),
              static_cast<std::uint8_t>(0));
  }

  auto modelExec = model.prepareExecution();

  double lastFitRoot = -1.0e30;
  Point3D prevTrajectory;
  bool havePrevTrajectory = false;

  for (int i = 0; i < clampedTrajCount; ++i)
  {
    const std::size_t stateIndex = stateBase + static_cast<std::size_t>(i);
    const std::size_t trajIndex = trajBase + static_cast<std::size_t>(i);

    outputs.trajectories[trajIndex] = model.reconstructTrajectoryXYZ(outputs.states[stateIndex]);
    const Point3D& currentTrajectory = outputs.trajectories[trajIndex];

    if (havePrevTrajectory)
    {
      const double dx = currentTrajectory.x - prevTrajectory.x;
      const double dy = currentTrajectory.y - prevTrajectory.y;
      const double dz = currentTrajectory.z - prevTrajectory.z;
      connectionLength += std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    prevTrajectory = currentTrajectory;
    havePrevTrajectory = true;
    trajCount = i + 1;

    if (maxPunctureCount > 0)
    {
      tryDetectPunctureOnLastSegment(modelExec,
                                     outputs.states,
                                     outputs.trajectories,
                                     stateBase,
                                     trajBase,
                                     i + 1,
                                     options.direction,
                                     maxPunctureCount,
                                     outputs.punctures,
                                     outputs.punctureValid,
                                     punctureBase,
                                     punctureCount,
                                     lastFitRoot);
    }
  }
}

} // namespace TracePostProcessor
