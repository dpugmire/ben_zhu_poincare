#include "FieldLineIntegrator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{

struct CrossingEval
{
  Point2D ind;
  double zvalue = 0.0;
  Point3D xyz;
};

int computeMaxStateCountFromOptions(const TraceOptions& options)
{
  constexpr int kDefaultMaxStepsPerPuncture = 200;
  const long long derivedMaxStepsLL = static_cast<long long>(kDefaultMaxStepsPerPuncture) *
                                      static_cast<long long>(std::max(1, options.npMax));
  const int derivedMaxSteps = (derivedMaxStepsLL > static_cast<long long>(std::numeric_limits<int>::max()))
      ? std::numeric_limits<int>::max()
      : static_cast<int>(std::max(1LL, derivedMaxStepsLL));

  const int maxStepsCap = (options.maxSteps > 0) ? options.maxSteps : derivedMaxSteps;
  const int nsteps = std::max(1, maxStepsCap);
  return (nsteps >= std::numeric_limits<int>::max())
      ? std::numeric_limits<int>::max()
      : (nsteps + 1);
}

CrossingEval evaluateCrossing(const AparFieldModel& model,
                              const std::vector<TrajectoryState>& states,
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
  else if (static_cast<int>(std::round(s0.ind.y)) == (d.nypf1 + 1) && direction == -1 && out.ind.x < static_cast<double>(d.ixsep) + 0.5)
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

void tryDetectPunctureOnLastSegment(const AparFieldModel& model,
                                    const std::vector<TrajectoryState>& states,
                                    const std::vector<Point3D>& trajectories,
                                    std::size_t stateBase,
                                    std::size_t trajBase,
                                    int stateCount,
                                    int direction,
                                    int maxPuncturesForSeed,
                                    std::vector<PuncturePoint>& punctures,
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

  punctures[punctureBase + static_cast<std::size_t>(punctureCount)] = puncture;
  ++punctureCount;
  lastFitRoot = fitRoot;
}

} // namespace

FieldLineIntegrator::FieldLineIntegrator(const AparFieldModel& model)
  : model_(model)
{
}

void FieldLineIntegrator::rk4Step(const XZPoint& start,
                                  int yStart,
                                  int region,
                                  int direction,
                                  XZPoint& end) const
{
  constexpr double h = 1.0;
  const double hh = 0.5 * h;
  const double h6 = h / 6.0;

  XZDeriv k1;
  XZDeriv k2;
  XZDeriv k3;
  XZDeriv k4;

  model_.evaluateStage(start, yStart, region, direction, 0, k1);
  const XZPoint p1{start.x + direction * hh * k1.dxdy,
                   start.z + direction * hh * k1.dzdy};

  model_.evaluateStage(p1, yStart, region, direction, 1, k2);
  const XZPoint p2{start.x + direction * hh * k2.dxdy,
                   start.z + direction * hh * k2.dzdy};

  model_.evaluateStage(p2, yStart, region, direction, 1, k3);
  const XZPoint p3{start.x + direction * k3.dxdy,
                   start.z + direction * k3.dzdy};

  model_.evaluateStage(p3, yStart, region, direction, 2, k4);

  end.x = start.x + direction * h6 * (k1.dxdy + 2.0 * k2.dxdy + 2.0 * k3.dxdy + k4.dxdy);
  end.z = start.z + direction * h6 * (k1.dzdy + 2.0 * k2.dzdy + 2.0 * k3.dzdy + k4.dzdy);
}

void FieldLineIntegrator::traceLine(const Point3D& seedInd,
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
                                    std::vector<double>& ilinePerSeed) const
{
  if (seedOffset >= ilinePerSeed.size() ||
      seedOffset >= endRegionPerSeed.size() ||
      seedOffset >= connectionLengthPerSeed.size() ||
      seedOffset >= stateCountPerSeed.size() ||
      seedOffset >= trajCountPerSeed.size() ||
      seedOffset >= punctureCountPerSeed.size())
  {
    throw std::runtime_error("seedOffset out of range for per-seed arrays");
  }

  const std::size_t stateBase = seedOffset * static_cast<std::size_t>(std::max(0, maxStatesPerSeed));
  const std::size_t trajBase = seedOffset * static_cast<std::size_t>(std::max(0, maxTrajPerSeed));
  const std::size_t punctureBase = seedOffset * static_cast<std::size_t>(std::max(0, maxPuncPerSeed));

  if (maxStatesPerSeed <= 0 || maxTrajPerSeed <= 0 || maxPuncPerSeed <= 0)
  {
    ilinePerSeed[seedOffset] = seedInd.x;
    endRegionPerSeed[seedOffset] = 98;
    connectionLengthPerSeed[seedOffset] = 0.0;
    stateCountPerSeed[seedOffset] = 0;
    trajCountPerSeed[seedOffset] = 0;
    punctureCountPerSeed[seedOffset] = 0;
    return;
  }

  if (stateBase + static_cast<std::size_t>(maxStatesPerSeed) > states.size() ||
      trajBase + static_cast<std::size_t>(maxTrajPerSeed) > trajectories.size() ||
      punctureBase + static_cast<std::size_t>(maxPuncPerSeed) > punctures.size())
  {
    throw std::runtime_error("traceLine output arrays are smaller than seeds.size() * maxValue");
  }

  const AparData& d = model_.data();

  ilinePerSeed[seedOffset] = seedInd.x;
  endRegionPerSeed[seedOffset] = 0;
  connectionLengthPerSeed[seedOffset] = 0.0;
  stateCountPerSeed[seedOffset] = 0;
  trajCountPerSeed[seedOffset] = 0;
  punctureCountPerSeed[seedOffset] = 0;

  const double xindSeed = seedInd.x;
  if (xindSeed < 1.0 || xindSeed > static_cast<double>(d.nx))
  {
    endRegionPerSeed[seedOffset] = 99;
    return;
  }

  const int maxStateCountFromOptions = computeMaxStateCountFromOptions(options);
  const int maxStateCount = std::min(std::min(maxStateCountFromOptions, maxStatesPerSeed), maxTrajPerSeed);
  const int maxPunctureCount = std::min(std::max(0, options.npMax), maxPuncPerSeed);
  if (maxStateCount <= 0)
  {
    endRegionPerSeed[seedOffset] = 98;
    return;
  }

  double xind = xindSeed;
  XZPoint current;
  current.x = model_.interp1(d.xiarray, d.xarray, xind);
  int yStart = static_cast<int>(std::llround(seedInd.y));
  if (yStart < 1)
  {
    yStart = 1;
  }
  else if (yStart > d.ny)
  {
    yStart = d.ny;
  }
  const double zind0 = seedInd.z;
  current.z = model_.interp1(d.ziarray, d.zarray, zind0);

  int region = 1;
  if (xind < static_cast<double>(d.ixsep) + 0.5)
  {
    region = 0;
    if (yStart < d.nypf1 + 1 || yStart > d.nypf2)
    {
      region = 2;
    }
  }

  int stateCount = 0;
  int trajCount = 0;
  int punctureCount = 0;

  TrajectoryState initial;
  initial.turn = 1;
  initial.ind.x = xind;
  initial.ind.y = static_cast<double>(yStart);
  initial.ind.z = zind0;
  initial.region = region;
  initial.segmentLength = 0.0;
  initial.rawZ = current.z;

  states[stateBase + static_cast<std::size_t>(stateCount)] = initial;
  trajectories[trajBase + static_cast<std::size_t>(trajCount)] = model_.reconstructTrajectoryXYZ(initial);
  ++stateCount;
  ++trajCount;

  double lastFitRoot = -1.0e30;

  int iturn = 1;
  while (region < 10 && stateCount < maxStateCount)
  {
    for (int iy = 0; iy < d.ny - 1; ++iy)
    {
      if (region >= 10 || stateCount >= maxStateCount)
      {
        break;
      }

      XZPoint next = current;
      int yEnd = yStart;

      if (region == 0 && yStart > d.nypf1 && yStart < d.nypf2 + 1)
      {
        rk4Step(current, yStart, region, options.direction, next);
        const double rawZEnd = next.z;

        yEnd = (options.direction == 1) ? (yStart + 1) : (yStart - 1);

        if (next.x > d.xMax)
        {
          region = 12;
        }
        else if (next.x < d.xMin)
        {
          region = 11;
        }
        else
        {
          xind = model_.interp1(d.xarray, d.xiarray, next.x);
          if (xind > static_cast<double>(d.ixsep) + 0.5)
          {
            region = 1;
          }
        }

        if (options.direction == 1 && yStart == d.nypf2 && region == 0)
        {
          const double shift = model_.interp1(d.xiarray, d.shiftAngle, xind);
          next.z += shift;
          yEnd = d.nypf1 + 1;
        }
        if (options.direction == -1 && yStart == (d.nypf1 + 1) && region == 0)
        {
          const double shift = model_.interp1(d.xiarray, d.shiftAngle, xind);
          next.z -= shift;
          yEnd = d.nypf2;
        }

        next.z = d.wrapZ(next.z);
        const double zind = model_.interp1(d.zarray, d.ziarray, next.z);

        TrajectoryState step;
        step.turn = iturn;
        step.ind.x = xind;
        step.ind.y = static_cast<double>(yEnd);
        step.ind.z = zind;
        step.region = region;
        step.segmentLength = 0.0;
        step.rawZ = rawZEnd;

        states[stateBase + static_cast<std::size_t>(stateCount)] = step;
        trajectories[trajBase + static_cast<std::size_t>(trajCount)] = model_.reconstructTrajectoryXYZ(step);
        ++stateCount;
        ++trajCount;

        if (maxPunctureCount > 0)
        {
          tryDetectPunctureOnLastSegment(model_,
                                         states,
                                         trajectories,
                                         stateBase,
                                         trajBase,
                                         stateCount,
                                         options.direction,
                                         maxPunctureCount,
                                         punctures,
                                         punctureBase,
                                         punctureCount,
                                         lastFitRoot);
        }

        current = next;
        yStart = yEnd;
      }
      else if (region == 1 || region == 2)
      {
        rk4Step(current, yStart, region, options.direction, next);
        const double rawZEnd = next.z;

        yEnd = (options.direction == 1) ? (yStart + 1) : (yStart - 1);

        if (options.direction == 1 && yStart == d.nypf1 && region == 2)
        {
          yEnd = d.nypf2 + 1;
        }
        else if (options.direction == -1 && yStart == d.nypf2 + 1 && region == 2)
        {
          yEnd = d.nypf1;
        }

        if (next.x > d.xMax)
        {
          region = 12;
        }
        else if (next.x < d.xMin)
        {
          region = 11;
        }
        else
        {
          xind = model_.interp1(d.xarray, d.xiarray, next.x);
          if (xind < static_cast<double>(d.ixsep) + 0.5 && yEnd > d.nypf1 && yEnd < d.nypf2 + 1)
          {
            region = 0;
          }
          else if (xind < static_cast<double>(d.ixsep) + 0.5 && (yEnd > d.nypf2 - 1 || yEnd < d.nypf1))
          {
            region = 2;
          }
        }

        if (options.direction == 1 && yEnd == d.ny)
        {
          region = 14;
        }
        else if (options.direction == -1 && yEnd == 1)
        {
          region = 13;
        }

        next.z = d.wrapZ(next.z);
        const double zind = model_.interp1(d.zarray, d.ziarray, next.z);

        TrajectoryState step;
        step.turn = iturn;
        step.ind.x = xind;
        step.ind.y = static_cast<double>(yEnd);
        step.ind.z = zind;
        step.region = region;
        step.segmentLength = 0.0;
        step.rawZ = rawZEnd;

        states[stateBase + static_cast<std::size_t>(stateCount)] = step;
        trajectories[trajBase + static_cast<std::size_t>(trajCount)] = model_.reconstructTrajectoryXYZ(step);
        ++stateCount;
        ++trajCount;

        if (maxPunctureCount > 0)
        {
          tryDetectPunctureOnLastSegment(model_,
                                         states,
                                         trajectories,
                                         stateBase,
                                         trajBase,
                                         stateCount,
                                         options.direction,
                                         maxPunctureCount,
                                         punctures,
                                         punctureBase,
                                         punctureCount,
                                         lastFitRoot);
        }

        current = next;
        yStart = yEnd;
      }
      else
      {
        break;
      }
    }

    ++iturn;
  }

  double length = 0.0;
  for (int i = 1; i < trajCount; ++i)
  {
    const Point3D& p1 = trajectories[trajBase + static_cast<std::size_t>(i - 1)];
    const Point3D& p2 = trajectories[trajBase + static_cast<std::size_t>(i)];
    const double dx = p2.x - p1.x;
    const double dy = p2.y - p1.y;
    const double dz = p2.z - p1.z;
    length += std::sqrt(dx * dx + dy * dy + dz * dz);
  }

  endRegionPerSeed[seedOffset] = region;
  connectionLengthPerSeed[seedOffset] = length;
  stateCountPerSeed[seedOffset] = stateCount;
  trajCountPerSeed[seedOffset] = trajCount;
  punctureCountPerSeed[seedOffset] = punctureCount;
}

void FieldLineIntegrator::traceLine(const Point3D& seedInd, const TraceOptions& options, LineTraceResult& out) const
{
  const int maxStatesPerSeed = computeMaxStateCountFromOptions(options);
  const int maxTrajPerSeed = maxStatesPerSeed;
  const int maxPuncPerSeed = std::max(1, options.npMax);

  std::vector<TrajectoryState> states(static_cast<std::size_t>(maxStatesPerSeed));
  std::vector<Point3D> trajectories(static_cast<std::size_t>(maxTrajPerSeed));
  std::vector<PuncturePoint> punctures(static_cast<std::size_t>(maxPuncPerSeed));

  std::vector<int> stateCountPerSeed(1, 0);
  std::vector<int> trajCountPerSeed(1, 0);
  std::vector<int> punctureCountPerSeed(1, 0);
  std::vector<int> endRegionPerSeed(1, 0);
  std::vector<double> connectionLengthPerSeed(1, 0.0);
  std::vector<double> ilinePerSeed(1, 0.0);

  traceLine(seedInd,
            0,
            options,
            maxStatesPerSeed,
            maxTrajPerSeed,
            maxPuncPerSeed,
            states,
            trajectories,
            punctures,
            stateCountPerSeed,
            trajCountPerSeed,
            punctureCountPerSeed,
            endRegionPerSeed,
            connectionLengthPerSeed,
            ilinePerSeed);

  out.iline = ilinePerSeed[0];
  out.endRegion = endRegionPerSeed[0];
  out.connectionLength = connectionLengthPerSeed[0];

  const int stateCount = std::max(0, std::min(stateCountPerSeed[0], maxStatesPerSeed));
  const int trajCount = std::max(0, std::min(trajCountPerSeed[0], maxTrajPerSeed));
  const int punctureCount = std::max(0, std::min(punctureCountPerSeed[0], maxPuncPerSeed));

  out.states.assign(states.begin(), states.begin() + stateCount);
  out.trajectoryXYZ.assign(trajectories.begin(), trajectories.begin() + trajCount);
  out.punctures.assign(punctures.begin(), punctures.begin() + punctureCount);
}
