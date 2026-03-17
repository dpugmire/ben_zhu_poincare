#include "FieldLineIntegrator.h"

#include <algorithm>
#include <cmath>

namespace
{

struct CrossingEval
{
  double xind = 0.0;
  double yind = 0.0;
  double zvalue = 0.0;
  Point3D xyz;
};

CrossingEval evaluateCrossing(const AparFieldModel& model, const LineTraceResult& line, int tc0, int tc1, int direction, double alpha)
{
  const AparData& d = model.data();

  alpha = std::max(0.0, std::min(1.0, alpha));
  const double beta = 1.0 - alpha;

  const TrajectoryState& s0 = line.states[tc0];
  const TrajectoryState& s1 = line.states[tc1];

  CrossingEval out;
  out.xind = beta * s0.xind + alpha * s1.xind;
  out.yind = beta * s0.yind + alpha * s1.yind;

  out.zvalue = beta * s0.rawZ + alpha * s1.rawZ;
  if (std::fabs(s0.rawZ - s1.rawZ) > 1.0)
  {
    const double z0 = d.wrapZ(s0.rawZ);
    const double z1 = d.wrapZ(s1.rawZ);
    out.zvalue = beta * z0 + alpha * z1;
  }

  if (static_cast<int>(std::round(s0.yind)) == d.nypf2 && direction == 1 && out.xind < static_cast<double>(d.ixsep) + 0.5)
  {
    out.yind = beta * s0.yind + alpha * static_cast<double>(d.nypf2 + 1);
  }
  else if (static_cast<int>(std::round(s0.yind)) == (d.nypf1 + 1) && direction == -1 && out.xind < static_cast<double>(d.ixsep) + 0.5)
  {
    out.yind = beta * static_cast<double>(d.nypf2 + 1) + alpha * s1.yind;
    const double shift = model.interp1(d.xiarray, d.shiftAngle, out.xind);
    out.zvalue = d.wrapZ(out.zvalue - shift);
  }
  else if (tc0 > 0)
  {
    const int yPrev = static_cast<int>(std::round(line.states[tc0 - 1].yind));
    if (yPrev == d.nypf2 || yPrev == (d.nypf1 + 1))
    {
      const double z0 = model.interp1(d.ziarray, d.zarray, s0.zind);
      const double z1 = model.interp1(d.ziarray, d.zarray, s1.zind);
      out.zvalue = beta * z0 + alpha * z1;
    }
  }

  out.zvalue = d.wrapZ(out.zvalue);
  out.xyz = model.reconstructPunctureXYZ(out.xind, out.yind, out.zvalue);
  return out;
}

bool hasSignChange(double a, double b)
{
  return (a <= 0.0 && b >= 0.0) || (a >= 0.0 && b <= 0.0);
}

void tryDetectPunctureOnLastSegment(const AparFieldModel& model, LineTraceResult& line, int direction, int npMax, double& lastFitRoot)
{
  if (static_cast<int>(line.punctures.size()) >= npMax)
  {
    return;
  }
  if (line.states.size() < 2 || line.trajectoryXYZ.size() < 2)
  {
    return;
  }

  const int tc1 = static_cast<int>(line.states.size()) - 1;
  const int tc0 = tc1 - 1;

  const double xPrev = line.trajectoryXYZ[tc0].x;
  const double xCurr = line.trajectoryXYZ[tc1].x;
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

  CrossingEval crossing = evaluateCrossing(model, line, tc0, tc1, direction, alpha);
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
    CrossingEval left = evaluateCrossing(model, line, tc0, tc1, direction, aLeft);
    CrossingEval right = evaluateCrossing(model, line, tc0, tc1, direction, aRight);
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

        CrossingEval next = evaluateCrossing(model, line, tc0, tc1, direction, aNext);
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
  puncture.thetaPsi.x = model.thetaFromY(crossing.yind);
  puncture.thetaPsi.y = model.psiFromX(crossing.xind);

  int step = static_cast<int>(std::floor(fitRoot));
  if (step < 1)
  {
    step = 1;
  }
  if (step >= static_cast<int>(line.states.size()))
  {
    step = static_cast<int>(line.states.size()) - 1;
  }
  puncture.step = step;

  line.punctures.push_back(puncture);
  lastFitRoot = fitRoot;
}

} // namespace

FieldLineIntegrator::FieldLineIntegrator(const AparFieldModel& model)
  : model_(model)
{
}

void FieldLineIntegrator::rk4Step(double xStart, int yStart, double zStart, int region, int direction, double& xEnd, double& zEnd) const
{
  constexpr double h = 1.0;
  const double hh = 0.5 * h;
  const double h6 = h / 6.0;

  double dxdy1 = 0.0;
  double dzdy1 = 0.0;
  double dxdy2 = 0.0;
  double dzdy2 = 0.0;
  double dxdy3 = 0.0;
  double dzdy3 = 0.0;
  double dxdy4 = 0.0;
  double dzdy4 = 0.0;

  model_.evaluateStage(xStart, zStart, yStart, region, direction, 0, dxdy1, dzdy1);
  const double x1 = xStart + direction * hh * dxdy1;
  const double z1 = zStart + direction * hh * dzdy1;

  model_.evaluateStage(x1, z1, yStart, region, direction, 1, dxdy2, dzdy2);
  const double x2 = xStart + direction * hh * dxdy2;
  const double z2 = zStart + direction * hh * dzdy2;

  model_.evaluateStage(x2, z2, yStart, region, direction, 1, dxdy3, dzdy3);
  const double x3 = xStart + direction * dxdy3;
  const double z3 = zStart + direction * dzdy3;

  model_.evaluateStage(x3, z3, yStart, region, direction, 2, dxdy4, dzdy4);

  xEnd = xStart + direction * h6 * (dxdy1 + 2.0 * dxdy2 + 2.0 * dxdy3 + dxdy4);
  zEnd = zStart + direction * h6 * (dzdy1 + 2.0 * dzdy2 + 2.0 * dzdy3 + dzdy4);
}

LineTraceResult FieldLineIntegrator::traceLine(double iline, const TraceOptions& options) const
{
  const AparData& d = model_.data();

  LineTraceResult out;
  out.iline = iline;

  if (iline < 1.0 || iline > static_cast<double>(d.nx))
  {
    out.endRegion = 99;
    return out;
  }

  const int nsteps = options.nturns * d.ny;
  out.states.reserve(static_cast<size_t>(std::max(0, nsteps)));
  out.trajectoryXYZ.reserve(static_cast<size_t>(std::max(0, nsteps)));
  out.punctures.reserve(static_cast<size_t>(std::max(0, options.npMax)));

  double xind = iline;
  double xStart = model_.interp1(d.xiarray, d.xarray, xind);
  int yStart = d.jyomp + 1;
  double zStart = d.zarray.empty() ? 0.0 : d.zarray.front();

  int region = 1;
  if (xind < static_cast<double>(d.ixsep) + 0.5)
  {
    region = 0;
    if (yStart < d.nypf1 + 1 || yStart > d.nypf2)
    {
      region = 2;
    }
  }

  const double zind0 = model_.interp1(d.zarray, d.ziarray, zStart);

  TrajectoryState initial;
  initial.turn = 1;
  initial.xind = xind;
  initial.yind = static_cast<double>(yStart);
  initial.zind = zind0;
  initial.region = region;
  initial.segmentLength = 0.0;
  initial.rawZ = zStart;
  out.states.push_back(initial);
  out.trajectoryXYZ.push_back(model_.reconstructTrajectoryXYZ(initial));

  double lastFitRoot = -1.0e30;

  int iturn = 1;
  while (region < 10 && iturn < options.nturns && static_cast<int>(out.states.size()) < nsteps)
  {
    for (int iy = 0; iy < d.ny - 1; ++iy)
    {
      if (region >= 10 || static_cast<int>(out.states.size()) >= nsteps)
      {
        break;
      }

      double xEnd = xStart;
      double zEnd = zStart;
      int yEnd = yStart;

      if (region == 0 && yStart > d.nypf1 && yStart < d.nypf2 + 1)
      {
        rk4Step(xStart, yStart, zStart, region, options.direction, xEnd, zEnd);
        const double rawZEnd = zEnd;

        yEnd = (options.direction == 1) ? (yStart + 1) : (yStart - 1);

        if (xEnd > d.xMax)
        {
          region = 12;
        }
        else if (xEnd < d.xMin)
        {
          region = 11;
        }
        else
        {
          xind = model_.interp1(d.xarray, d.xiarray, xEnd);
          if (xind > static_cast<double>(d.ixsep) + 0.5)
          {
            region = 1;
          }
        }

        if (options.direction == 1 && yStart == d.nypf2 && region == 0)
        {
          const double shift = model_.interp1(d.xiarray, d.shiftAngle, xind);
          zEnd += shift;
          yEnd = d.nypf1 + 1;
        }
        if (options.direction == -1 && yStart == (d.nypf1 + 1) && region == 0)
        {
          const double shift = model_.interp1(d.xiarray, d.shiftAngle, xind);
          zEnd -= shift;
          yEnd = d.nypf2;
        }

        zEnd = d.wrapZ(zEnd);
        const double zind = model_.interp1(d.zarray, d.ziarray, zEnd);

        TrajectoryState step;
        step.turn = iturn;
        step.xind = xind;
        step.yind = static_cast<double>(yEnd);
        step.zind = zind;
        step.region = region;
        step.segmentLength = 0.0;
        step.rawZ = rawZEnd;
        out.states.push_back(step);
        out.trajectoryXYZ.push_back(model_.reconstructTrajectoryXYZ(step));
        tryDetectPunctureOnLastSegment(model_, out, options.direction, options.npMax, lastFitRoot);

        xStart = xEnd;
        yStart = yEnd;
        zStart = zEnd;
      }
      else if (region == 1 || region == 2)
      {
        rk4Step(xStart, yStart, zStart, region, options.direction, xEnd, zEnd);
        const double rawZEnd = zEnd;

        yEnd = (options.direction == 1) ? (yStart + 1) : (yStart - 1);

        if (options.direction == 1 && yStart == d.nypf1 && region == 2)
        {
          yEnd = d.nypf2 + 1;
        }
        else if (options.direction == -1 && yStart == d.nypf2 + 1 && region == 2)
        {
          yEnd = d.nypf1;
        }

        if (xEnd > d.xMax)
        {
          region = 12;
        }
        else if (xEnd < d.xMin)
        {
          region = 11;
        }
        else
        {
          xind = model_.interp1(d.xarray, d.xiarray, xEnd);
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

        zEnd = d.wrapZ(zEnd);
        const double zind = model_.interp1(d.zarray, d.ziarray, zEnd);

        TrajectoryState step;
        step.turn = iturn;
        step.xind = xind;
        step.yind = static_cast<double>(yEnd);
        step.zind = zind;
        step.region = region;
        step.segmentLength = 0.0;
        step.rawZ = rawZEnd;
        out.states.push_back(step);
        out.trajectoryXYZ.push_back(model_.reconstructTrajectoryXYZ(step));
        tryDetectPunctureOnLastSegment(model_, out, options.direction, options.npMax, lastFitRoot);

        xStart = xEnd;
        yStart = yEnd;
        zStart = zEnd;
      }
      else
      {
        break;
      }
    }

    ++iturn;
  }

  out.endRegion = region;

  double length = 0.0;
  for (size_t i = 1; i < out.trajectoryXYZ.size(); ++i)
  {
    const double dx = out.trajectoryXYZ[i].x - out.trajectoryXYZ[i - 1].x;
    const double dy = out.trajectoryXYZ[i].y - out.trajectoryXYZ[i - 1].y;
    const double dz = out.trajectoryXYZ[i].z - out.trajectoryXYZ[i - 1].z;
    const double ds = std::sqrt(dx * dx + dy * dy + dz * dz);
    length += ds;
  }
  out.connectionLength = length;

  return out;
}
