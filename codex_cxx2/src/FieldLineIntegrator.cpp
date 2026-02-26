#include "FieldLineIntegrator.h"

#include <cmath>

FieldLineIntegrator::FieldLineIntegrator(const AparFieldModel& model) : model_(model) {}

void FieldLineIntegrator::rk4Step(double xStart,
                                  int yStart,
                                  double zStart,
                                  int region,
                                  int direction,
                                  double& xEnd,
                                  double& zEnd) const {
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

LineTraceResult FieldLineIntegrator::traceLine(int iline, const TraceOptions& options) const {
    const AparData& d = model_.data();

    LineTraceResult out;
    out.iline = iline;

    if (iline < 1 || iline > d.nx) {
        out.endRegion = 99;
        return out;
    }

    const int nsteps = options.nturns * d.ny;

    double xind = static_cast<double>(iline);
    double xStart = model_.interp1(d.xiarray, d.xarray, xind);
    int yStart = d.jyomp + 1;
    double zStart = d.zarray.empty() ? 0.0 : d.zarray.front();

    int region = 1;
    if (xind < static_cast<double>(d.ixsep) + 0.5) {
        region = 0;
        if (yStart < d.nypf1 + 1 || yStart > d.nypf2) {
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

    int iturn = 1;
    while (region < 10 && iturn < options.nturns && static_cast<int>(out.states.size()) < nsteps) {
        for (int iy = 0; iy < d.ny - 1; ++iy) {
            if (region >= 10 || static_cast<int>(out.states.size()) >= nsteps) {
                break;
            }

            double xEnd = xStart;
            double zEnd = zStart;
            int yEnd = yStart;

            if (region == 0 && yStart > d.nypf1 && yStart < d.nypf2 + 1) {
                rk4Step(xStart, yStart, zStart, region, options.direction, xEnd, zEnd);
                const double rawZEnd = zEnd;

                yEnd = (options.direction == 1) ? (yStart + 1) : (yStart - 1);

                if (xEnd > d.xMax) {
                    region = 12;
                } else if (xEnd < d.xMin) {
                    region = 11;
                } else {
                    xind = model_.interp1(d.xarray, d.xiarray, xEnd);
                    if (xind > static_cast<double>(d.ixsep) + 0.5) {
                        region = 1;
                    }
                }

                if (options.direction == 1 && yStart == d.nypf2 && region == 0) {
                    const double shift = model_.interp1(d.xiarray, d.shiftAngle, xind);
                    zEnd += shift;
                    yEnd = d.nypf1 + 1;
                }
                if (options.direction == -1 && yStart == (d.nypf1 + 1) && region == 0) {
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

                xStart = xEnd;
                yStart = yEnd;
                zStart = zEnd;
            } else if (region == 1 || region == 2) {
                rk4Step(xStart, yStart, zStart, region, options.direction, xEnd, zEnd);
                const double rawZEnd = zEnd;

                yEnd = (options.direction == 1) ? (yStart + 1) : (yStart - 1);

                if (options.direction == 1 && yStart == d.nypf1 && region == 2) {
                    yEnd = d.nypf2 + 1;
                } else if (options.direction == -1 && yStart == d.nypf2 + 1 && region == 2) {
                    yEnd = d.nypf1;
                }

                if (xEnd > d.xMax) {
                    region = 12;
                } else if (xEnd < d.xMin) {
                    region = 11;
                } else {
                    xind = model_.interp1(d.xarray, d.xiarray, xEnd);
                    if (xind < static_cast<double>(d.ixsep) + 0.5 && yEnd > d.nypf1 && yEnd < d.nypf2 + 1) {
                        region = 0;
                    } else if (xind < static_cast<double>(d.ixsep) + 0.5 &&
                               (yEnd > d.nypf2 - 1 || yEnd < d.nypf1)) {
                        region = 2;
                    }
                }

                if (options.direction == 1 && yEnd == d.ny) {
                    region = 14;
                } else if (options.direction == -1 && yEnd == 1) {
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

                xStart = xEnd;
                yStart = yEnd;
                zStart = zEnd;
            } else {
                break;
            }
        }

        ++iturn;
    }

    out.endRegion = region;

    out.trajectoryXYZ.reserve(out.states.size());
    for (const TrajectoryState& state : out.states) {
        out.trajectoryXYZ.push_back(model_.reconstructTrajectoryXYZ(state));
    }

    double length = 0.0;
    for (size_t i = 1; i < out.trajectoryXYZ.size(); ++i) {
        const double dx = out.trajectoryXYZ[i].x - out.trajectoryXYZ[i - 1].x;
        const double dy = out.trajectoryXYZ[i].y - out.trajectoryXYZ[i - 1].y;
        const double dz = out.trajectoryXYZ[i].z - out.trajectoryXYZ[i - 1].z;
        const double ds = std::sqrt(dx * dx + dy * dy + dz * dz);
        length += ds;
    }
    out.connectionLength = length;

    return out;
}
