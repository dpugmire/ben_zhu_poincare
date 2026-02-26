#include "PunctureDetector.h"

#include <algorithm>
#include <cmath>

#include "Interpolator.h"

namespace {

struct CrossingEval {
    double xind = 0.0;
    double yind = 0.0;
    double zvalue = 0.0;
    Point3D xyz;
};

CrossingEval evaluateCrossing(const AparFieldModel& model,
                              const LineTraceResult& line,
                              int tc0,
                              int tc1,
                              int direction,
                              double alpha) {
    const AparData& d = model.data();

    alpha = std::max(0.0, std::min(1.0, alpha));
    const double beta = 1.0 - alpha;

    const TrajectoryState& s0 = line.states[tc0];
    const TrajectoryState& s1 = line.states[tc1];

    CrossingEval out;
    out.xind = beta * s0.xind + alpha * s1.xind;
    out.yind = beta * s0.yind + alpha * s1.yind;

    out.zvalue = beta * s0.rawZ + alpha * s1.rawZ;
    if (std::fabs(s0.rawZ - s1.rawZ) > 1.0) {
        const double z0 = d.wrapZ(s0.rawZ);
        const double z1 = d.wrapZ(s1.rawZ);
        out.zvalue = beta * z0 + alpha * z1;
    }

    if (static_cast<int>(std::round(s0.yind)) == d.nypf2 &&
        direction == 1 &&
        out.xind < static_cast<double>(d.ixsep) + 0.5) {
        out.yind = beta * s0.yind + alpha * static_cast<double>(d.nypf2 + 1);
    } else if (static_cast<int>(std::round(s0.yind)) == (d.nypf1 + 1) &&
               direction == -1 &&
               out.xind < static_cast<double>(d.ixsep) + 0.5) {
        out.yind = beta * static_cast<double>(d.nypf2 + 1) + alpha * s1.yind;
        const double shift = model.interp1(d.xiarray, d.shiftAngle, out.xind);
        out.zvalue = d.wrapZ(out.zvalue - shift);
    } else if (tc0 > 0) {
        const int yPrev = static_cast<int>(std::round(line.states[tc0 - 1].yind));
        if (yPrev == d.nypf2 || yPrev == (d.nypf1 + 1)) {
            const double z0 = model.interp1(d.ziarray, d.zarray, s0.zind);
            const double z1 = model.interp1(d.ziarray, d.zarray, s1.zind);
            out.zvalue = beta * z0 + alpha * z1;
        }
    }

    out.zvalue = d.wrapZ(out.zvalue);
    out.xyz = model.reconstructPunctureXYZ(out.xind, out.yind, out.zvalue);
    return out;
}

}  // namespace

PunctureDetector::PunctureDetector(const AparFieldModel& model) : model_(model) {}

void PunctureDetector::detect(LineTraceResult& line,
                              int direction,
                              int npMax) const {
    line.punctures.clear();

    if (line.states.size() < 2 || line.trajectoryXYZ.size() < 2) {
        return;
    }

    const size_t n = line.trajectoryXYZ.size();

    std::vector<double> fitX(n, 0.0);
    std::vector<double> fitIt(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        fitX[i] = line.trajectoryXYZ[i].x;
        fitIt[i] = static_cast<double>(i + 1);
    }

    NaturalCubicSpline spline;
    spline.build(fitIt, fitX);

    double lastFitRoot = -1.0e30;
    constexpr double endpointEps = 1.0e-6;
    constexpr double dedupEps = 1.0e-5;

    for (size_t seg = 0; seg + 1 < n && static_cast<int>(line.punctures.size()) < npMax; ++seg) {
        const std::vector<double> roots = Interpolator::solveCubicSegmentRoots(
            spline.values()[seg], spline.b()[seg], spline.c()[seg], spline.d()[seg]);

        std::vector<double> localRoots = roots;
        if (localRoots.empty()) {
            const double x0 = fitX[seg];
            const double x1 = fitX[seg + 1];
            if (x0 * x1 <= 0.0) {
                double alpha = 0.5;
                const double denom = x1 - x0;
                if (std::fabs(denom) > 1.0e-20) {
                    alpha = -x0 / denom;
                }
                alpha = std::max(0.0, std::min(1.0, alpha));
                localRoots.push_back(alpha);
            }
        }

        for (double alpha : localRoots) {
            if (alpha <= endpointEps && seg > 0) {
                continue;
            }
            if (alpha >= 1.0 - endpointEps && seg + 2 < n) {
                continue;
            }

            const double fitRoot = static_cast<double>(seg + 1) + alpha;
            if (std::fabs(fitRoot - lastFitRoot) < dedupEps) {
                continue;
            }
            lastFitRoot = fitRoot;

            const int tc0 = static_cast<int>(seg);
            const int tc1 = static_cast<int>(seg + 1);
            CrossingEval crossing = evaluateCrossing(model_, line, tc0, tc1, direction, alpha);

            if (crossing.xyz.y > 0.0) {
                PuncturePoint puncture;
                puncture.xyz = crossing.xyz;
                puncture.thetaPsi.x = model_.thetaFromY(crossing.yind);
                puncture.thetaPsi.y = model_.psiFromX(crossing.xind);

                int step = static_cast<int>(std::floor(fitRoot));
                if (step < 1) {
                    step = 1;
                }
                if (step >= static_cast<int>(line.states.size())) {
                    step = static_cast<int>(line.states.size()) - 1;
                }
                puncture.step = step;

                line.punctures.push_back(puncture);

                if (static_cast<int>(line.punctures.size()) >= npMax) {
                    break;
                }
            }
        }
    }
}
