#ifndef CODEX_CXX2_TYPES_H
#define CODEX_CXX2_TYPES_H

#include <cstddef>
#include <string>
#include <vector>

struct Point2D {
    double x = 0.0;
    double y = 0.0;
};

struct Point3D {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct XZPoint {
    double x = 0.0;
    double z = 0.0;
};

struct XZDeriv {
    double dxdy = 0.0;
    double dzdy = 0.0;
};

struct TrajectoryState {
    int turn = 0;
    Point3D ind;         // MATLAB 1-based continuous index-space point (xind,yind,zind)
    int region = 0;
    double segmentLength = 0.0;
    double rawZ = 0.0;   // unwrapped zEnd for branch-cut-safe puncture interpolation
};

struct PuncturePoint {
    int step = 0;
    Point3D xyz;
    Point2D thetaPsi;
};

struct LineTraceResult {
    double iline = 0.0;
    int endRegion = 0;
    double connectionLength = 0.0;
    std::vector<TrajectoryState> states;
    std::vector<Point3D> trajectoryXYZ;
    std::vector<PuncturePoint> punctures;
};

struct TraceOptions {
    int direction = 1;
    int npMax = 100;
    int maxSteps = 0;   // <=0 => auto cap: kDefaultMaxStepsPerPuncture * npMax
};

struct PackedLineTraceBatch {
    int maxStatesPerLine = 0;
    int maxPuncturesPerLine = 0;

    std::vector<Point3D> seeds;
    std::vector<double> iline;
    std::vector<int> endRegion;
    std::vector<double> connectionLength;
    std::vector<int> stateCount;
    std::vector<int> trajectoryCount;
    std::vector<int> punctureCount;

    std::vector<TrajectoryState> states;
    std::vector<Point3D> trajectoryXYZ;
    std::vector<PuncturePoint> punctures;

    std::vector<bool> stateSet;
    std::vector<bool> trajectorySet;
    std::vector<bool> punctureSet;

    size_t lineCount() const {
        return iline.size();
    }

    size_t stateOffset(size_t lineIndex) const {
        return lineIndex * static_cast<size_t>(maxStatesPerLine);
    }

    size_t punctureOffset(size_t lineIndex) const {
        return lineIndex * static_cast<size_t>(maxPuncturesPerLine);
    }

    void resize(size_t lines, int maxStates, int maxPunctures) {
        maxStatesPerLine = (maxStates > 0) ? maxStates : 1;
        maxPuncturesPerLine = (maxPunctures > 0) ? maxPunctures : 1;

        seeds.assign(lines, Point3D{});
        iline.assign(lines, 0.0);
        endRegion.assign(lines, 0);
        connectionLength.assign(lines, 0.0);
        stateCount.assign(lines, 0);
        trajectoryCount.assign(lines, 0);
        punctureCount.assign(lines, 0);

        const size_t totalStateSlots = lines * static_cast<size_t>(maxStatesPerLine);
        const size_t totalPunctureSlots = lines * static_cast<size_t>(maxPuncturesPerLine);

        states.assign(totalStateSlots, TrajectoryState{});
        trajectoryXYZ.assign(totalStateSlots, Point3D{});
        punctures.assign(totalPunctureSlots, PuncturePoint{});

        stateSet.assign(totalStateSlots, false);
        trajectorySet.assign(totalStateSlots, false);
        punctureSet.assign(totalPunctureSlots, false);
    }
};

struct ValidationCase {
    std::string divertorTag;
    std::string aparFile;
    int line = 0;
    std::string refIpPath;
    std::string refTrajPath;
};

struct CompareSummary {
    std::string label;
    size_t comparedRows = 0;
    size_t aRows = 0;
    size_t bRows = 0;
    double maxAbsError = 0.0;
    double l2Error = 0.0;
    bool pass = false;
};

struct ValidationResult {
    ValidationCase testCase;
    CompareSummary ipSummary;
    CompareSummary trajSummary;
    bool pass = false;
};

#endif
