#ifndef CODEX_CXX2_TYPES_H
#define CODEX_CXX2_TYPES_H

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

struct TrajectoryState {
    int turn = 0;
    double xind = 0.0;   // MATLAB 1-based continuous x-index
    double yind = 0.0;   // MATLAB 1-based y-index (integer-valued in stepping)
    double zind = 0.0;   // MATLAB 1-based continuous z-index
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
    int nturns = 100;
    int npMax = 1250;
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
