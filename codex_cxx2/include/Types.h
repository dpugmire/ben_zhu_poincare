#ifndef CODEX_CXX2_TYPES_H
#define CODEX_CXX2_TYPES_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "ViskoresCompat.h"

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

enum class TraceStatus {
    Ok = 0,
    InvalidSeed = 1,
    InvalidConfiguration = 2,
    OutputTooSmall = 3,
    MaxStepLimitReached = 4
};

inline const char* traceStatusName(TraceStatus status) {
    switch (status) {
        case TraceStatus::Ok:
            return "ok";
        case TraceStatus::InvalidSeed:
            return "invalid-seed";
        case TraceStatus::InvalidConfiguration:
            return "invalid-configuration";
        case TraceStatus::OutputTooSmall:
            return "output-too-small";
        case TraceStatus::MaxStepLimitReached:
            return "max-step-limit-reached";
    }
    return "unknown";
}

struct TraceOutputViews {
    TrajectoryState* states = nullptr;
    std::size_t statesSize = 0;

    Point3D* trajectories = nullptr;
    std::size_t trajectoriesSize = 0;

    PuncturePoint* punctures = nullptr;
    std::size_t puncturesSize = 0;

    std::uint8_t* punctureValid = nullptr;  // optional
    std::size_t punctureValidSize = 0;
};

inline TraceOutputViews makeTraceOutputViews(std::vector<TrajectoryState>& states,
                                             std::vector<Point3D>& trajectories,
                                             std::vector<PuncturePoint>& punctures,
                                             std::vector<std::uint8_t>* punctureValid = nullptr) {
    TraceOutputViews out;
    out.states = states.empty() ? nullptr : states.data();
    out.statesSize = states.size();
    out.trajectories = trajectories.empty() ? nullptr : trajectories.data();
    out.trajectoriesSize = trajectories.size();
    out.punctures = punctures.empty() ? nullptr : punctures.data();
    out.puncturesSize = punctures.size();
    if (punctureValid != nullptr) {
        out.punctureValid = punctureValid->empty() ? nullptr : punctureValid->data();
        out.punctureValidSize = punctureValid->size();
    }
    return out;
}

#endif
