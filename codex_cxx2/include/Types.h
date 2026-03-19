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
    int maxStatesPerSeed = 0;
    int maxTrajPerSeed = 0;
    int maxPuncPerSeed = 0;

    std::vector<Point3D> seeds;
    std::vector<double> ilinePerSeed;
    std::vector<int> endRegionPerSeed;
    std::vector<double> connectionLengthPerSeed;
    std::vector<int> stateCountPerSeed;
    std::vector<int> trajCountPerSeed;
    std::vector<int> punctureCountPerSeed;

    std::vector<TrajectoryState> states;
    std::vector<Point3D> trajectories;
    std::vector<PuncturePoint> punctures;

    size_t seedCount() const {
        return ilinePerSeed.size();
    }

    size_t stateOffset(size_t seedIndex) const {
        return seedIndex * static_cast<size_t>(maxStatesPerSeed);
    }

    size_t trajOffset(size_t seedIndex) const {
        return seedIndex * static_cast<size_t>(maxTrajPerSeed);
    }

    size_t punctureOffset(size_t seedIndex) const {
        return seedIndex * static_cast<size_t>(maxPuncPerSeed);
    }

    void resize(size_t seedsCount, int maxStates, int maxTraj, int maxPunc) {
        maxStatesPerSeed = (maxStates > 0) ? maxStates : 1;
        maxTrajPerSeed = (maxTraj > 0) ? maxTraj : 1;
        maxPuncPerSeed = (maxPunc > 0) ? maxPunc : 1;

        seeds.resize(seedsCount);
        ilinePerSeed.resize(seedsCount, 0.0);
        endRegionPerSeed.resize(seedsCount, 0);
        connectionLengthPerSeed.resize(seedsCount, 0.0);
        stateCountPerSeed.resize(seedsCount, 0);
        trajCountPerSeed.resize(seedsCount, 0);
        punctureCountPerSeed.resize(seedsCount, 0);

        const size_t totalStateSlots = seedsCount * static_cast<size_t>(maxStatesPerSeed);
        const size_t totalTrajSlots = seedsCount * static_cast<size_t>(maxTrajPerSeed);
        const size_t totalPunctureSlots = seedsCount * static_cast<size_t>(maxPuncPerSeed);

        states.resize(totalStateSlots);
        trajectories.resize(totalTrajSlots);
        punctures.resize(totalPunctureSlots);
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
