#include "PoincareOutput.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace {

std::string formatLineToken(double line) {
    if (std::isfinite(line) && std::fabs(line - std::round(line)) < 1.0e-12) {
        return std::to_string(static_cast<long long>(std::llround(line)));
    }
    std::ostringstream oss;
    oss << std::setprecision(12) << std::defaultfloat << line;
    return oss.str();
}

void writeLineToStreams(const LineTraceResult& line,
                        std::ofstream& ipOut,
                        std::ofstream& tpOut,
                        std::ofstream& trajOut) {
    const std::string lineToken = formatLineToken(line.iline);
    for (size_t i = 0; i < line.trajectoryXYZ.size(); ++i) {
        const Point3D& p = line.trajectoryXYZ[i];
        trajOut << lineToken << " " << (i + 1) << " "
                << p.x << " " << p.y << " " << p.z << "\n";
    }

    for (const PuncturePoint& p : line.punctures) {
        ipOut << lineToken << " " << p.step << " "
              << p.xyz.x << " " << p.xyz.y << " " << p.xyz.z << "\n";
        tpOut << lineToken << " " << p.step << " "
              << 0.0 << " " << p.thetaPsi.x << " " << p.thetaPsi.y << "\n";
    }
}

void writeFlatLineToStreams(double iline,
                            std::size_t seedIndex,
                            int maxTrajPerSeed,
                            int maxPuncPerSeed,
                            int trajCount,
                            int punctureCount,
                            const std::vector<Point3D>& trajectories,
                            const std::vector<PuncturePoint>& punctures,
                            std::ofstream& ipOut,
                            std::ofstream& tpOut,
                            std::ofstream& trajOut) {
    if (maxTrajPerSeed <= 0 || maxPuncPerSeed <= 0) {
        return;
    }

    const std::size_t trajBase = seedIndex * static_cast<std::size_t>(maxTrajPerSeed);
    const std::size_t punctureBase = seedIndex * static_cast<std::size_t>(maxPuncPerSeed);
    const std::string lineToken = formatLineToken(iline);

    const int clampedTrajCount = std::max(0, std::min(trajCount, maxTrajPerSeed));
    const int clampedPunctureCount = std::max(0, std::min(punctureCount, maxPuncPerSeed));

    if (trajBase + static_cast<std::size_t>(clampedTrajCount) > trajectories.size() ||
        punctureBase + static_cast<std::size_t>(clampedPunctureCount) > punctures.size()) {
        throw std::runtime_error("Output arrays are smaller than numSeeds * maxXPerSeed");
    }

    for (int i = 0; i < clampedTrajCount; ++i) {
        const std::size_t idx = trajBase + static_cast<std::size_t>(i);
        const Point3D& p = trajectories[idx];
        trajOut << lineToken << " " << (i + 1) << " "
                << p.x << " " << p.y << " " << p.z << "\n";
    }

    for (int i = 0; i < clampedPunctureCount; ++i) {
        const std::size_t idx = punctureBase + static_cast<std::size_t>(i);
        const PuncturePoint& p = punctures[idx];
        ipOut << lineToken << " " << p.step << " "
              << p.xyz.x << " " << p.xyz.y << " " << p.xyz.z << "\n";
        tpOut << lineToken << " " << p.step << " "
              << 0.0 << " " << p.thetaPsi.x << " " << p.thetaPsi.y << "\n";
    }
}

}  // namespace

void PoincareOutput::writeLineOutputs(const LineTraceResult& line,
                                      const std::string& outputDir,
                                      const std::string& divertorTag) const {
    std::filesystem::create_directories(outputDir);

    const std::string lineToken = formatLineToken(line.iline);
    const std::string ipPath = outputDir + "/ip_cxx." + divertorTag + "." + lineToken + ".txt";
    const std::string tpPath = outputDir + "/ip_cxx." + divertorTag + "." + lineToken + ".TP.txt";
    const std::string trajPath = outputDir + "/traj_cxx." + divertorTag + "." + lineToken + ".txt";

    std::ofstream ipOut(ipPath);
    std::ofstream tpOut(tpPath);
    std::ofstream trajOut(trajPath);
    if (!ipOut || !tpOut || !trajOut) {
        throw std::runtime_error("Failed to open output files for writing");
    }

    ipOut << "iline it ipx ipy ipz\n";
    tpOut << "iline it dummy theta psi\n";
    trajOut << "iline it x y z\n";

    ipOut << std::setprecision(16);
    tpOut << std::setprecision(16);
    trajOut << std::setprecision(16);

    writeLineToStreams(line, ipOut, tpOut, trajOut);
}

void PoincareOutput::writeLineOutputsFlat(double iline,
                                          std::size_t seedIndex,
                                          int maxTrajPerSeed,
                                          int maxPuncPerSeed,
                                          int trajCount,
                                          int punctureCount,
                                          const std::vector<Point3D>& trajectories,
                                          const std::vector<PuncturePoint>& punctures,
                                          const std::string& outputDir,
                                          const std::string& divertorTag) const {
    std::filesystem::create_directories(outputDir);

    const std::string lineToken = formatLineToken(iline);
    const std::string ipPath = outputDir + "/ip_cxx." + divertorTag + "." + lineToken + ".txt";
    const std::string tpPath = outputDir + "/ip_cxx." + divertorTag + "." + lineToken + ".TP.txt";
    const std::string trajPath = outputDir + "/traj_cxx." + divertorTag + "." + lineToken + ".txt";

    std::ofstream ipOut(ipPath);
    std::ofstream tpOut(tpPath);
    std::ofstream trajOut(trajPath);
    if (!ipOut || !tpOut || !trajOut) {
        throw std::runtime_error("Failed to open output files for writing");
    }

    ipOut << "iline it ipx ipy ipz\n";
    tpOut << "iline it dummy theta psi\n";
    trajOut << "iline it x y z\n";

    ipOut << std::setprecision(16);
    tpOut << std::setprecision(16);
    trajOut << std::setprecision(16);

    writeFlatLineToStreams(iline,
                           seedIndex,
                           maxTrajPerSeed,
                           maxPuncPerSeed,
                           trajCount,
                           punctureCount,
                           trajectories,
                           punctures,
                           ipOut,
                           tpOut,
                           trajOut);
}

void PoincareOutput::writeCombinedOutputs(const std::vector<LineTraceResult>& lines,
                                          const std::string& outputDir) const {
    std::filesystem::create_directories(outputDir);

    const std::string ipPath = outputDir + "/ip_cxx.txt";
    const std::string tpPath = outputDir + "/ip_cxx.TP.txt";
    const std::string trajPath = outputDir + "/traj_cxx.txt";

    std::ofstream ipOut(ipPath);
    std::ofstream tpOut(tpPath);
    std::ofstream trajOut(trajPath);
    if (!ipOut || !tpOut || !trajOut) {
        throw std::runtime_error("Failed to open combined output files for writing");
    }

    ipOut << "iline it ipx ipy ipz\n";
    tpOut << "iline it dummy theta psi\n";
    trajOut << "iline it x y z\n";

    ipOut << std::setprecision(16);
    tpOut << std::setprecision(16);
    trajOut << std::setprecision(16);

    for (const LineTraceResult& line : lines) {
        writeLineToStreams(line, ipOut, tpOut, trajOut);
    }
}

void PoincareOutput::writeCombinedOutputsFlat(const std::vector<double>& ilinePerSeed,
                                              const std::vector<int>& trajCountPerSeed,
                                              const std::vector<int>& punctureCountPerSeed,
                                              int maxTrajPerSeed,
                                              int maxPuncPerSeed,
                                              const std::vector<Point3D>& trajectories,
                                              const std::vector<PuncturePoint>& punctures,
                                              const std::string& outputDir) const {
    if (ilinePerSeed.size() != trajCountPerSeed.size() || ilinePerSeed.size() != punctureCountPerSeed.size()) {
        throw std::runtime_error("Per-seed metadata arrays must all have the same size");
    }

    std::filesystem::create_directories(outputDir);

    const std::string ipPath = outputDir + "/ip_cxx.txt";
    const std::string tpPath = outputDir + "/ip_cxx.TP.txt";
    const std::string trajPath = outputDir + "/traj_cxx.txt";

    std::ofstream ipOut(ipPath);
    std::ofstream tpOut(tpPath);
    std::ofstream trajOut(trajPath);
    if (!ipOut || !tpOut || !trajOut) {
        throw std::runtime_error("Failed to open combined output files for writing");
    }

    ipOut << "iline it ipx ipy ipz\n";
    tpOut << "iline it dummy theta psi\n";
    trajOut << "iline it x y z\n";

    ipOut << std::setprecision(16);
    tpOut << std::setprecision(16);
    trajOut << std::setprecision(16);

    for (std::size_t i = 0; i < ilinePerSeed.size(); ++i) {
        writeFlatLineToStreams(ilinePerSeed[i],
                               i,
                               maxTrajPerSeed,
                               maxPuncPerSeed,
                               trajCountPerSeed[i],
                               punctureCountPerSeed[i],
                               trajectories,
                               punctures,
                               ipOut,
                               tpOut,
                               trajOut);
    }
}
