#include "PoincareOutput.h"

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

void writeBatchLineToStreams(const PackedLineTraceBatch& batch,
                             size_t lineIndex,
                             std::ofstream& ipOut,
                             std::ofstream& tpOut,
                             std::ofstream& trajOut) {
    const std::string lineToken = formatLineToken(batch.iline[lineIndex]);
    const size_t stateBase = batch.stateOffset(lineIndex);
    const size_t punctureBase = batch.punctureOffset(lineIndex);

    const int trajCount = batch.trajectoryCount[lineIndex];
    for (int i = 0; i < trajCount; ++i) {
        const size_t idx = stateBase + static_cast<size_t>(i);
        if (!batch.trajectorySet[idx]) {
            continue;
        }
        const Point3D& p = batch.trajectoryXYZ[idx];
        trajOut << lineToken << " " << (i + 1) << " "
                << p.x << " " << p.y << " " << p.z << "\n";
    }

    const int punctureCount = batch.punctureCount[lineIndex];
    for (int i = 0; i < punctureCount; ++i) {
        const size_t idx = punctureBase + static_cast<size_t>(i);
        if (!batch.punctureSet[idx]) {
            continue;
        }
        const PuncturePoint& p = batch.punctures[idx];
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

void PoincareOutput::writeCombinedOutputs(const PackedLineTraceBatch& batch,
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

    for (size_t i = 0; i < batch.lineCount(); ++i) {
        writeBatchLineToStreams(batch, i, ipOut, tpOut, trajOut);
    }
}
