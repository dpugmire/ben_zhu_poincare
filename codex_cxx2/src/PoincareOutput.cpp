#include "PoincareOutput.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

void PoincareOutput::writeLineOutputs(const LineTraceResult& line,
                                      const std::string& outputDir,
                                      const std::string& divertorTag) const {
    std::filesystem::create_directories(outputDir);

    const std::string ipPath = outputDir + "/ip_cxx." + divertorTag + "." + std::to_string(line.iline) + ".txt";
    const std::string trajPath = outputDir + "/traj_cxx." + divertorTag + "." + std::to_string(line.iline) + ".txt";

    std::ofstream ipOut(ipPath);
    std::ofstream trajOut(trajPath);
    if (!ipOut || !trajOut) {
        throw std::runtime_error("Failed to open output files for writing");
    }

    ipOut << "iline it ipx ipy ipz\n";
    trajOut << "iline it x y z\n";

    ipOut << std::setprecision(16);
    trajOut << std::setprecision(16);

    for (size_t i = 0; i < line.trajectoryXYZ.size(); ++i) {
        const Point3D& p = line.trajectoryXYZ[i];
        trajOut << line.iline << " " << (i + 1) << " "
                << p.x << " " << p.y << " " << p.z << "\n";
    }

    for (const PuncturePoint& p : line.punctures) {
        ipOut << line.iline << " " << p.step << " "
              << p.xyz.x << " " << p.xyz.y << " " << p.xyz.z << "\n";
    }
}
