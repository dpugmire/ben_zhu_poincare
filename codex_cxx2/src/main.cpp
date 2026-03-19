#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <limits>
#include <vector>

#include "AparData.h"
#include "AparFieldModel.h"
#include "FieldLineIntegrator.h"
#include "PoincareOutput.h"
#include "ValidationSuite.h"

namespace {

void printUsage(const char* prog) {
    std::cout
        << "Usage: " << prog << " [options]\n"
        << "Options:\n"
        << "  --reference-dir DIR   MATLAB reference directory (default: /Users/dpn/proj/bout++/ben_zhu_poincare/zperiod_5)\n"
        << "  --output-dir DIR      Output directory for generated files (default: ./outputs)\n"
        << "  --single-apar FILE    apar.single.nc path\n"
        << "  --circ-apar FILE      apar.circ.nc path\n"
        << "  --divertor TAG        all|single|circ (default: all)\n"
        << "  --lines LIST          comma-separated x-index values (double supported)\n"
        << "  --nlines X0 X1 N      generate N evenly-spaced lines from X0 to X1 (double supported; writes ip_cxx.txt/traj_cxx.txt)\n"
        << "  --direction DIR       +1 or -1 (default: 1)\n"
        << "  --np-max N            max punctures per line (default: 100)\n"
        << "  --max-steps N         max integration steps per line (default: 200 * np-max)\n"
        << "  --tol VALUE           max-abs tolerance (default: 1e-8)\n"
        << "  --compare             enable MATLAB comparison mode (default is trace-only)\n"
        << "  --help                show this help\n";
}

std::vector<double> parseLinesCsvDoubles(const std::string& text) {
    std::vector<double> out;
    std::stringstream ss(text);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (token.empty()) {
            continue;
        }
        out.push_back(std::stod(token));
    }
    return out;
}

std::vector<double> buildLineRange(double x0, double x1, int n) {
    std::vector<double> out;
    if (n <= 0) {
        return out;
    }
    out.reserve(static_cast<size_t>(n));

    if (n == 1) {
        out.push_back(x0);
        return out;
    }

    for (int i = 0; i < n; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(n - 1);
        out.push_back(x0 + t * (x1 - x0));
    }
    return out;
}

int computeMaxStateCount(const TraceOptions& options) {
    constexpr int kDefaultMaxStepsPerPuncture = 200;
    const long long derivedMaxStepsLL = static_cast<long long>(kDefaultMaxStepsPerPuncture) *
                                        static_cast<long long>(std::max(1, options.npMax));
    const int derivedMaxSteps = (derivedMaxStepsLL > static_cast<long long>(std::numeric_limits<int>::max()))
        ? std::numeric_limits<int>::max()
        : static_cast<int>(std::max(1LL, derivedMaxStepsLL));

    const int maxStepsCap = (options.maxSteps > 0) ? options.maxSteps : derivedMaxSteps;
    const int clampedSteps = std::max(1, maxStepsCap);
    return (clampedSteps >= std::numeric_limits<int>::max())
        ? std::numeric_limits<int>::max()
        : (clampedSteps + 1);
}

bool isIntegerVal(double v) {
    return std::isfinite(v) && std::fabs(v - std::round(v)) < 1.0e-12;
}

enum class LineSpecMode {
    none,
    csv,
    range
};

std::vector<Point3D> buildSeedPoints(const std::vector<double>& requestedLines, const AparData& data) {
    std::vector<Point3D> seeds;
    seeds.reserve(requestedLines.size());

    const double ySeed = static_cast<double>(data.jyomp + 1);
    const double zSeed = data.ziarray.empty() ? 1.0 : data.ziarray.front();
    for (double line : requestedLines) {
        Point3D seedInd;
        seedInd.x = line;
        seedInd.y = ySeed;
        seedInd.z = zSeed;
        seeds.push_back(seedInd);
    }
    return seeds;
}

void runDivertorTrace(const std::string& tag,
                      const std::string& aparPath,
                      const std::vector<double>& requestedLines,
                      const ValidationConfig& config,
                      bool writePerLineOutput,
                      PoincareOutput& output,
                      std::vector<Point3D>& seedsAll,
                      std::vector<double>& ilinePerSeed,
                      std::vector<int>& endRegionPerSeed,
                      std::vector<double>& connectionLengthPerSeed,
                      std::vector<int>& stateCountPerSeed,
                      std::vector<int>& trajCountPerSeed,
                      std::vector<int>& punctureCountPerSeed,
                      std::vector<TrajectoryState>& states,
                      std::vector<Point3D>& trajectories,
                      std::vector<PuncturePoint>& punctures,
                      int maxStatesPerSeed,
                      int maxTrajPerSeed,
                      int maxPuncPerSeed,
                      size_t seedIndexOffset) {
    AparData data;
    data.load(aparPath);

    AparFieldModel model(data);
    FieldLineIntegrator integrator(model, config.traceOptions);
    const std::vector<Point3D> seeds = buildSeedPoints(requestedLines, data);
    if (integrator.maxStatesPerSeed() != maxStatesPerSeed ||
        integrator.maxTrajPerSeed() != maxTrajPerSeed ||
        integrator.maxPuncPerSeed() != maxPuncPerSeed) {
        throw std::runtime_error("Integrator max-per-seed caps differ from preallocated array caps");
    }
    if (seedIndexOffset + seeds.size() > ilinePerSeed.size()) {
        throw std::runtime_error("Per-seed output arrays are smaller than number of seeds");
    }

    for (size_t i = 0; i < seeds.size(); ++i) {
        const Point3D& seedInd = seeds[i];
        const size_t seedIndex = seedIndexOffset + i;
        seedsAll[seedIndex] = seedInd;
        integrator.traceLine(seedInd,
                             seedIndex,
                             states,
                             trajectories,
                             punctures,
                             stateCountPerSeed[seedIndex],
                             trajCountPerSeed[seedIndex],
                             punctureCountPerSeed[seedIndex],
                             endRegionPerSeed[seedIndex],
                             connectionLengthPerSeed[seedIndex],
                             ilinePerSeed[seedIndex]);
        std::cout << "Traced " << tag << " line " << seedInd.x
                  << ": traj=" << trajCountPerSeed[seedIndex]
                  << ", punctures=" << punctureCountPerSeed[seedIndex] << "\n";
    }

    if (writePerLineOutput) {
        for (size_t i = 0; i < seeds.size(); ++i) {
            const size_t seedIndex = seedIndexOffset + i;
            output.writeLineOutputsFlat(ilinePerSeed[seedIndex],
                                        seedIndex,
                                        maxTrajPerSeed,
                                        maxPuncPerSeed,
                                        trajCountPerSeed[seedIndex],
                                        punctureCountPerSeed[seedIndex],
                                        trajectories,
                                        punctures,
                                        config.outputDir,
                                        tag);
        }
        std::cout << "Wrote per-line outputs for " << tag
                  << " (" << seeds.size() << " lines)\n";
    }
}

}  // namespace

int main(int argc, char** argv) {
    ValidationConfig config;
    config.referenceDir = "/Users/dpn/proj/bout++/ben_zhu_poincare/zperiod_5";
    config.outputDir = "./outputs";
    config.aparSinglePath = "/Users/dpn/proj/bout++/ben_zhu_poincare/apar.single.nc";
    config.aparCircPath = "/Users/dpn/proj/bout++/ben_zhu_poincare/apar.circ.nc";
    config.traceOptions.direction = 1;
    config.traceOptions.npMax = 100;
    config.tolerance = 1.0e-8;

    bool doCompare = false;
    LineSpecMode lineSpecMode = LineSpecMode::none;
    std::vector<double> requestedLines;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help") {
            printUsage(argv[0]);
            return 0;
        }
        if (arg == "--reference-dir" && i + 1 < argc) {
            config.referenceDir = argv[++i];
            continue;
        }
        if (arg == "--output-dir" && i + 1 < argc) {
            config.outputDir = argv[++i];
            continue;
        }
        if (arg == "--single-apar" && i + 1 < argc) {
            config.aparSinglePath = argv[++i];
            continue;
        }
        if (arg == "--circ-apar" && i + 1 < argc) {
            config.aparCircPath = argv[++i];
            continue;
        }
        if (arg == "--divertor" && i + 1 < argc) {
            config.divertorFilter = argv[++i];
            continue;
        }
        if (arg == "--lines" && i + 1 < argc) {
            if (lineSpecMode == LineSpecMode::range) {
                std::cerr << "Error: cannot use both --lines and --nlines\n";
                return 1;
            }
            lineSpecMode = LineSpecMode::csv;
            requestedLines = parseLinesCsvDoubles(argv[++i]);
            continue;
        }
        if (arg == "--nlines" && i + 3 < argc) {
            if (lineSpecMode == LineSpecMode::csv) {
                std::cerr << "Error: cannot use both --lines and --nlines\n";
                return 1;
            }
            lineSpecMode = LineSpecMode::range;
            const double x0 = std::stod(argv[++i]);
            const double x1 = std::stod(argv[++i]);
            const int n = std::atoi(argv[++i]);
            if (n <= 0) {
                std::cerr << "Error: --nlines requires N > 0\n";
                return 1;
            }
            requestedLines = buildLineRange(x0, x1, n);
            continue;
        }
        if (arg == "--direction" && i + 1 < argc) {
            config.traceOptions.direction = std::atoi(argv[++i]);
            continue;
        }
        if (arg == "--np-max" && i + 1 < argc) {
            config.traceOptions.npMax = std::atoi(argv[++i]);
            continue;
        }
        if (arg == "--max-steps" && i + 1 < argc) {
            config.traceOptions.maxSteps = std::atoi(argv[++i]);
            continue;
        }
        if (arg == "--tol" && i + 1 < argc) {
            config.tolerance = std::atof(argv[++i]);
            continue;
        }
        if (arg == "--compare") {
            doCompare = true;
            continue;
        }

        std::cerr << "Unknown or incomplete argument: " << arg << "\n";
        printUsage(argv[0]);
        return 1;
    }

    if (config.traceOptions.direction != 1 && config.traceOptions.direction != -1) {
        std::cerr << "Error: --direction must be 1 or -1\n";
        return 1;
    }

    if (config.traceOptions.npMax <= 0) {
        std::cerr << "Error: --np-max must be positive\n";
        return 1;
    }

    if (config.traceOptions.maxSteps < 0) {
        std::cerr << "Error: --max-steps must be non-negative\n";
        return 1;
    }

    try {
        if (doCompare) {
            if (lineSpecMode == LineSpecMode::range) {
                std::cerr << "Error: --nlines is only supported in trace-only mode (no --compare)\n";
                return 1;
            }

            if (lineSpecMode == LineSpecMode::csv) {
                config.linesFilter.clear();
                config.linesFilter.reserve(requestedLines.size());
                for (double line : requestedLines) {
                    if (!isIntegerVal(line)) {
                        std::cerr << "Error: --compare requires integer line indices; got " << line << "\n";
                        return 1;
                    }
                    config.linesFilter.push_back(static_cast<int>(std::llround(line)));
                }
            }

            ValidationSuite suite;
            const std::vector<ValidationResult> results = suite.run(config);

            int passed = 0;
            for (const auto& r : results) {
                if (r.pass) {
                    ++passed;
                }
            }

            std::cout << "\nValidation summary: " << passed << "/" << results.size() << " cases passed\n";
            return (passed == static_cast<int>(results.size())) ? 0 : 2;
        }

        // Default mode: trace only, no comparisons.
        if (requestedLines.empty()) {
            std::cerr << "Error: trace-only mode requires --lines or --nlines\n";
            return 1;
        }

        const bool doSingle = (config.divertorFilter == "all" || config.divertorFilter == "single");
        const bool doCirc = (config.divertorFilter == "all" || config.divertorFilter == "circ");
        const bool useCombinedOutput = (lineSpecMode == LineSpecMode::range);

        if (!doSingle && !doCirc) {
            std::cerr << "Error: --divertor must be all|single|circ\n";
            return 1;
        }

        PoincareOutput output;
        const size_t tracesPerDivertor = requestedLines.size();
        const size_t divertorCount = static_cast<size_t>((doSingle ? 1 : 0) + (doCirc ? 1 : 0));
        const size_t totalTraceCount = tracesPerDivertor * divertorCount;
        const int maxStatesPerSeed = computeMaxStateCount(config.traceOptions);
        const int maxTrajPerSeed = maxStatesPerSeed;
        const int maxPuncPerSeed = std::max(1, config.traceOptions.npMax);

        std::vector<Point3D> seedsAll(totalTraceCount);
        std::vector<double> ilinePerSeed(totalTraceCount, 0.0);
        std::vector<int> endRegionPerSeed(totalTraceCount, 0);
        std::vector<double> connectionLengthPerSeed(totalTraceCount, 0.0);
        std::vector<int> stateCountPerSeed(totalTraceCount, 0);
        std::vector<int> trajCountPerSeed(totalTraceCount, 0);
        std::vector<int> punctureCountPerSeed(totalTraceCount, 0);
        std::vector<TrajectoryState> states(totalTraceCount * static_cast<size_t>(maxStatesPerSeed));
        std::vector<Point3D> trajectories(totalTraceCount * static_cast<size_t>(maxTrajPerSeed));
        std::vector<PuncturePoint> punctures(totalTraceCount * static_cast<size_t>(maxPuncPerSeed));

        size_t seedIndexOffset = 0;

        if (doSingle) {
            runDivertorTrace("single",
                             config.aparSinglePath,
                             requestedLines,
                             config,
                             !useCombinedOutput,
                             output,
                             seedsAll,
                             ilinePerSeed,
                             endRegionPerSeed,
                             connectionLengthPerSeed,
                             stateCountPerSeed,
                             trajCountPerSeed,
                             punctureCountPerSeed,
                             states,
                             trajectories,
                             punctures,
                             maxStatesPerSeed,
                             maxTrajPerSeed,
                             maxPuncPerSeed,
                             seedIndexOffset);
            seedIndexOffset += tracesPerDivertor;
        }
        if (doCirc) {
            runDivertorTrace("circ",
                             config.aparCircPath,
                             requestedLines,
                             config,
                             !useCombinedOutput,
                             output,
                             seedsAll,
                             ilinePerSeed,
                             endRegionPerSeed,
                             connectionLengthPerSeed,
                             stateCountPerSeed,
                             trajCountPerSeed,
                             punctureCountPerSeed,
                             states,
                             trajectories,
                             punctures,
                             maxStatesPerSeed,
                             maxTrajPerSeed,
                             maxPuncPerSeed,
                             seedIndexOffset);
            seedIndexOffset += tracesPerDivertor;
        }

        if (useCombinedOutput) {
            output.writeCombinedOutputsFlat(ilinePerSeed,
                                            trajCountPerSeed,
                                            punctureCountPerSeed,
                                            maxTrajPerSeed,
                                            maxPuncPerSeed,
                                            trajectories,
                                            punctures,
                                            config.outputDir);
            std::cout << "Wrote combined outputs: " << config.outputDir
                      << "/ip_cxx.txt, " << config.outputDir << "/ip_cxx.TP.txt"
                      << " and " << config.outputDir << "/traj_cxx.txt\n";
        }

        std::cout << "\nTrace-only run complete.\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
}
