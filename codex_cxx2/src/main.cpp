#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "AparData.h"
#include "AparFieldModel.h"
#include "FieldLineIntegrator.h"
#include "PoincareOutput.h"
#include "PunctureDetector.h"
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
        << "  --nturns N            number of turns (default: 100)\n"
        << "  --direction DIR       +1 or -1 (default: 1)\n"
        << "  --np-max N            max punctures per line (default: 1250)\n"
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

bool isIntegerVal(double v) {
    return std::isfinite(v) && std::fabs(v - std::round(v)) < 1.0e-12;
}

enum class LineSpecMode {
    none,
    csv,
    range
};

void runDivertorTrace(const std::string& tag,
                      const std::string& aparPath,
                      const std::vector<double>& requestedLines,
                      const ValidationConfig& config,
                      bool useCombinedOutput,
                      PoincareOutput& output,
                      std::vector<LineTraceResult>& combinedLines) {
    AparData data;
    data.load(aparPath);

    AparFieldModel model(data);
    FieldLineIntegrator integrator(model);
    PunctureDetector punctureDetector(model);

    for (double line : requestedLines) {
        LineTraceResult traced = integrator.traceLine(line, config.traceOptions);
        punctureDetector.detect(traced, config.traceOptions.direction, config.traceOptions.npMax);
        const size_t trajCount = traced.trajectoryXYZ.size();
        const size_t punctureCount = traced.punctures.size();
        if (useCombinedOutput) {
            combinedLines.push_back(traced);
        } else {
            output.writeLineOutputs(traced, config.outputDir, tag);
        }
        std::cout << "Wrote " << tag << " line " << line
                  << ": traj=" << trajCount
                  << ", punctures=" << punctureCount << "\n";
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
    config.traceOptions.nturns = 100;
    config.traceOptions.npMax = 1250;
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
        if (arg == "--nturns" && i + 1 < argc) {
            config.traceOptions.nturns = std::atoi(argv[++i]);
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

    if (config.traceOptions.nturns <= 0) {
        std::cerr << "Error: --nturns must be positive\n";
        return 1;
    }

    if (config.traceOptions.npMax <= 0) {
        std::cerr << "Error: --np-max must be positive\n";
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
        std::vector<LineTraceResult> combinedLines;

        if (doSingle) {
            runDivertorTrace("single",
                             config.aparSinglePath,
                             requestedLines,
                             config,
                             useCombinedOutput,
                             output,
                             combinedLines);
        }
        if (doCirc) {
            runDivertorTrace("circ",
                             config.aparCircPath,
                             requestedLines,
                             config,
                             useCombinedOutput,
                             output,
                             combinedLines);
        }

        if (useCombinedOutput) {
            output.writeCombinedOutputs(combinedLines, config.outputDir);
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
