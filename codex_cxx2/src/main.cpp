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
        << "  --lines LIST          comma-separated line numbers (default: all discovered)\n"
        << "  --nturns N            number of turns (default: 100)\n"
        << "  --direction DIR       +1 or -1 (default: 1)\n"
        << "  --np-max N            max punctures per line (default: 1250)\n"
        << "  --tol VALUE           max-abs tolerance (default: 1e-8)\n"
        << "  --no-compare          run tracing only, skip MATLAB comparisons\n"
        << "  --help                show this help\n";
}

std::vector<int> parseLinesCsv(const std::string& text) {
    std::vector<int> out;
    std::stringstream ss(text);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (token.empty()) {
            continue;
        }
        out.push_back(std::stoi(token));
    }
    return out;
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
    bool noCompare = false;

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
            config.linesFilter = parseLinesCsv(argv[++i]);
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
        if (arg == "--no-compare") {
            noCompare = true;
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
        if (noCompare) {
            if (config.linesFilter.empty()) {
                std::cerr << "Error: --no-compare requires --lines\n";
                return 1;
            }

            const bool doSingle = (config.divertorFilter == "all" || config.divertorFilter == "single");
            const bool doCirc = (config.divertorFilter == "all" || config.divertorFilter == "circ");

            if (!doSingle && !doCirc) {
                std::cerr << "Error: --divertor must be all|single|circ\n";
                return 1;
            }

            auto runDivertor = [&](const std::string& tag, const std::string& aparPath) {
                AparData data;
                data.load(aparPath);

                AparFieldModel model(data);
                FieldLineIntegrator integrator(model);
                PunctureDetector punctureDetector(model);
                PoincareOutput output;

                for (int line : config.linesFilter) {
                    LineTraceResult traced = integrator.traceLine(line, config.traceOptions);
                    punctureDetector.detect(traced, config.traceOptions.direction, config.traceOptions.npMax);
                    output.writeLineOutputs(traced, config.outputDir, tag);
                    std::cout << "Wrote " << tag << " line " << line
                              << ": traj=" << traced.trajectoryXYZ.size()
                              << ", punctures=" << traced.punctures.size() << "\n";
                }
            };

            if (doSingle) {
                runDivertor("single", config.aparSinglePath);
            }
            if (doCirc) {
                runDivertor("circ", config.aparCircPath);
            }

            std::cout << "\nTrace-only run complete.\n";
            return 0;
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
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
}
