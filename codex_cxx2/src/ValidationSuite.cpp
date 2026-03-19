#include "ValidationSuite.h"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <regex>
#include <set>
#include <stdexcept>
#include <tuple>

#include "AparData.h"
#include "AparFieldModel.h"
#include "FieldLineIntegrator.h"
#include "MatlabComparator.h"
#include "PoincareOutput.h"

namespace {

struct CaseAccumulator {
    std::string ipPath;
    std::string trajPath;
};

bool lineIsSelected(const std::vector<int>& filter, int line) {
    if (filter.empty()) {
        return true;
    }
    for (int v : filter) {
        if (v == line) {
            return true;
        }
    }
    return false;
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

bool compareValidationCase(const ValidationCase& a, const ValidationCase& b) {
    return std::tie(a.divertorTag, a.line) < std::tie(b.divertorTag, b.line);
}

}  // namespace

std::vector<ValidationCase> ValidationSuite::discoverCases(const ValidationConfig& config) {
    namespace fs = std::filesystem;

    std::map<std::pair<std::string, int>, CaseAccumulator> found;
    const std::regex pattern(R"((ip|traj)_matlab\.(single|circ)\.(\d+)\.txt)");

    for (const auto& entry : fs::directory_iterator(config.referenceDir)) {
        if (!entry.is_regular_file()) {
            continue;
        }

        const std::string name = entry.path().filename().string();
        std::smatch m;
        if (!std::regex_match(name, m, pattern)) {
            continue;
        }

        const std::string type = m[1].str();
        const std::string tag = m[2].str();
        const int line = std::stoi(m[3].str());

        auto& acc = found[{tag, line}];
        if (type == "ip") {
            acc.ipPath = entry.path().string();
        } else {
            acc.trajPath = entry.path().string();
        }
    }

    std::vector<ValidationCase> out;
    out.reserve(found.size());

    for (const auto& kv : found) {
        const std::string& tag = kv.first.first;
        const int line = kv.first.second;
        const CaseAccumulator& acc = kv.second;

        if (acc.ipPath.empty() || acc.trajPath.empty()) {
            continue;
        }

        if (config.divertorFilter != "all" && config.divertorFilter != tag) {
            continue;
        }

        if (!lineIsSelected(config.linesFilter, line)) {
            continue;
        }

        ValidationCase testCase;
        testCase.divertorTag = tag;
        testCase.line = line;
        testCase.refIpPath = acc.ipPath;
        testCase.refTrajPath = acc.trajPath;
        if (tag == "single") {
            testCase.aparFile = config.aparSinglePath;
        } else if (tag == "circ") {
            testCase.aparFile = config.aparCircPath;
        } else {
            continue;
        }

        out.push_back(testCase);
    }

    std::sort(out.begin(), out.end(), compareValidationCase);

    return out;
}

std::vector<ValidationResult> ValidationSuite::run(const ValidationConfig& config) const {
    std::vector<ValidationCase> cases = discoverCases(config);
    if (cases.empty()) {
        throw std::runtime_error("No validation cases discovered");
    }

    std::map<std::string, std::shared_ptr<AparData>> dataCache;

    PoincareOutput output;
    MatlabComparator comparator;

    std::vector<ValidationResult> results;
    results.reserve(cases.size());
    const int maxStateCount = computeMaxStateCount(config.traceOptions);
    LineTraceResult lineResult;
    lineResult.states.reserve(static_cast<size_t>(std::max(0, maxStateCount)));
    lineResult.trajectoryXYZ.reserve(static_cast<size_t>(std::max(0, maxStateCount)));
    lineResult.punctures.reserve(static_cast<size_t>(std::max(0, config.traceOptions.npMax)));

    for (const ValidationCase& testCase : cases) {
        auto& slot = dataCache[testCase.divertorTag];
        if (!slot) {
            auto data = std::make_shared<AparData>();
            data->load(testCase.aparFile);
            slot = std::move(data);
        }

        AparFieldModel model(*slot);
        FieldLineIntegrator integrator(model, config.traceOptions);

        std::cout << "Running case: divertor=" << testCase.divertorTag
                  << ", line=" << testCase.line << std::endl;

        Point3D seedInd;
        seedInd.x = static_cast<double>(testCase.line);
        seedInd.y = static_cast<double>(slot->jyomp + 1);
        seedInd.z = slot->ziarray.empty() ? 1.0 : slot->ziarray.front();
        integrator.traceLine(seedInd, lineResult);

        output.writeLineOutputs(lineResult, config.outputDir, testCase.divertorTag);

        const std::string genIp = config.outputDir + "/ip_cxx." + testCase.divertorTag + "." +
                                  std::to_string(testCase.line) + ".txt";
        const std::string genTraj = config.outputDir + "/traj_cxx." + testCase.divertorTag + "." +
                                    std::to_string(testCase.line) + ".txt";

        ValidationResult vr;
        vr.testCase = testCase;
        vr.ipSummary = comparator.compareFiles(
            genIp,
            testCase.refIpPath,
            "ip " + testCase.divertorTag + " line " + std::to_string(testCase.line),
            config.tolerance);
        vr.trajSummary = comparator.compareFiles(
            genTraj,
            testCase.refTrajPath,
            "traj " + testCase.divertorTag + " line " + std::to_string(testCase.line),
            config.tolerance);
        vr.pass = vr.ipSummary.pass && vr.trajSummary.pass;

        std::cout << "  ip:   rows(gen/ref)=" << vr.ipSummary.aRows << "/" << vr.ipSummary.bRows
                  << ", maxAbs=" << vr.ipSummary.maxAbsError
                  << ", L2=" << vr.ipSummary.l2Error
                  << ", pass=" << (vr.ipSummary.pass ? "yes" : "no") << std::endl;

        std::cout << "  traj: rows(gen/ref)=" << vr.trajSummary.aRows << "/" << vr.trajSummary.bRows
                  << ", maxAbs=" << vr.trajSummary.maxAbsError
                  << ", L2=" << vr.trajSummary.l2Error
                  << ", pass=" << (vr.trajSummary.pass ? "yes" : "no") << std::endl;

        std::cout << "  case result: " << (vr.pass ? "PASS" : "FAIL") << std::endl;

        results.push_back(vr);
    }

    return results;
}
