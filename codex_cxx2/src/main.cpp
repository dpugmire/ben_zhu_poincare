#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(CODEX_USE_MPI)
#include <mpi.h>
#endif

#if defined(CODEX_USE_VISKORES)
#include <viskores/cont/RuntimeDeviceInformation.h>
#include <viskores/cont/RuntimeDeviceTracker.h>
#endif

#include "AparData.h"
#include "AparFieldModel.h"
#include "FieldLineIntegrator.h"
#include "PoincareOutput.h"
#include "ValidationSuite.h"
#if defined(CODEX_USE_VISKORES)
#include "ViskoresFieldLineTracer.h"
#endif

namespace
{

#if defined(CODEX_USE_MPI)
class MpiRuntime
{
public:
  MpiRuntime(int* argc, char*** argv)
  {
    MPI_Init(argc, argv);
    initialized_ = true;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }

  ~MpiRuntime()
  {
    if (initialized_)
    {
      MPI_Finalize();
    }
  }

  void barrier() const { MPI_Barrier(MPI_COMM_WORLD); }

  int allreduceMaxInt(int localValue) const
  {
    int globalValue = localValue;
    MPI_Allreduce(&localValue, &globalValue, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    return globalValue;
  }

  double allreduceMaxDouble(double localValue) const
  {
    double globalValue = localValue;
    MPI_Allreduce(&localValue, &globalValue, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return globalValue;
  }

  double allreduceSumDouble(double localValue) const
  {
    double globalValue = localValue;
    MPI_Allreduce(&localValue, &globalValue, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return globalValue;
  }

  void abortAll(int code) const { MPI_Abort(MPI_COMM_WORLD, code); }

  int rank = 0;
  int size = 1;

private:
  bool initialized_ = false;
};
#else
class MpiRuntime
{
public:
  explicit MpiRuntime(int*, char***) {}

  void barrier() const {}

  int allreduceMaxInt(int localValue) const { return localValue; }

  double allreduceMaxDouble(double localValue) const { return localValue; }

  double allreduceSumDouble(double localValue) const { return localValue; }

  void abortAll(int code) const { std::exit(code); }

  int rank = 0;
  int size = 1;
};
#endif

struct TraceTask
{
  std::string divertorTag;
  Point3D seedInd;
};

void printUsage(const char* prog)
{
  std::cout << "Usage: " << prog << " [options]\n"
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
            << "  --trace-engine ENG    tracing backend: cpu|viskores (default: viskores)\n"
            << "  --viskores-device DEV viskores device: serial|openmp|kokkos (default: serial)\n"
            << "  --compare             enable MATLAB comparison mode (default is trace-only)\n"
            << "  --help                show this help\n";
}

bool hasHelpArgument(int argc, char** argv)
{
  for (int i = 1; i < argc; ++i)
  {
    const std::string arg = argv[i];
    if (arg == "--help" || arg == "-h")
    {
      return true;
    }
  }
  return false;
}

std::vector<double> parseLinesCsvDoubles(const std::string& text)
{
  std::vector<double> out;
  std::stringstream ss(text);
  std::string token;
  while (std::getline(ss, token, ','))
  {
    if (token.empty())
    {
      continue;
    }
    out.push_back(std::stod(token));
  }
  return out;
}

std::vector<double> buildLineRange(double x0, double x1, int n)
{
  std::vector<double> out;
  if (n <= 0)
  {
    return out;
  }
  out.reserve(static_cast<size_t>(n));

  if (n == 1)
  {
    out.push_back(x0);
    return out;
  }

  for (int i = 0; i < n; ++i)
  {
    const double t = static_cast<double>(i) / static_cast<double>(n - 1);
    out.push_back(x0 + t * (x1 - x0));
  }
  return out;
}

int computeMaxStateCount(const TraceOptions& options)
{
  constexpr int kDefaultMaxStepsPerPuncture = 200;
  const long long derivedMaxStepsLL = static_cast<long long>(kDefaultMaxStepsPerPuncture) * static_cast<long long>(std::max(1, options.npMax));
  const int derivedMaxSteps = (derivedMaxStepsLL > static_cast<long long>(std::numeric_limits<int>::max()))
    ? std::numeric_limits<int>::max()
    : static_cast<int>(std::max(1LL, derivedMaxStepsLL));

  const int maxStepsCap = (options.maxSteps > 0) ? options.maxSteps : derivedMaxSteps;
  const int clampedSteps = std::max(1, maxStepsCap);
  return (clampedSteps >= std::numeric_limits<int>::max()) ? std::numeric_limits<int>::max() : (clampedSteps + 1);
}

bool isIntegerVal(double v)
{
  return std::isfinite(v) && std::fabs(v - std::round(v)) < 1.0e-12;
}

enum class LineSpecMode
{
  none,
  csv,
  range
};

enum class TraceEngine
{
  Cpu,
  Viskores
};

bool isFatalTraceStatus(TraceStatus status)
{
  return status == TraceStatus::InvalidConfiguration || status == TraceStatus::OutputTooSmall;
}

std::string toLowerAscii(std::string text)
{
  std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return text;
}

bool isSupportedViskoresDeviceChoice(const std::string& valueLower)
{
  return valueLower == "serial" || valueLower == "openmp" || valueLower == "kokkos";
}

bool isSupportedTraceEngineChoice(const std::string& valueLower)
{
  return valueLower == "cpu" || valueLower == "viskores";
}

TraceEngine parseTraceEngineChoice(const std::string& valueLower)
{
  return (valueLower == "viskores") ? TraceEngine::Viskores : TraceEngine::Cpu;
}

const char* traceEngineName(TraceEngine traceEngine)
{
  return (traceEngine == TraceEngine::Viskores) ? "viskores" : "cpu";
}

void printTraceBackendSelection(const MpiRuntime& mpi, TraceEngine traceEngine, const std::string& viskoresDevice)
{
  const bool viskoresActive = (traceEngine == TraceEngine::Viskores);

  mpi.barrier();
  for (int printer = 0; printer < mpi.size; ++printer)
  {
    if (mpi.rank == printer)
    {
      std::cout << "Rank " << mpi.rank << " runtime selection: trace-engine=" << traceEngineName(traceEngine)
                << ", viskores-device=" << viskoresDevice << ", viskores-active=" << (viskoresActive ? "yes" : "no");
      if (!viskoresActive)
      {
        std::cout << " (device ignored unless --trace-engine viskores)";
      }
      std::cout << std::endl;
    }
    mpi.barrier();
  }
}


using SteadyClock = std::chrono::steady_clock;

struct TimingSummary
{
  double avgSeconds = 0.0;
  double maxSeconds = 0.0;
};

double elapsedSeconds(const SteadyClock::time_point& start, const SteadyClock::time_point& end)
{
  return std::chrono::duration<double>(end - start).count();
}

TimingSummary summarizeTiming(const MpiRuntime& mpi, double localSeconds)
{
  TimingSummary summary;
  summary.avgSeconds = mpi.allreduceSumDouble(localSeconds) / static_cast<double>(std::max(1, mpi.size));
  summary.maxSeconds = mpi.allreduceMaxDouble(localSeconds);
  return summary;
}

void printTimingLine(const char* label, const TimingSummary& summary)
{
  std::cout << "  " << label << ": avg=" << summary.avgSeconds << " s, max=" << summary.maxSeconds << " s\n";
}

void printTraceTimingSummary(const MpiRuntime& mpi,
                             double localLoadSeconds,
                             double localTraceSeconds,
                             double localOutputSeconds,
                             double localTotalSeconds)
{
  const TimingSummary load = summarizeTiming(mpi, localLoadSeconds);
  const TimingSummary trace = summarizeTiming(mpi, localTraceSeconds);
  const TimingSummary output = summarizeTiming(mpi, localOutputSeconds);
  const TimingSummary total = summarizeTiming(mpi, localTotalSeconds);

  if (mpi.rank != 0)
  {
    return;
  }

  const std::ios::fmtflags oldFlags = std::cout.flags();
  const std::streamsize oldPrecision = std::cout.precision();
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "\nTiming summary (trace-only, seconds)\n";
  printTimingLine("load", load);
  printTimingLine("trace", trace);
  printTimingLine("output", output);
  printTimingLine("total", total);
  std::cout.flags(oldFlags);
  std::cout.precision(oldPrecision);
}

void printCompareTimingSummary(const MpiRuntime& mpi, double localCompareSeconds, double localTotalSeconds)
{
  const TimingSummary compare = summarizeTiming(mpi, localCompareSeconds);
  const TimingSummary total = summarizeTiming(mpi, localTotalSeconds);

  if (mpi.rank != 0)
  {
    return;
  }

  const std::ios::fmtflags oldFlags = std::cout.flags();
  const std::streamsize oldPrecision = std::cout.precision();
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "\nTiming summary (compare, seconds)\n";
  printTimingLine("compare", compare);
  printTimingLine("total", total);
  std::cout.flags(oldFlags);
  std::cout.precision(oldPrecision);
}

#if defined(CODEX_USE_VISKORES)
const char* toViskoresDeviceName(const std::string& valueLower)
{
  if (valueLower == "serial")
  {
    return "Serial";
  }
  if (valueLower == "openmp")
  {
    return "OpenMP";
  }
  if (valueLower == "kokkos")
  {
    return "Kokkos";
  }
  return "";
}

bool configureViskoresDevice(const std::string& valueLower, std::string& errorMsg)
{
  const std::string viskoresName = toViskoresDeviceName(valueLower);
  if (viskoresName.empty())
  {
    errorMsg = "Unsupported --viskores-device value '" + valueLower + "' (allowed: serial|openmp|kokkos)";
    return false;
  }

  viskores::cont::RuntimeDeviceInformation info;
  const viskores::cont::DeviceAdapterId deviceId = info.GetId(viskoresName);
  if (!deviceId.IsValueValid())
  {
    errorMsg = "Viskores device name '" + viskoresName + "' is not recognized";
    return false;
  }
  if (!info.Exists(deviceId))
  {
    errorMsg = "Viskores device '" + valueLower + "' is not available in this build/runtime";
    return false;
  }

  viskores::cont::GetRuntimeDeviceTracker().ForceDevice(deviceId);
  return true;
}
#endif

std::vector<Point3D> buildSeedPoints(const std::vector<double>& requestedLines, const AparData& data)
{
  std::vector<Point3D> seeds;
  seeds.reserve(requestedLines.size());

  const double ySeed = static_cast<double>(data.jyomp + 1);
  const double zSeed = data.ziarray.empty() ? 1.0 : data.ziarray.front();
  for (double line : requestedLines)
  {
    Point3D seedInd;
    seedInd.x = line;
    seedInd.y = ySeed;
    seedInd.z = zSeed;
    seeds.push_back(seedInd);
  }
  return seeds;
}

std::vector<TraceTask> buildAllTasks(const std::vector<double>& requestedLines,
                                     bool doSingle,
                                     bool doCirc,
                                     const AparData* singleData,
                                     const AparData* circData)
{
  std::vector<TraceTask> allTasks;
  const size_t divertorCount = static_cast<size_t>((doSingle ? 1 : 0) + (doCirc ? 1 : 0));
  allTasks.reserve(requestedLines.size() * divertorCount);

  if (doSingle && singleData != nullptr)
  {
    const std::vector<Point3D> singleSeeds = buildSeedPoints(requestedLines, *singleData);
    for (const Point3D& seed : singleSeeds)
    {
      TraceTask t;
      t.divertorTag = "single";
      t.seedInd = seed;
      allTasks.push_back(t);
    }
  }

  if (doCirc && circData != nullptr)
  {
    const std::vector<Point3D> circSeeds = buildSeedPoints(requestedLines, *circData);
    for (const Point3D& seed : circSeeds)
    {
      TraceTask t;
      t.divertorTag = "circ";
      t.seedInd = seed;
      allTasks.push_back(t);
    }
  }

  return allTasks;
}

size_t partitionBegin(size_t totalTasks, int rank, int nranks)
{
  return (totalTasks * static_cast<size_t>(rank)) / static_cast<size_t>(nranks);
}

size_t partitionEnd(size_t totalTasks, int rank, int nranks)
{
  return (totalTasks * static_cast<size_t>(rank + 1)) / static_cast<size_t>(nranks);
}

void traceLocalTasksForDivertor(const std::string& tag,
                                const AparData& data,
                                const ValidationConfig& config,
                                const std::vector<TraceTask>& localTasks,
                                int rank,
                                int maxStatesPerSeed,
                                int maxTrajPerSeed,
                                int maxPuncPerSeed,
                                std::vector<double>& ilinePerSeed,
                                std::vector<int>& endRegionPerSeed,
                                std::vector<double>& connectionLengthPerSeed,
                                std::vector<int>& stateCountPerSeed,
                                std::vector<int>& trajCountPerSeed,
                                std::vector<int>& punctureCountPerSeed,
                                std::vector<TrajectoryState>& states,
                                std::vector<Point3D>& trajectories,
                                std::vector<PuncturePoint>& punctures,
                                std::vector<std::uint8_t>& punctureValid,
                                bool& localFatal,
                                std::string& localFatalMsg)
{
  AparFieldModel model(data);
  FieldLineIntegrator integrator(model, config.traceOptions);

  if (integrator.maxStatesPerSeed() != maxStatesPerSeed || integrator.maxTrajPerSeed() != maxTrajPerSeed ||
      integrator.maxPuncPerSeed() != maxPuncPerSeed)
  {
    localFatal = true;
    localFatalMsg = "Integrator max-per-seed caps differ from preallocated array caps";
    return;
  }

  TraceOutputViews outputViews = makeTraceOutputViews(states, trajectories, punctures, &punctureValid);

  for (size_t localIdx = 0; localIdx < localTasks.size(); ++localIdx)
  {
    if (localTasks[localIdx].divertorTag != tag)
    {
      continue;
    }

    const TraceStatus traceStatus = integrator.traceLine(localTasks[localIdx].seedInd,
                                                         localIdx,
                                                         outputViews,
                                                         stateCountPerSeed[localIdx],
                                                         trajCountPerSeed[localIdx],
                                                         punctureCountPerSeed[localIdx],
                                                         endRegionPerSeed[localIdx],
                                                         connectionLengthPerSeed[localIdx],
                                                         ilinePerSeed[localIdx]);
    if (isFatalTraceStatus(traceStatus))
    {
      localFatal = true;
      localFatalMsg = "traceLine failed on rank " + std::to_string(rank) + " for local seed " + std::to_string(localIdx) + " with status " +
        traceStatusName(traceStatus);
      return;
    }

    std::cout << "Rank " << rank << " traced " << tag << " line " << localTasks[localIdx].seedInd.x << ": traj=" << trajCountPerSeed[localIdx]
              << ", punctures=" << punctureCountPerSeed[localIdx];
    if (traceStatus != TraceStatus::Ok)
    {
      std::cout << ", status=" << traceStatusName(traceStatus);
    }
    std::cout << std::endl;
  }
}

#if defined(CODEX_USE_VISKORES)
void traceLocalTasksForDivertorViskores(const std::string& tag,
                                        const AparData& data,
                                        const ValidationConfig& config,
                                        const std::vector<TraceTask>& localTasks,
                                        int rank,
                                        int maxStatesPerSeed,
                                        int maxTrajPerSeed,
                                        int maxPuncPerSeed,
                                        std::vector<double>& ilinePerSeed,
                                        std::vector<int>& endRegionPerSeed,
                                        std::vector<double>& connectionLengthPerSeed,
                                        std::vector<int>& stateCountPerSeed,
                                        std::vector<int>& trajCountPerSeed,
                                        std::vector<int>& punctureCountPerSeed,
                                        std::vector<TrajectoryState>& states,
                                        std::vector<Point3D>& trajectories,
                                        std::vector<PuncturePoint>& punctures,
                                        std::vector<std::uint8_t>& punctureValid,
                                        bool& localFatal,
                                        std::string& localFatalMsg)
{
  ViskoresFieldLineTracer tracer(data, config.traceOptions);

  if (tracer.maxStatesPerSeed() != maxStatesPerSeed || tracer.maxTrajPerSeed() != maxTrajPerSeed ||
      tracer.maxPuncPerSeed() != maxPuncPerSeed)
  {
    localFatal = true;
    localFatalMsg = "Viskores tracer max-per-seed caps differ from preallocated array caps";
    return;
  }

  std::vector<Point3D> batchSeeds;
  std::vector<CodeXId> batchGlobalSeedIndices;
  batchSeeds.reserve(localTasks.size());
  batchGlobalSeedIndices.reserve(localTasks.size());

  for (std::size_t localIdx = 0; localIdx < localTasks.size(); ++localIdx)
  {
    if (localTasks[localIdx].divertorTag != tag)
    {
      continue;
    }
    batchSeeds.push_back(localTasks[localIdx].seedInd);
    batchGlobalSeedIndices.push_back(static_cast<CodeXId>(localIdx));
  }

  if (batchSeeds.empty())
  {
    return;
  }

  std::vector<TraceStatus> traceStatuses;
  tracer.traceLines(batchSeeds,
                    batchGlobalSeedIndices,
                    makeTraceOutputViews(states, trajectories, punctures, &punctureValid),
                    ilinePerSeed,
                    endRegionPerSeed,
                    connectionLengthPerSeed,
                    stateCountPerSeed,
                    trajCountPerSeed,
                    punctureCountPerSeed,
                    traceStatuses);

  for (std::size_t batchIdx = 0; batchIdx < batchSeeds.size(); ++batchIdx)
  {
    const std::size_t localIdx = static_cast<std::size_t>(batchGlobalSeedIndices[batchIdx]);
    const TraceStatus traceStatus = traceStatuses[batchIdx];
    if (isFatalTraceStatus(traceStatus))
    {
      localFatal = true;
      localFatalMsg = "viskores trace failed on rank " + std::to_string(rank) + " for local seed " + std::to_string(localIdx) +
        " with status " + traceStatusName(traceStatus);
      return;
    }

    std::cout << "Rank " << rank << " traced " << tag << " line " << localTasks[localIdx].seedInd.x << ": traj=" << trajCountPerSeed[localIdx]
              << ", punctures=" << punctureCountPerSeed[localIdx];
    if (traceStatus != TraceStatus::Ok)
    {
      std::cout << ", status=" << traceStatusName(traceStatus);
    }
    std::cout << std::endl;
  }
}
#endif

} // namespace

int main(int argc, char** argv)
{
  if (hasHelpArgument(argc, argv))
  {
    printUsage(argv[0]);
    return 0;
  }

  MpiRuntime mpi(&argc, &argv);

#if defined(CODEX_USE_MPI)
  std::cout << "Rank " << mpi.rank << " MPI initialized (size=" << mpi.size << ")" << std::endl;
#endif

  ValidationConfig config;
  config.referenceDir = "/Users/dpn/proj/bout++/ben_zhu_poincare/zperiod_5";
  config.outputDir = "./outputs";
  config.aparSinglePath = "/Users/dpn/proj/bout++/ben_zhu_poincare/apar.single.nc";
  config.aparCircPath = "/Users/dpn/proj/bout++/ben_zhu_poincare/apar.circ.nc";
  config.traceOptions.direction = 1;
  config.traceOptions.npMax = 100;
  config.tolerance = 1.0e-8;

  bool doCompare = false;
  std::string traceEngineChoice = "viskores";
  std::string viskoresDevice = "serial";
  bool viskoresDeviceSpecified = false;
  LineSpecMode lineSpecMode = LineSpecMode::none;
  std::vector<double> requestedLines;

  for (int i = 1; i < argc; ++i)
  {
    const std::string arg = argv[i];
    if (arg == "--reference-dir" && i + 1 < argc)
    {
      config.referenceDir = argv[++i];
      continue;
    }
    if (arg == "--output-dir" && i + 1 < argc)
    {
      config.outputDir = argv[++i];
      continue;
    }
    if (arg == "--single-apar" && i + 1 < argc)
    {
      config.aparSinglePath = argv[++i];
      continue;
    }
    if (arg == "--circ-apar" && i + 1 < argc)
    {
      config.aparCircPath = argv[++i];
      continue;
    }
    if (arg == "--divertor" && i + 1 < argc)
    {
      config.divertorFilter = argv[++i];
      continue;
    }
    if (arg == "--lines" && i + 1 < argc)
    {
      if (lineSpecMode == LineSpecMode::range)
      {
        if (mpi.rank == 0)
        {
          std::cerr << "Error: cannot use both --lines and --nlines\n";
        }
        return 1;
      }
      lineSpecMode = LineSpecMode::csv;
      requestedLines = parseLinesCsvDoubles(argv[++i]);
      continue;
    }
    if (arg == "--nlines" && i + 3 < argc)
    {
      if (lineSpecMode == LineSpecMode::csv)
      {
        if (mpi.rank == 0)
        {
          std::cerr << "Error: cannot use both --lines and --nlines\n";
        }
        return 1;
      }
      lineSpecMode = LineSpecMode::range;
      const double x0 = std::stod(argv[++i]);
      const double x1 = std::stod(argv[++i]);
      const int n = std::atoi(argv[++i]);
      if (n <= 0)
      {
        if (mpi.rank == 0)
        {
          std::cerr << "Error: --nlines requires N > 0\n";
        }
        return 1;
      }
      requestedLines = buildLineRange(x0, x1, n);
      continue;
    }
    if (arg == "--direction" && i + 1 < argc)
    {
      config.traceOptions.direction = std::atoi(argv[++i]);
      continue;
    }
    if (arg == "--np-max" && i + 1 < argc)
    {
      config.traceOptions.npMax = std::atoi(argv[++i]);
      continue;
    }
    if (arg == "--max-steps" && i + 1 < argc)
    {
      config.traceOptions.maxSteps = std::atoi(argv[++i]);
      continue;
    }
    if (arg == "--tol" && i + 1 < argc)
    {
      config.tolerance = std::atof(argv[++i]);
      continue;
    }
    if (arg == "--trace-engine" && i + 1 < argc)
    {
      traceEngineChoice = toLowerAscii(argv[++i]);
      continue;
    }
    if (arg == "--viskores-device" && i + 1 < argc)
    {
      viskoresDevice = toLowerAscii(argv[++i]);
      viskoresDeviceSpecified = true;
      continue;
    }
    if (arg == "--compare")
    {
      doCompare = true;
      continue;
    }

    if (mpi.rank == 0)
    {
      std::cerr << "Unknown or incomplete argument: " << arg << "\n";
      printUsage(argv[0]);
    }
    return 1;
  }

  if (config.traceOptions.direction != 1 && config.traceOptions.direction != -1)
  {
    if (mpi.rank == 0)
    {
      std::cerr << "Error: --direction must be 1 or -1\n";
    }
    return 1;
  }

  if (config.traceOptions.npMax <= 0)
  {
    if (mpi.rank == 0)
    {
      std::cerr << "Error: --np-max must be positive\n";
    }
    return 1;
  }

  if (config.traceOptions.maxSteps < 0)
  {
    if (mpi.rank == 0)
    {
      std::cerr << "Error: --max-steps must be non-negative\n";
    }
    return 1;
  }

  if (!isSupportedViskoresDeviceChoice(viskoresDevice))
  {
    if (mpi.rank == 0)
    {
      std::cerr << "Error: --viskores-device must be one of serial|openmp|kokkos\n";
    }
    return 1;
  }

  if (!isSupportedTraceEngineChoice(traceEngineChoice))
  {
    if (mpi.rank == 0)
    {
      std::cerr << "Error: --trace-engine must be one of cpu|viskores\n";
    }
    return 1;
  }

  const TraceEngine traceEngine = parseTraceEngineChoice(traceEngineChoice);

#if defined(CODEX_USE_VISKORES)
  if (traceEngine == TraceEngine::Viskores)
  {
    std::string viskoresDeviceError;
    if (!configureViskoresDevice(viskoresDevice, viskoresDeviceError))
    {
      if (mpi.rank == 0)
      {
        std::cerr << "Error: " << viskoresDeviceError << "\n";
      }
      return 1;
    }
  }
#else
  if (traceEngine == TraceEngine::Viskores)
  {
    if (mpi.rank == 0)
    {
      std::cerr << "Error: --trace-engine viskores was requested, but this binary was built without Viskores support\n";
    }
    return 1;
  }
  if (viskoresDeviceSpecified)
  {
    if (mpi.rank == 0)
    {
      std::cerr << "Error: --viskores-device was specified, but this binary was built without Viskores support\n";
    }
    return 1;
  }
#endif

  try
  {
    const SteadyClock::time_point totalStart = SteadyClock::now();
    printTraceBackendSelection(mpi, traceEngine, viskoresDevice);

    if (doCompare && mpi.size > 1)
    {
      if (mpi.rank == 0)
      {
        std::cerr << "Error: --compare is currently only supported with a single MPI rank\n";
      }
      return 1;
    }

    if (doCompare)
    {
      if (traceEngine == TraceEngine::Viskores)
      {
        if (mpi.rank == 0)
        {
          std::cerr << "Error: --compare currently uses the CPU tracer only; run trace-only mode for viskores output comparison\n";
        }
        return 1;
      }

      if (lineSpecMode == LineSpecMode::range)
      {
        if (mpi.rank == 0)
        {
          std::cerr << "Error: --nlines is only supported in trace-only mode (no --compare)\n";
        }
        return 1;
      }

      if (lineSpecMode == LineSpecMode::csv)
      {
        config.linesFilter.clear();
        config.linesFilter.reserve(requestedLines.size());
        for (double line : requestedLines)
        {
          if (!isIntegerVal(line))
          {
            if (mpi.rank == 0)
            {
              std::cerr << "Error: --compare requires integer line indices; got " << line << "\n";
            }
            return 1;
          }
          config.linesFilter.push_back(static_cast<int>(std::llround(line)));
        }
      }

      const SteadyClock::time_point compareStart = SteadyClock::now();
      ValidationSuite suite;
      const std::vector<ValidationResult> results = suite.run(config);
      const SteadyClock::time_point compareEnd = SteadyClock::now();

      int passed = 0;
      for (const auto& r : results)
      {
        if (r.pass)
        {
          ++passed;
        }
      }

      if (mpi.rank == 0)
      {
        std::cout << "\nValidation summary: " << passed << "/" << results.size() << " cases passed\n";
      }
      return (passed == static_cast<int>(results.size())) ? 0 : 2;
    }

    if (requestedLines.empty())
    {
      if (mpi.rank == 0)
      {
        std::cerr << "Error: trace-only mode requires --lines or --nlines\n";
      }
      return 1;
    }

    const bool doSingle = (config.divertorFilter == "all" || config.divertorFilter == "single");
    const bool doCirc = (config.divertorFilter == "all" || config.divertorFilter == "circ");
    const bool useCombinedOutput = (lineSpecMode == LineSpecMode::range);

    if (!doSingle && !doCirc)
    {
      if (mpi.rank == 0)
      {
        std::cerr << "Error: --divertor must be all|single|circ\n";
      }
      return 1;
    }

    const SteadyClock::time_point loadStart = SteadyClock::now();
    AparData singleData;
    AparData circData;
    if (doSingle)
    {
      singleData.load(config.aparSinglePath);
    }
    if (doCirc)
    {
      circData.load(config.aparCircPath);
    }
    const SteadyClock::time_point loadEnd = SteadyClock::now();
    const double localLoadSeconds = elapsedSeconds(loadStart, loadEnd);

    const SteadyClock::time_point traceStart = SteadyClock::now();
    const std::vector<TraceTask> allTasks =
      buildAllTasks(requestedLines, doSingle, doCirc, doSingle ? &singleData : nullptr, doCirc ? &circData : nullptr);
    const size_t localBegin = partitionBegin(allTasks.size(), mpi.rank, mpi.size);
    const size_t localEnd = partitionEnd(allTasks.size(), mpi.rank, mpi.size);
    std::vector<TraceTask> localTasks(allTasks.begin() + static_cast<std::ptrdiff_t>(localBegin),
                                      allTasks.begin() + static_cast<std::ptrdiff_t>(localEnd));

    mpi.barrier();
    for (int printer = 0; printer < mpi.size; ++printer)
    {
      if (mpi.rank == printer)
      {
        std::cout << "Rank " << mpi.rank << " seed IDs:";
        if (localBegin == localEnd)
        {
          std::cout << " (none)";
        }
        else
        {
          for (size_t seedId = localBegin; seedId < localEnd; ++seedId)
          {
            std::cout << " " << seedId;
          }
        }
        std::cout << std::endl;
      }
      mpi.barrier();
    }

    const int maxStatesPerSeed = computeMaxStateCount(config.traceOptions);
    const int maxTrajPerSeed = maxStatesPerSeed;
    const int maxPuncPerSeed = std::max(1, config.traceOptions.npMax);

    const size_t localTaskCount = localTasks.size();
    std::vector<double> ilinePerSeed(localTaskCount, 0.0);
    std::vector<int> endRegionPerSeed(localTaskCount, 0);
    std::vector<double> connectionLengthPerSeed(localTaskCount, 0.0);
    std::vector<int> stateCountPerSeed(localTaskCount, 0);
    std::vector<int> trajCountPerSeed(localTaskCount, 0);
    std::vector<int> punctureCountPerSeed(localTaskCount, 0);
    std::vector<TrajectoryState> states(localTaskCount * static_cast<size_t>(maxStatesPerSeed));
    std::vector<Point3D> trajectories(localTaskCount * static_cast<size_t>(maxTrajPerSeed));
    std::vector<PuncturePoint> punctures(localTaskCount * static_cast<size_t>(maxPuncPerSeed));
    std::vector<std::uint8_t> punctureValid(localTaskCount * static_cast<size_t>(maxPuncPerSeed), static_cast<std::uint8_t>(0));

    bool localFatal = false;
    std::string localFatalMsg;

    if (doSingle && !localFatal)
    {
      if (traceEngine == TraceEngine::Viskores)
      {
#if defined(CODEX_USE_VISKORES)
        traceLocalTasksForDivertorViskores("single",
                                           singleData,
                                           config,
                                           localTasks,
                                           mpi.rank,
                                           maxStatesPerSeed,
                                           maxTrajPerSeed,
                                           maxPuncPerSeed,
                                           ilinePerSeed,
                                           endRegionPerSeed,
                                           connectionLengthPerSeed,
                                           stateCountPerSeed,
                                           trajCountPerSeed,
                                           punctureCountPerSeed,
                                           states,
                                           trajectories,
                                           punctures,
                                           punctureValid,
                                           localFatal,
                                           localFatalMsg);
#endif
      }
      else
      {
        traceLocalTasksForDivertor("single",
                                   singleData,
                                   config,
                                   localTasks,
                                   mpi.rank,
                                   maxStatesPerSeed,
                                   maxTrajPerSeed,
                                   maxPuncPerSeed,
                                   ilinePerSeed,
                                   endRegionPerSeed,
                                   connectionLengthPerSeed,
                                   stateCountPerSeed,
                                   trajCountPerSeed,
                                   punctureCountPerSeed,
                                   states,
                                   trajectories,
                                   punctures,
                                   punctureValid,
                                   localFatal,
                                   localFatalMsg);
      }
    }

    if (doCirc && !localFatal)
    {
      if (traceEngine == TraceEngine::Viskores)
      {
#if defined(CODEX_USE_VISKORES)
        traceLocalTasksForDivertorViskores("circ",
                                           circData,
                                           config,
                                           localTasks,
                                           mpi.rank,
                                           maxStatesPerSeed,
                                           maxTrajPerSeed,
                                           maxPuncPerSeed,
                                           ilinePerSeed,
                                           endRegionPerSeed,
                                           connectionLengthPerSeed,
                                           stateCountPerSeed,
                                           trajCountPerSeed,
                                           punctureCountPerSeed,
                                           states,
                                           trajectories,
                                           punctures,
                                           punctureValid,
                                           localFatal,
                                           localFatalMsg);
#endif
      }
      else
      {
        traceLocalTasksForDivertor("circ",
                                   circData,
                                   config,
                                   localTasks,
                                   mpi.rank,
                                   maxStatesPerSeed,
                                   maxTrajPerSeed,
                                   maxPuncPerSeed,
                                   ilinePerSeed,
                                   endRegionPerSeed,
                                   connectionLengthPerSeed,
                                   stateCountPerSeed,
                                   trajCountPerSeed,
                                   punctureCountPerSeed,
                                   states,
                                   trajectories,
                                   punctures,
                                   punctureValid,
                                   localFatal,
                                   localFatalMsg);
      }
    }

    const int anyFatal = mpi.allreduceMaxInt(localFatal ? 1 : 0);
    if (anyFatal != 0)
    {
      if (localFatal)
      {
        std::cerr << "Error: " << localFatalMsg << "\n";
      }
      return 1;
    }

    const SteadyClock::time_point traceEnd = SteadyClock::now();
    const double localTraceSeconds = elapsedSeconds(traceStart, traceEnd);

    const SteadyClock::time_point outputStart = SteadyClock::now();
    PoincareOutput output;
    mpi.barrier();

    for (int writer = 0; writer < mpi.size; ++writer)
    {
      if (mpi.rank == writer)
      {
        if (useCombinedOutput)
        {
          output.writeCombinedOutputsFlat(ilinePerSeed,
                                          trajCountPerSeed,
                                          punctureCountPerSeed,
                                          maxTrajPerSeed,
                                          maxPuncPerSeed,
                                          trajectories,
                                          punctures,
                                          &punctureValid,
                                          config.outputDir,
                                          writer != 0,
                                          writer == 0);
        }
        else
        {
          for (size_t i = 0; i < localTasks.size(); ++i)
          {
            output.writeLineOutputsFlat(ilinePerSeed[i],
                                        i,
                                        maxTrajPerSeed,
                                        maxPuncPerSeed,
                                        trajCountPerSeed[i],
                                        punctureCountPerSeed[i],
                                        trajectories,
                                        punctures,
                                        &punctureValid,
                                        config.outputDir,
                                        localTasks[i].divertorTag);
          }
        }
      }
      mpi.barrier();
    }

    if (mpi.rank == 0)
    {
      if (useCombinedOutput)
      {
        std::cout << "Wrote combined outputs: " << config.outputDir << "/ip_cxx.txt, " << config.outputDir << "/ip_cxx.TP.txt"
                  << " and " << config.outputDir << "/traj_cxx.txt\n";
      }
      else
      {
        std::cout << "Wrote per-line outputs in rank order\n";
      }
      std::cout << "\nTrace-only run complete.\n";
    }
    return 0;
  }
  catch (const std::exception& ex)
  {
    std::cerr << "Error: " << ex.what() << "\n";
#if defined(CODEX_USE_MPI)
    if (mpi.size > 1)
    {
      mpi.abortAll(1);
    }
#endif
    return 1;
  }
}
