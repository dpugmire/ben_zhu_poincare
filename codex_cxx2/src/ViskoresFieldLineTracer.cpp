#include "ViskoresFieldLineTracer.h"

#if defined(CODEX_USE_VISKORES)

#include <algorithm>
#include <chrono>
#include <limits>
#include <stdexcept>

#include <viskores/cont/ArrayHandle.h>
#include <viskores/cont/Invoker.h>
#include <viskores/cont/Timer.h>

#include "AparFieldModel.h"
#include "TracePostProcessor.h"
#include "ViskoresTraceWorklet.h"

namespace
{

int computeMaxStateCountFromOptions(const TraceOptions& options)
{
  constexpr int kDefaultMaxStepsPerPuncture = 200;
  const long long derivedMaxStepsLL = static_cast<long long>(kDefaultMaxStepsPerPuncture) * static_cast<long long>(std::max(1, options.npMax));
  const int derivedMaxSteps = (derivedMaxStepsLL > static_cast<long long>(std::numeric_limits<int>::max()))
    ? std::numeric_limits<int>::max()
    : static_cast<int>(std::max(1LL, derivedMaxStepsLL));

  const int maxStepsCap = (options.maxSteps > 0) ? options.maxSteps : derivedMaxSteps;
  const int nsteps = std::max(1, maxStepsCap);
  return (nsteps >= std::numeric_limits<int>::max()) ? std::numeric_limits<int>::max() : (nsteps + 1);
}

void clearSeedPunctureValid(const TraceOutputViews& outputs, std::size_t seedIndex, int maxPuncPerSeed)
{
  if (outputs.punctureValid == nullptr || maxPuncPerSeed <= 0)
  {
    return;
  }

  const std::size_t punctureBase = seedIndex * static_cast<std::size_t>(maxPuncPerSeed);
  if (punctureBase + static_cast<std::size_t>(maxPuncPerSeed) > outputs.punctureValidSize)
  {
    throw std::runtime_error("Puncture valid-mask array is smaller than numSeeds * maxPuncPerSeed");
  }

  std::fill(outputs.punctureValid + punctureBase,
            outputs.punctureValid + punctureBase + static_cast<std::size_t>(maxPuncPerSeed),
            static_cast<std::uint8_t>(0));
}

} // namespace

ViskoresFieldLineTracer::ViskoresFieldLineTracer(const AparData& data, const TraceOptions& options)
  : data_(data)
  , options_(options)
  , maxStatesPerSeed_(computeMaxStateCountFromOptions(options))
  , maxTrajPerSeed_(maxStatesPerSeed_)
  , maxPuncPerSeed_(std::max(1, options.npMax))
{
}

void ViskoresFieldLineTracer::traceLines(const std::vector<Point3D>& seeds,
                                         const std::vector<CodeXId>& globalSeedIndices,
                                         const TraceOutputViews& outputs,
                                         std::vector<double>& ilinePerSeed,
                                         std::vector<int>& endRegionPerSeed,
                                         std::vector<double>& connectionLengthPerSeed,
                                         std::vector<int>& stateCountPerSeed,
                                         std::vector<int>& trajCountPerSeed,
                                         std::vector<int>& punctureCountPerSeed,
                                         std::vector<TraceStatus>& traceStatuses,
                                         double* deviceInvokeSeconds,
                                         double* hostPostprocessSeconds) const
{
  traceStatuses.assign(seeds.size(), TraceStatus::Ok);
  if (deviceInvokeSeconds != nullptr)
  {
    *deviceInvokeSeconds = 0.0;
  }
  if (hostPostprocessSeconds != nullptr)
  {
    *hostPostprocessSeconds = 0.0;
  }

  if (seeds.empty())
  {
    return;
  }

  if (globalSeedIndices.size() != seeds.size())
  {
    throw std::runtime_error("Viskores trace batch requires one global seed index per seed");
  }
  if (options_.direction != 1 && options_.direction != -1)
  {
    throw std::runtime_error("Viskores trace requires direction to be +1 or -1");
  }
  if (outputs.states == nullptr || outputs.trajectories == nullptr || outputs.punctures == nullptr)
  {
    throw std::runtime_error("Viskores trace requires states, trajectories, and puncture output arrays");
  }

  const std::size_t totalSeeds = ilinePerSeed.size();
  if (endRegionPerSeed.size() != totalSeeds || connectionLengthPerSeed.size() != totalSeeds || stateCountPerSeed.size() != totalSeeds ||
      trajCountPerSeed.size() != totalSeeds || punctureCountPerSeed.size() != totalSeeds)
  {
    throw std::runtime_error("Per-seed metadata arrays must all have the same size");
  }
  if (outputs.statesSize < totalSeeds * static_cast<std::size_t>(maxStatesPerSeed_) ||
      outputs.trajectoriesSize < totalSeeds * static_cast<std::size_t>(maxTrajPerSeed_) ||
      outputs.puncturesSize < totalSeeds * static_cast<std::size_t>(maxPuncPerSeed_))
  {
    throw std::runtime_error("Output arrays are smaller than numSeeds * maxXPerSeed");
  }
  if (outputs.punctureValid != nullptr && outputs.punctureValidSize < totalSeeds * static_cast<std::size_t>(maxPuncPerSeed_))
  {
    throw std::runtime_error("Puncture valid-mask array is smaller than numSeeds * maxPuncPerSeed");
  }

  ViskoresAparField field(data_);
  ViskoresTraceStatesWorklet worklet(maxStatesPerSeed_, options_.direction);
  viskores::cont::Invoker invoker;

  auto seedHandle = viskores::cont::make_ArrayHandle(seeds, viskores::CopyFlag::On);

  viskores::cont::ArrayHandle<viskores::Id> batchStateCounts;
  viskores::cont::ArrayHandle<viskores::Id> batchEndRegions;
  viskores::cont::ArrayHandle<viskores::Id> batchStatusCodes;
  viskores::cont::ArrayHandle<TrajectoryState> batchStates;

  batchStateCounts.Allocate(seedHandle.GetNumberOfValues());
  batchEndRegions.Allocate(seedHandle.GetNumberOfValues());
  batchStatusCodes.Allocate(seedHandle.GetNumberOfValues());
  batchStates.Allocate(seedHandle.GetNumberOfValues() * static_cast<viskores::Id>(maxStatesPerSeed_));

  using SteadyClock = std::chrono::steady_clock;

  viskores::cont::Timer deviceTimer(invoker.GetDevice());
  deviceTimer.Start();
  invoker(worklet, seedHandle, field, batchStateCounts, batchEndRegions, batchStatusCodes, batchStates);
  deviceTimer.Stop();
  if (deviceInvokeSeconds != nullptr)
  {
    *deviceInvokeSeconds = static_cast<double>(deviceTimer.GetElapsedTime());
  }
  std::cout << "Device invoke seconds: " << *deviceInvokeSeconds << std::endl;

  const SteadyClock::time_point hostPostprocessStart = SteadyClock::now();
  const auto stateCountPortal = batchStateCounts.ReadPortal();
  const auto endRegionPortal = batchEndRegions.ReadPortal();
  const auto statusPortal = batchStatusCodes.ReadPortal();
  const auto statePortal = batchStates.ReadPortal();

  AparFieldModel model(data_);

  for (std::size_t batchIdx = 0; batchIdx < seeds.size(); ++batchIdx)
  {
    const std::size_t globalSeedIndex = static_cast<std::size_t>(globalSeedIndices[batchIdx]);
    if (globalSeedIndex >= totalSeeds)
    {
      throw std::runtime_error("Global seed index is outside the output metadata arrays");
    }

    ilinePerSeed[globalSeedIndex] = seeds[batchIdx].x;
    endRegionPerSeed[globalSeedIndex] = static_cast<int>(endRegionPortal.Get(static_cast<viskores::Id>(batchIdx)));
    stateCountPerSeed[globalSeedIndex] =
      std::max(0, std::min(static_cast<int>(stateCountPortal.Get(static_cast<viskores::Id>(batchIdx))), maxStatesPerSeed_));
    trajCountPerSeed[globalSeedIndex] = 0;
    punctureCountPerSeed[globalSeedIndex] = 0;
    connectionLengthPerSeed[globalSeedIndex] = 0.0;
    clearSeedPunctureValid(outputs, globalSeedIndex, maxPuncPerSeed_);

    const int statusCode = static_cast<int>(statusPortal.Get(static_cast<viskores::Id>(batchIdx)));
    traceStatuses[batchIdx] = static_cast<TraceStatus>(statusCode);

    const std::size_t globalStateBase = globalSeedIndex * static_cast<std::size_t>(maxStatesPerSeed_);
    const std::size_t batchStateBase = batchIdx * static_cast<std::size_t>(maxStatesPerSeed_);
    for (int i = 0; i < stateCountPerSeed[globalSeedIndex]; ++i)
    {
      outputs.states[globalStateBase + static_cast<std::size_t>(i)] =
        statePortal.Get(static_cast<viskores::Id>(batchStateBase + static_cast<std::size_t>(i)));
    }

    if (traceStatuses[batchIdx] == TraceStatus::Ok || traceStatuses[batchIdx] == TraceStatus::MaxStepLimitReached)
    {
      TracePostProcessor::rebuildSeedOutputs(model,
                                             options_,
                                             globalSeedIndex,
                                             maxStatesPerSeed_,
                                             maxTrajPerSeed_,
                                             maxPuncPerSeed_,
                                             outputs,
                                             stateCountPerSeed[globalSeedIndex],
                                             trajCountPerSeed[globalSeedIndex],
                                             punctureCountPerSeed[globalSeedIndex],
                                             connectionLengthPerSeed[globalSeedIndex]);
    }
  }

  if (hostPostprocessSeconds != nullptr)
  {
    *hostPostprocessSeconds = std::chrono::duration<double>(SteadyClock::now() - hostPostprocessStart).count();
  }
}

#endif
