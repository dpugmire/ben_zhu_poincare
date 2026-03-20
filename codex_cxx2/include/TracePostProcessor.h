#ifndef CODEX_CXX2_TRACEPOSTPROCESSOR_H
#define CODEX_CXX2_TRACEPOSTPROCESSOR_H

#include <cstddef>

#include "AparFieldModel.h"
#include "Types.h"

namespace TracePostProcessor
{

void rebuildSeedOutputs(const AparFieldModel& model,
                        const TraceOptions& options,
                        std::size_t seedIndex,
                        int maxStatesPerSeed,
                        int maxTrajPerSeed,
                        int maxPuncPerSeed,
                        const TraceOutputViews& outputs,
                        int stateCount,
                        int& trajCount,
                        int& punctureCount,
                        double& connectionLength);

}

#endif
