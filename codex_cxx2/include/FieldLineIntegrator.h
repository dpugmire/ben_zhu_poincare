#ifndef CODEX_CXX2_FIELDLINEINTEGRATOR_H
#define CODEX_CXX2_FIELDLINEINTEGRATOR_H

#include "AparFieldModel.h"
#include "Types.h"

class FieldLineIntegrator {
public:
    explicit FieldLineIntegrator(const AparFieldModel& model);

    LineTraceResult traceLine(int iline, const TraceOptions& options) const;

private:
    const AparFieldModel& model_;

    void rk4Step(double xStart,
                 int yStart,
                 double zStart,
                 int region,
                 int direction,
                 double& xEnd,
                 double& zEnd) const;
};

#endif
