#ifndef CODEX_CXX2_FIELDLINEINTEGRATOR_H
#define CODEX_CXX2_FIELDLINEINTEGRATOR_H

#include "AparFieldModel.h"
#include "Types.h"

class FieldLineIntegrator
{
public:
  explicit FieldLineIntegrator(const AparFieldModel& model);

  LineTraceResult traceLine(double iline, const TraceOptions& options) const;

private:
  const AparFieldModel& model_;

  void rk4Step(const XZPoint& start,
               int yStart,
               int region,
               int direction,
               XZPoint& end) const;
};

#endif
