#ifndef TRACE_H
#define TRACE_H

#include <vector>

#include "field_data.h"

struct TraceOptions {
    int direction;
    int nlines;
    int nturns;
    std::vector<int> lines;
};

void trace_field_lines(const FieldData& field, const TraceOptions& options);

#endif
