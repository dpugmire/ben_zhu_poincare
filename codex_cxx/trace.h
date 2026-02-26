#ifndef TRACE_H
#define TRACE_H

#include <vector>

#include "field_data.h"

struct Point2D {
    float x;
    float y;
};

struct Point3D {
    float x;
    float y;
    float z;
};

struct TraceOptions {
    int direction;
    int nlines;
    int nturns;
    std::vector<float> lines;
};

struct InitialPoint {
    float iline;
    float xind;
    float x_start;
    int y_start;
    float z_start;
};

struct LineTraceResult {
    float iline;
    int end_region;
    float connection_length;
    std::vector<Point3D> trajectory_xyz;
    std::vector<int> puncture_steps;
    std::vector<Point3D> puncture_xyz;
    std::vector<Point2D> puncture_theta_psi;
};

std::vector<InitialPoint> build_initial_points(const FieldData& field, const TraceOptions& options);
std::vector<LineTraceResult> trace_initial_points(const FieldData& field,
                                                  const TraceOptions& options,
                                                  const std::vector<InitialPoint>& initial_points);

void trace_field_lines(const FieldData& field, const TraceOptions& options);

#endif
