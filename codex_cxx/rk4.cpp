#include "rk4.h"

void rk4_step(const FieldData& field,
              float x_start,
              int y_start,
              float z_start,
              int region,
              int direction,
              float& x_end,
              float& z_end) {
    const float h = 1.0f;
    const float hh = 0.5f * h;
    const float h6 = h / 6.0f;

    float dxdy1 = 0.0f;
    float dzdy1 = 0.0f;
    float dxdy2 = 0.0f;
    float dzdy2 = 0.0f;
    float dxdy3 = 0.0f;
    float dzdy3 = 0.0f;
    float dxdy4 = 0.0f;
    float dzdy4 = 0.0f;

    field.evaluate_stage(x_start, z_start, y_start, region, direction, 0, dxdy1, dzdy1);
    float x1 = x_start + direction * hh * dxdy1;
    float z1 = z_start + direction * hh * dzdy1;

    field.evaluate_stage(x1, z1, y_start, region, direction, 1, dxdy2, dzdy2);
    float x2 = x_start + direction * hh * dxdy2;
    float z2 = z_start + direction * hh * dzdy2;

    field.evaluate_stage(x2, z2, y_start, region, direction, 1, dxdy3, dzdy3);
    float x3 = x_start + direction * dxdy3;
    float z3 = z_start + direction * dzdy3;

    field.evaluate_stage(x3, z3, y_start, region, direction, 2, dxdy4, dzdy4);

    x_end = x_start + direction * h6 * (dxdy1 + 2.0f * dxdy2 + 2.0f * dxdy3 + dxdy4);
    z_end = z_start + direction * h6 * (dzdy1 + 2.0f * dzdy2 + 2.0f * dzdy3 + dzdy4);
}
