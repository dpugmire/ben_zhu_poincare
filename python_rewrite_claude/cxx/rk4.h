#ifndef RK4_H
#define RK4_H

#include "field_data.h"

void rk4_step(const FieldData& field,
              float x_start,
              int y_start,
              float z_start,
              int region,
              int direction,
              float& x_end,
              float& z_end);

#endif
