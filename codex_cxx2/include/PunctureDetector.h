#ifndef CODEX_CXX2_PUNCTUREDETECTOR_H
#define CODEX_CXX2_PUNCTUREDETECTOR_H

#include "AparFieldModel.h"
#include "Types.h"

class PunctureDetector {
public:
    explicit PunctureDetector(const AparFieldModel& model);

    void detect(LineTraceResult& line,
                int direction,
                int npMax) const;

private:
    const AparFieldModel& model_;
};

#endif
