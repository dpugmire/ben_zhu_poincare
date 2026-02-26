#ifndef CODEX_CXX2_APARFIELDMODEL_H
#define CODEX_CXX2_APARFIELDMODEL_H

#include <vector>

#include "AparData.h"
#include "Types.h"

class AparFieldModel {
public:
    explicit AparFieldModel(const AparData& data);

    const AparData& data() const { return data_; }

    void evaluateStage(double x,
                       double z,
                       int yStart1b,
                       int region,
                       int direction,
                       int stage,
                       double& dxdy,
                       double& dzdy) const;

    double interp1(const std::vector<double>& xp,
                   const std::vector<double>& fp,
                   double x) const;

    Point3D reconstructTrajectoryXYZ(const TrajectoryState& state) const;

    Point3D reconstructPunctureXYZ(double xind,
                                   double yind,
                                   double zvalue) const;

    double thetaFromY(double yind) const;
    double psiFromX(double xind) const;

private:
    const AparData& data_;

    double interpXIndex2D(const std::vector<double>& data2d,
                          int yidx0,
                          double xind1b) const;

    double interpPeriodicRow3D(const std::vector<double>& data3d,
                               int ix,
                               int iy,
                               double z) const;

    double interpPeriodicRow3DSpline(const std::vector<double>& data3d,
                                     int ix,
                                     int iy,
                                     double z) const;

    double interpXZ3DAtY(const std::vector<double>& data3d,
                         int y0,
                         double x,
                         double z) const;

    double interpXZ3DAtYSpline(const std::vector<double>& data3d,
                               int y0,
                               double x,
                               double z) const;

    double interpXZ2D(const std::vector<double>& data2d,
                      double x,
                      double z) const;

    double interpXZ2DSpline(const std::vector<double>& data2d,
                            double x,
                            double z) const;
};

#endif
