#ifndef CODEX_CXX2_APARFIELDMODEL_H
#define CODEX_CXX2_APARFIELDMODEL_H

#include <vector>

#include "AparData.h"
#include "Types.h"

class AparFieldModel {
public:
    explicit AparFieldModel(const AparData& data);

    class ExecutionAccessor {
    public:
        explicit ExecutionAccessor(const AparFieldModel& model) : model_(model) {}

        const AparData& data() const { return model_.data(); }

        void evaluateStage(const XZPoint& point,
                           int yStart1b,
                           int region,
                           int direction,
                           int stage,
                           XZDeriv& deriv) const {
            model_.evaluateStage(point, yStart1b, region, direction, stage, deriv);
        }

        double interp1(const std::vector<double>& xp,
                       const std::vector<double>& fp,
                       double x) const {
            return model_.interp1(xp, fp, x);
        }

        Point3D reconstructTrajectoryXYZ(const TrajectoryState& state) const {
            return model_.reconstructTrajectoryXYZ(state);
        }

        Point3D reconstructPunctureXYZ(const Point2D& ind,
                                       double zvalue) const {
            return model_.reconstructPunctureXYZ(ind, zvalue);
        }

        double thetaFromY(double yind) const { return model_.thetaFromY(yind); }
        double psiFromX(double xind) const { return model_.psiFromX(xind); }

    private:
        const AparFieldModel& model_;
    };

    ExecutionAccessor prepareExecution() const { return ExecutionAccessor(*this); }

    const AparData& data() const { return data_; }

    void evaluateStage(const XZPoint& point,
                       int yStart1b,
                       int region,
                       int direction,
                       int stage,
                       XZDeriv& deriv) const;

    double interp1(const std::vector<double>& xp,
                   const std::vector<double>& fp,
                   double x) const;

    Point3D reconstructTrajectoryXYZ(const TrajectoryState& state) const;

    Point3D reconstructPunctureXYZ(const Point2D& ind,
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
