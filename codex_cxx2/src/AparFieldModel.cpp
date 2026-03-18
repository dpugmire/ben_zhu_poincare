#include "AparFieldModel.h"

#include <cmath>
#include <stdexcept>

#include "Interpolator.h"

namespace {

double samplePeriodicRow3D(const AparData& data,
                           const std::vector<double>& data3d,
                           int ix,
                           int iy,
                           int izExt) {
    if (izExt <= 0) {
        return data3d[data.idx3(ix, iy, 0)];
    }
    if (izExt >= data.nzG) {
        return data3d[data.idx3(ix, iy, 0)];
    }
    return data3d[data.idx3(ix, iy, izExt)];
}

double sampleClampedRow3D(const AparData& data,
                          const std::vector<double>& data3d,
                          int ix,
                          int iy,
                          int izExt) {
    if (izExt <= 0) {
        return data3d[data.idx3(ix, iy, 0)];
    }
    if (izExt >= data.nzG) {
        return data3d[data.idx3(ix, iy, data.nzG - 1)];
    }
    return data3d[data.idx3(ix, iy, izExt)];
}

double sampleClampedXZ(const AparData& data,
                       const std::vector<double>& data2d,
                       int ix,
                       int izExt) {
    if (izExt <= 0) {
        return data2d[data.idx_xz(ix, 0)];
    }
    if (izExt >= data.nzG) {
        return data2d[data.idx_xz(ix, data.nzG - 1)];
    }
    return data2d[data.idx_xz(ix, izExt)];
}

}  // namespace

AparFieldModel::AparFieldModel(const AparData& data) : data_(data) {}

double AparFieldModel::interp1(const std::vector<double>& xp,
                               const std::vector<double>& fp,
                               double x) const {
    return Interpolator::linear1D(xp, fp, x);
}

double AparFieldModel::interpXIndex2D(const std::vector<double>& data2d,
                                      int yidx0,
                                      double xind1b) const {
    yidx0 = Interpolator::clampInt(yidx0, 0, data_.ny - 1);

    if (data_.nx < 2) {
        return data2d[data_.idx2(0, yidx0)];
    }

    double xf = xind1b - 1.0;
    if (xf <= 0.0) {
        return data2d[data_.idx2(0, yidx0)];
    }
    if (xf >= static_cast<double>(data_.nx - 1)) {
        return data2d[data_.idx2(data_.nx - 1, yidx0)];
    }

    int ix0 = static_cast<int>(std::floor(xf));
    ix0 = Interpolator::clampInt(ix0, 0, data_.nx - 2);
    int ix1 = ix0 + 1;
    double tx = xf - static_cast<double>(ix0);

    const double v0 = data2d[data_.idx2(ix0, yidx0)];
    const double v1 = data2d[data_.idx2(ix1, yidx0)];
    return v0 + tx * (v1 - v0);
}

double AparFieldModel::interpPeriodicRow3D(const std::vector<double>& data3d,
                                            int ix,
                                            int iy,
                                            double z) const {
    ix = Interpolator::clampInt(ix, 0, data_.nx - 1);
    iy = Interpolator::clampInt(iy, 0, data_.ny - 1);

    const double zWrapped = data_.wrapZ(z);
    const double zf = zWrapped / data_.dz_torus;
    int iz0 = static_cast<int>(std::floor(zf));
    iz0 = Interpolator::clampInt(iz0, 0, data_.nzG - 1);
    const int iz1 = (iz0 + 1) % data_.nzG;
    double t = zf - static_cast<double>(iz0);
    t = std::max(0.0, std::min(1.0, t));

    const double v0 = data3d[data_.idx3(ix, iy, iz0)];
    const double v1 = data3d[data_.idx3(ix, iy, iz1)];
    return v0 + t * (v1 - v0);
}

double AparFieldModel::interpPeriodicRow3DSpline(const std::vector<double>& data3d,
                                                  int ix,
                                                  int iy,
                                                  double z) const {
    ix = Interpolator::clampInt(ix, 0, data_.nx - 1);
    iy = Interpolator::clampInt(iy, 0, data_.ny - 1);
    if (data_.nzG < 4) {
        return interpPeriodicRow3D(data3d, ix, iy, z);
    }

    const double zq = data_.wrapZ(z);
    int izBase = static_cast<int>(std::floor(zq / data_.dz_torus));
    izBase = Interpolator::clampInt(izBase, 0, data_.nzG - 1);

    int izs[4] = {0, 1, 2, 3};
    Interpolator::cubicStencilCentered(izBase, data_.nzG + 1, izs);

    const double zc0 = static_cast<double>(izs[0]) * data_.dz_torus;
    const double zc1 = static_cast<double>(izs[1]) * data_.dz_torus;
    const double zc2 = static_cast<double>(izs[2]) * data_.dz_torus;
    const double zc3 = static_cast<double>(izs[3]) * data_.dz_torus;

    const double zv0 = samplePeriodicRow3D(data_, data3d, ix, iy, izs[0]);
    const double zv1 = samplePeriodicRow3D(data_, data3d, ix, iy, izs[1]);
    const double zv2 = samplePeriodicRow3D(data_, data3d, ix, iy, izs[2]);
    const double zv3 = samplePeriodicRow3D(data_, data3d, ix, iy, izs[3]);

    return Interpolator::lagrange4(zq, zc0, zc1, zc2, zc3, zv0, zv1, zv2, zv3);
}

double AparFieldModel::interpXZ3DAtY(const std::vector<double>& data3d,
                                     int y0,
                                     double x,
                                     double z) const {
    if (data_.nx < 2 || data_.nzG < 1) {
        return 0.0;
    }

    y0 = Interpolator::clampInt(y0, 0, data_.ny - 1);

    int ix0 = 0;
    double tx = 0.0;
    if (x <= data_.xarray.front()) {
        ix0 = 0;
        tx = 0.0;
    } else if (x >= data_.xarray.back()) {
        ix0 = data_.nx - 2;
        tx = 1.0;
    } else {
        ix0 = Interpolator::lowerBracket(data_.xarray, x);
        const double denom = data_.xarray[ix0 + 1] - data_.xarray[ix0];
        tx = (std::fabs(denom) > 1.0e-20) ? (x - data_.xarray[ix0]) / denom : 0.0;
    }
    tx = std::max(0.0, std::min(1.0, tx));

    const double zWrapped = data_.wrapZ(z);
    const double zf = zWrapped / data_.dz_torus;
    int iz0 = static_cast<int>(std::floor(zf));
    iz0 = Interpolator::clampInt(iz0, 0, data_.nzG - 1);
    const int iz1 = (iz0 + 1) % data_.nzG;
    double tz = zf - static_cast<double>(iz0);
    tz = std::max(0.0, std::min(1.0, tz));

    int ix1 = ix0 + 1;
    if (ix1 >= data_.nx) {
        ix1 = data_.nx - 1;
    }

    const double v00 = data3d[data_.idx3(ix0, y0, iz0)];
    const double v01 = data3d[data_.idx3(ix0, y0, iz1)];
    const double v10 = data3d[data_.idx3(ix1, y0, iz0)];
    const double v11 = data3d[data_.idx3(ix1, y0, iz1)];

    const double v0 = v00 + tz * (v01 - v00);
    const double v1 = v10 + tz * (v11 - v10);
    return v0 + tx * (v1 - v0);
}

double AparFieldModel::interpXZ3DAtYSpline(const std::vector<double>& data3d,
                                           int y0,
                                           double x,
                                           double z) const {
    if (data_.nx < 4 || data_.nzG < 4) {
        return interpXZ3DAtY(data3d, y0, x, z);
    }

    y0 = Interpolator::clampInt(y0, 0, data_.ny - 1);

    double xq = x;
    xq = std::max(data_.xarray.front(), std::min(data_.xarray.back(), xq));

    const int ixBase = Interpolator::lowerBracket(data_.xarray, xq);
    int ixs[4] = {0, 1, 2, 3};
    Interpolator::cubicStencilCentered(ixBase, data_.nx, ixs);

    const double zq = data_.wrapZ(z);
    int izBase = static_cast<int>(std::floor(zq / data_.dz_torus));
    izBase = Interpolator::clampInt(izBase, 0, data_.nzG - 1);

    int izs[4] = {0, 1, 2, 3};
    Interpolator::cubicStencilCentered(izBase, data_.nzG + 1, izs);

    const double xvals[4] = {
        data_.xarray[ixs[0]], data_.xarray[ixs[1]], data_.xarray[ixs[2]], data_.xarray[ixs[3]]};

    double zinterp[4] = {0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < 4; ++i) {
        const int ix = ixs[i];

        const double zc0 = static_cast<double>(izs[0]) * data_.dz_torus;
        const double zc1 = static_cast<double>(izs[1]) * data_.dz_torus;
        const double zc2 = static_cast<double>(izs[2]) * data_.dz_torus;
        const double zc3 = static_cast<double>(izs[3]) * data_.dz_torus;

        const double zv0 = sampleClampedRow3D(data_, data3d, ix, y0, izs[0]);
        const double zv1 = sampleClampedRow3D(data_, data3d, ix, y0, izs[1]);
        const double zv2 = sampleClampedRow3D(data_, data3d, ix, y0, izs[2]);
        const double zv3 = sampleClampedRow3D(data_, data3d, ix, y0, izs[3]);

        zinterp[i] = Interpolator::lagrange4(zq, zc0, zc1, zc2, zc3, zv0, zv1, zv2, zv3);
    }

    return Interpolator::lagrange4(
        xq,
        xvals[0], xvals[1], xvals[2], xvals[3],
        zinterp[0], zinterp[1], zinterp[2], zinterp[3]);
}

double AparFieldModel::interpXZ2D(const std::vector<double>& data2d,
                                  double x,
                                  double z) const {
    if (data_.nx < 2 || data_.nzG < 1) {
        return 0.0;
    }

    int ix0 = 0;
    double tx = 0.0;
    if (x <= data_.xarray.front()) {
        ix0 = 0;
        tx = 0.0;
    } else if (x >= data_.xarray.back()) {
        ix0 = data_.nx - 2;
        tx = 1.0;
    } else {
        ix0 = Interpolator::lowerBracket(data_.xarray, x);
        const double denom = data_.xarray[ix0 + 1] - data_.xarray[ix0];
        tx = (std::fabs(denom) > 1.0e-20) ? (x - data_.xarray[ix0]) / denom : 0.0;
    }
    tx = std::max(0.0, std::min(1.0, tx));

    const double zWrapped = data_.wrapZ(z);
    const double zf = zWrapped / data_.dz_torus;
    int iz0 = static_cast<int>(std::floor(zf));
    iz0 = Interpolator::clampInt(iz0, 0, data_.nzG - 1);
    const int iz1 = (iz0 + 1) % data_.nzG;
    double tz = zf - static_cast<double>(iz0);
    tz = std::max(0.0, std::min(1.0, tz));

    int ix1 = ix0 + 1;
    if (ix1 >= data_.nx) {
        ix1 = data_.nx - 1;
    }

    const double v00 = data2d[data_.idx_xz(ix0, iz0)];
    const double v01 = data2d[data_.idx_xz(ix0, iz1)];
    const double v10 = data2d[data_.idx_xz(ix1, iz0)];
    const double v11 = data2d[data_.idx_xz(ix1, iz1)];

    const double v0 = v00 + tz * (v01 - v00);
    const double v1 = v10 + tz * (v11 - v10);
    return v0 + tx * (v1 - v0);
}

double AparFieldModel::interpXZ2DSpline(const std::vector<double>& data2d,
                                        double x,
                                        double z) const {
    if (data_.nx < 4 || data_.nzG < 4) {
        return interpXZ2D(data2d, x, z);
    }

    double xq = std::max(data_.xarray.front(), std::min(data_.xarray.back(), x));
    const int ixBase = Interpolator::lowerBracket(data_.xarray, xq);
    int ixs[4] = {0, 1, 2, 3};
    Interpolator::cubicStencilCentered(ixBase, data_.nx, ixs);

    const double zq = data_.wrapZ(z);
    int izBase = static_cast<int>(std::floor(zq / data_.dz_torus));
    izBase = Interpolator::clampInt(izBase, 0, data_.nzG - 1);
    int izs[4] = {0, 1, 2, 3};
    Interpolator::cubicStencilCentered(izBase, data_.nzG + 1, izs);

    const double xvals[4] = {
        data_.xarray[ixs[0]], data_.xarray[ixs[1]], data_.xarray[ixs[2]], data_.xarray[ixs[3]]};

    double zinterp[4] = {0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < 4; ++i) {
        const int ix = ixs[i];

        const double zc0 = static_cast<double>(izs[0]) * data_.dz_torus;
        const double zc1 = static_cast<double>(izs[1]) * data_.dz_torus;
        const double zc2 = static_cast<double>(izs[2]) * data_.dz_torus;
        const double zc3 = static_cast<double>(izs[3]) * data_.dz_torus;

        const double zv0 = sampleClampedXZ(data_, data2d, ix, izs[0]);
        const double zv1 = sampleClampedXZ(data_, data2d, ix, izs[1]);
        const double zv2 = sampleClampedXZ(data_, data2d, ix, izs[2]);
        const double zv3 = sampleClampedXZ(data_, data2d, ix, izs[3]);

        zinterp[i] = Interpolator::lagrange4(zq, zc0, zc1, zc2, zc3, zv0, zv1, zv2, zv3);
    }

    return Interpolator::lagrange4(
        xq,
        xvals[0], xvals[1], xvals[2], xvals[3],
        zinterp[0], zinterp[1], zinterp[2], zinterp[3]);
}

void AparFieldModel::evaluateStage(const XZPoint& point,
                                   int yStart1b,
                                   int region,
                                   int direction,
                                   int stage,
                                   XZDeriv& deriv) const {
    stage = Interpolator::clampInt(stage, 0, 2);

    int yp = yStart1b - 1;
    yp = Interpolator::clampInt(yp, 0, data_.ny - 1);

    bool useTwist = false;
    bool usePlus = false;
    int yn = yp;

    if (direction == 1) {
        if (region == 0 && yStart1b == data_.nypf2) {
            useTwist = true;
            usePlus = true;
        } else {
            yn = yStart1b;
            yn = Interpolator::clampInt(yn, 0, data_.ny - 1);
        }
    } else if (direction == -1) {
        if (region == 0 && yStart1b == (data_.nypf1 + 1)) {
            useTwist = true;
            usePlus = false;
        } else {
            yn = yStart1b - 2;
            yn = Interpolator::clampInt(yn, 0, data_.ny - 1);
        }
    } else {
        throw std::runtime_error("direction must be +1 or -1");
    }

    if (stage == 0) {
        deriv.dxdy = interpXZ3DAtYSpline(data_.dxdy, yp, point.x, point.z);
        deriv.dzdy = interpXZ3DAtYSpline(data_.dzdy, yp, point.x, point.z);
        return;
    }

    if (stage == 2) {
        if (useTwist) {
            if (usePlus) {
                deriv.dxdy = interpXZ2DSpline(data_.dxdy_p1, point.x, point.z);
                deriv.dzdy = interpXZ2DSpline(data_.dzdy_p1, point.x, point.z);
            } else {
                deriv.dxdy = interpXZ2DSpline(data_.dxdy_m1, point.x, point.z);
                deriv.dzdy = interpXZ2DSpline(data_.dzdy_m1, point.x, point.z);
            }
        } else {
            deriv.dxdy = interpXZ3DAtYSpline(data_.dxdy, yn, point.x, point.z);
            deriv.dzdy = interpXZ3DAtYSpline(data_.dzdy, yn, point.x, point.z);
        }
        return;
    }

    const double dxP = interpXZ3DAtYSpline(data_.dxdy, yp, point.x, point.z);
    const double dzP = interpXZ3DAtYSpline(data_.dzdy, yp, point.x, point.z);

    double dxN = 0.0;
    double dzN = 0.0;
    if (useTwist) {
        if (usePlus) {
            dxN = interpXZ2DSpline(data_.dxdy_p1, point.x, point.z);
            dzN = interpXZ2DSpline(data_.dzdy_p1, point.x, point.z);
        } else {
            dxN = interpXZ2DSpline(data_.dxdy_m1, point.x, point.z);
            dzN = interpXZ2DSpline(data_.dzdy_m1, point.x, point.z);
        }
    } else {
        dxN = interpXZ3DAtYSpline(data_.dxdy, yn, point.x, point.z);
        dzN = interpXZ3DAtYSpline(data_.dzdy, yn, point.x, point.z);
    }

    deriv.dxdy = 0.5 * (dxP + dxN);
    deriv.dzdy = 0.5 * (dzP + dzN);
}

Point3D AparFieldModel::reconstructTrajectoryXYZ(const TrajectoryState& state) const {
    int yidx = static_cast<int>(std::round(state.ind.y)) - 1;
    yidx = Interpolator::clampInt(yidx, 0, data_.ny - 1);

    const double rxyValue = interpXIndex2D(data_.rxy, yidx, state.ind.x);
    const double zsValue = interpXIndex2D(data_.zShift, yidx, state.ind.x);
    const double zxyValue = interpXIndex2D(data_.zxy, yidx, state.ind.x);
    const double zValue = interp1(data_.ziarray, data_.zarray, state.ind.z);

    const double x3dTmp = rxyValue * std::cos(zsValue);
    const double y3dTmp = rxyValue * std::sin(zsValue);

    Point3D p;
    p.x = x3dTmp * std::cos(zValue) - y3dTmp * std::sin(zValue);
    p.y = x3dTmp * std::sin(zValue) + y3dTmp * std::cos(zValue);
    p.z = zxyValue;
    return p;
}

Point3D AparFieldModel::reconstructPunctureXYZ(const Point2D& ind,
                                                double zvalue) const {
    Point3D p;

    double rxyValue = 0.0;
    double zxyValue = 0.0;
    double zsValue = 0.0;

    if (ind.x < static_cast<double>(data_.ixsep) + 0.5) {
        rxyValue = Interpolator::spline2D(
            data_.xiarray_cfr, data_.yiarray_cfr, data_.rxy_cfr, data_.nx_cfr, data_.ny_cfr, ind.x, ind.y);
        zxyValue = Interpolator::spline2D(
            data_.xiarray_cfr, data_.yiarray_cfr, data_.zxy_cfr, data_.nx_cfr, data_.ny_cfr, ind.x, ind.y);
        zsValue = Interpolator::spline2D(
            data_.xiarray_cfr, data_.yiarray_cfr, data_.zShift_cfr, data_.nx_cfr, data_.ny_cfr, ind.x, ind.y);
    } else {
        rxyValue = Interpolator::spline2D(
            data_.xiarray, data_.yiarray, data_.rxy, data_.nx, data_.ny, ind.x, ind.y);
        zxyValue = Interpolator::spline2D(
            data_.xiarray, data_.yiarray, data_.zxy, data_.nx, data_.ny, ind.x, ind.y);
        zsValue = Interpolator::spline2D(
            data_.xiarray, data_.yiarray, data_.zShift, data_.nx, data_.ny, ind.x, ind.y);
    }

    const double ipx3dTmp = rxyValue * std::cos(zsValue);
    const double ipy3dTmp = rxyValue * std::sin(zsValue);

    p.x = ipx3dTmp * std::cos(zvalue) - ipy3dTmp * std::sin(zvalue);
    p.y = ipx3dTmp * std::sin(zvalue) + ipy3dTmp * std::cos(zvalue);
    p.z = zxyValue;

    return p;
}

double AparFieldModel::thetaFromY(double yind) const {
    return interp1(data_.yiarray_cfr, data_.theta_cfr, yind);
}

double AparFieldModel::psiFromX(double xind) const {
    return interp1(data_.xiarray, data_.xarray, xind);
}
