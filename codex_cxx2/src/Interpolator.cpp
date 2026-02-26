#include "Interpolator.h"

#include <algorithm>
#include <cmath>

namespace {

constexpr double kTwoPi = 6.283185307179586476925286766559;

void addRootIfValid(std::vector<double>& roots, double r) {
    constexpr double eps = 1.0e-6;
    if (r < -eps || r > 1.0 + eps) {
        return;
    }
    if (r < 0.0) {
        r = 0.0;
    }
    if (r > 1.0) {
        r = 1.0;
    }
    for (double existing : roots) {
        if (std::fabs(existing - r) < 1.0e-5) {
            return;
        }
    }
    roots.push_back(r);
}

}  // namespace

void NaturalCubicSpline::build(const std::vector<double>& x, const std::vector<double>& y) {
    x_ = x;
    y_ = y;
    const int n = static_cast<int>(x_.size());

    b_.assign(n, 0.0);
    c_.assign(n, 0.0);
    d_.assign(n, 0.0);

    if (n <= 2 || static_cast<int>(y_.size()) != n) {
        if (n == 2) {
            const double dx = x_[1] - x_[0];
            if (std::fabs(dx) > 1.0e-20) {
                b_[0] = (y_[1] - y_[0]) / dx;
            }
        }
        return;
    }

    const int nm1 = n - 1;
    std::vector<double> h(nm1, 0.0);
    for (int i = 0; i < nm1; ++i) {
        h[i] = x_[i + 1] - x_[i];
    }

    std::vector<double> alpha(n, 0.0);
    for (int i = 1; i < nm1; ++i) {
        if (std::fabs(h[i]) < 1.0e-20 || std::fabs(h[i - 1]) < 1.0e-20) {
            continue;
        }
        alpha[i] = 3.0 / h[i] * (y_[i + 1] - y_[i]) -
                   3.0 / h[i - 1] * (y_[i] - y_[i - 1]);
    }

    std::vector<double> l(n, 0.0), mu(n, 0.0), z(n, 0.0);
    l[0] = 1.0;

    for (int i = 1; i < nm1; ++i) {
        l[i] = 2.0 * (x_[i + 1] - x_[i - 1]) - h[i - 1] * mu[i - 1];
        if (std::fabs(l[i]) < 1.0e-20) {
            l[i] = 1.0e-20;
        }
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[nm1] = 1.0;
    c_[nm1] = 0.0;

    for (int j = nm1 - 1; j >= 0; --j) {
        c_[j] = z[j] - mu[j] * c_[j + 1];
        if (std::fabs(h[j]) < 1.0e-20) {
            b_[j] = 0.0;
            d_[j] = 0.0;
        } else {
            b_[j] = (y_[j + 1] - y_[j]) / h[j] - h[j] * (c_[j + 1] + 2.0 * c_[j]) / 3.0;
            d_[j] = (c_[j + 1] - c_[j]) / (3.0 * h[j]);
        }
    }
}

double NaturalCubicSpline::eval(double x) const {
    const int n = static_cast<int>(x_.size());
    if (n <= 0) {
        return 0.0;
    }
    if (n == 1) {
        return y_[0];
    }
    if (x <= x_.front()) {
        return y_.front();
    }
    if (x >= x_.back()) {
        return y_.back();
    }

    int lo = 0;
    int hi = n - 1;
    while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;
        if (x_[mid] <= x) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    const double t = x - x_[lo];
    return evalSegment(lo, t);
}

double NaturalCubicSpline::evalSegment(int segment, double t) const {
    return y_[segment] + b_[segment] * t + c_[segment] * t * t + d_[segment] * t * t * t;
}

int Interpolator::clampInt(int value, int lo, int hi) {
    if (value < lo) {
        return lo;
    }
    if (value > hi) {
        return hi;
    }
    return value;
}

int Interpolator::lowerBracket(const std::vector<double>& xp, double x) {
    const int n = static_cast<int>(xp.size());
    if (n < 2) {
        return 0;
    }
    if (x <= xp.front()) {
        return 0;
    }
    if (x >= xp.back()) {
        return n - 2;
    }

    int lo = 0;
    int hi = n - 1;
    while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;
        if (xp[mid] <= x) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return lo;
}

double Interpolator::linear1D(const std::vector<double>& xp,
                              const std::vector<double>& fp,
                              double x) {
    const int n = static_cast<int>(xp.size());
    if (n == 0 || static_cast<int>(fp.size()) != n) {
        return 0.0;
    }
    if (n == 1) {
        return fp[0];
    }
    if (x <= xp.front()) {
        return fp.front();
    }
    if (x >= xp.back()) {
        return fp.back();
    }

    const int i0 = lowerBracket(xp, x);
    const int i1 = i0 + 1;
    const double denom = xp[i1] - xp[i0];
    if (std::fabs(denom) < 1.0e-20) {
        return fp[i0];
    }
    const double t = (x - xp[i0]) / denom;
    return fp[i0] + t * (fp[i1] - fp[i0]);
}

double Interpolator::linear1DStride(const std::vector<double>& xp,
                                    const double* fp,
                                    int n,
                                    int stride,
                                    double x) {
    if (n <= 0 || fp == nullptr) {
        return 0.0;
    }
    if (n == 1) {
        return fp[0];
    }
    if (x <= xp.front()) {
        return fp[0];
    }
    if (x >= xp.back()) {
        return fp[(n - 1) * stride];
    }

    int lo = 0;
    int hi = n - 1;
    while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;
        if (xp[mid] <= x) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    const double denom = xp[hi] - xp[lo];
    if (std::fabs(denom) < 1.0e-20) {
        return fp[lo * stride];
    }
    const double t = (x - xp[lo]) / denom;
    const double v0 = fp[lo * stride];
    const double v1 = fp[hi * stride];
    return v0 + t * (v1 - v0);
}

double Interpolator::bilinear(const std::vector<double>& xcoords,
                              const std::vector<double>& ycoords,
                              const std::vector<double>& data,
                              int nx,
                              int ny,
                              double x,
                              double y) {
    if (nx < 2 || ny < 2) {
        return 0.0;
    }
    if (static_cast<int>(xcoords.size()) != nx ||
        static_cast<int>(ycoords.size()) != ny ||
        static_cast<int>(data.size()) != nx * ny) {
        return 0.0;
    }

    int ix0 = 0;
    if (x <= xcoords[0]) {
        ix0 = 0;
    } else if (x >= xcoords[nx - 1]) {
        ix0 = nx - 2;
    } else {
        ix0 = lowerBracket(xcoords, x);
    }
    const int ix1 = ix0 + 1;

    int iy0 = 0;
    if (y <= ycoords[0]) {
        iy0 = 0;
    } else if (y >= ycoords[ny - 1]) {
        iy0 = ny - 2;
    } else {
        iy0 = lowerBracket(ycoords, y);
    }
    const int iy1 = iy0 + 1;

    const double dx = xcoords[ix1] - xcoords[ix0];
    const double dy = ycoords[iy1] - ycoords[iy0];
    double tx = 0.0;
    double ty = 0.0;
    if (std::fabs(dx) > 1.0e-20) {
        tx = (x - xcoords[ix0]) / dx;
    }
    if (std::fabs(dy) > 1.0e-20) {
        ty = (y - ycoords[iy0]) / dy;
    }

    tx = std::max(0.0, std::min(1.0, tx));
    ty = std::max(0.0, std::min(1.0, ty));

    const double v00 = data[ix0 * ny + iy0];
    const double v10 = data[ix1 * ny + iy0];
    const double v01 = data[ix0 * ny + iy1];
    const double v11 = data[ix1 * ny + iy1];

    const double v0 = v00 + tx * (v10 - v00);
    const double v1 = v01 + tx * (v11 - v01);
    return v0 + ty * (v1 - v0);
}

double Interpolator::spline2D(const std::vector<double>& xcoords,
                              const std::vector<double>& ycoords,
                              const std::vector<double>& data,
                              int nx,
                              int ny,
                              double x,
                              double y) {
    if (nx < 2 || ny < 2) {
        return 0.0;
    }
    if (static_cast<int>(xcoords.size()) != nx ||
        static_cast<int>(ycoords.size()) != ny ||
        static_cast<int>(data.size()) != nx * ny) {
        return 0.0;
    }

    std::vector<double> alongX(nx, 0.0);
    std::vector<double> row(ny, 0.0);
    NaturalCubicSpline spline;

    for (int ix = 0; ix < nx; ++ix) {
        const int base = ix * ny;
        for (int iy = 0; iy < ny; ++iy) {
            row[iy] = data[base + iy];
        }
        spline.build(ycoords, row);
        alongX[ix] = spline.eval(y);
    }

    spline.build(xcoords, alongX);
    return spline.eval(x);
}

double Interpolator::trilinear(const std::vector<double>& xcoords,
                               const std::vector<double>& ycoords,
                               const std::vector<double>& zcoords,
                               const std::vector<double>& data,
                               int nx,
                               int ny,
                               int nz,
                               double x,
                               double y,
                               double z) {
    if (nx < 2 || ny < 2 || nz < 2) {
        return 0.0;
    }
    if (static_cast<int>(xcoords.size()) != nx ||
        static_cast<int>(ycoords.size()) != ny ||
        static_cast<int>(zcoords.size()) != nz ||
        static_cast<int>(data.size()) != nx * ny * nz) {
        return 0.0;
    }

    int ix0 = (x <= xcoords.front()) ? 0 : ((x >= xcoords.back()) ? nx - 2 : lowerBracket(xcoords, x));
    int iy0 = (y <= ycoords.front()) ? 0 : ((y >= ycoords.back()) ? ny - 2 : lowerBracket(ycoords, y));
    int iz0 = (z <= zcoords.front()) ? 0 : ((z >= zcoords.back()) ? nz - 2 : lowerBracket(zcoords, z));

    const int ix1 = ix0 + 1;
    const int iy1 = iy0 + 1;
    const int iz1 = iz0 + 1;

    const double dx = xcoords[ix1] - xcoords[ix0];
    const double dy = ycoords[iy1] - ycoords[iy0];
    const double dz = zcoords[iz1] - zcoords[iz0];

    double tx = (std::fabs(dx) > 1.0e-20) ? (x - xcoords[ix0]) / dx : 0.0;
    double ty = (std::fabs(dy) > 1.0e-20) ? (y - ycoords[iy0]) / dy : 0.0;
    double tz = (std::fabs(dz) > 1.0e-20) ? (z - zcoords[iz0]) / dz : 0.0;

    tx = std::max(0.0, std::min(1.0, tx));
    ty = std::max(0.0, std::min(1.0, ty));
    tz = std::max(0.0, std::min(1.0, tz));

    const auto at = [ny, nz, &data](int ix, int iy, int iz) {
        return data[(ix * ny + iy) * nz + iz];
    };

    const double c000 = at(ix0, iy0, iz0);
    const double c001 = at(ix0, iy0, iz1);
    const double c010 = at(ix0, iy1, iz0);
    const double c011 = at(ix0, iy1, iz1);
    const double c100 = at(ix1, iy0, iz0);
    const double c101 = at(ix1, iy0, iz1);
    const double c110 = at(ix1, iy1, iz0);
    const double c111 = at(ix1, iy1, iz1);

    const double c00 = c000 + tx * (c100 - c000);
    const double c01 = c001 + tx * (c101 - c001);
    const double c10 = c010 + tx * (c110 - c010);
    const double c11 = c011 + tx * (c111 - c011);

    const double c0 = c00 + ty * (c10 - c00);
    const double c1 = c01 + ty * (c11 - c01);

    return c0 + tz * (c1 - c0);
}

double Interpolator::lagrange4(double x,
                               double x0,
                               double x1,
                               double x2,
                               double x3,
                               double y0,
                               double y1,
                               double y2,
                               double y3) {
    const double xs[4] = {x0, x1, x2, x3};
    const double ys[4] = {y0, y1, y2, y3};

    double out = 0.0;
    for (int i = 0; i < 4; ++i) {
        double li = 1.0;
        for (int j = 0; j < 4; ++j) {
            if (i == j) {
                continue;
            }
            const double den = xs[i] - xs[j];
            if (std::fabs(den) < 1.0e-20) {
                return ys[1];
            }
            li *= (x - xs[j]) / den;
        }
        out += ys[i] * li;
    }
    return out;
}

void Interpolator::cubicStencilCentered(int baseIdx, int nNodes, int outIdx[4]) {
    if (nNodes < 4) {
        outIdx[0] = 0;
        outIdx[1] = clampInt(baseIdx, 0, std::max(0, nNodes - 1));
        outIdx[2] = outIdx[1];
        outIdx[3] = outIdx[1];
        return;
    }

    if (baseIdx <= 1) {
        outIdx[0] = 0;
        outIdx[1] = 1;
        outIdx[2] = 2;
        outIdx[3] = 3;
        return;
    }

    if (baseIdx >= nNodes - 3) {
        outIdx[0] = nNodes - 4;
        outIdx[1] = nNodes - 3;
        outIdx[2] = nNodes - 2;
        outIdx[3] = nNodes - 1;
        return;
    }

    outIdx[0] = baseIdx - 1;
    outIdx[1] = baseIdx;
    outIdx[2] = baseIdx + 1;
    outIdx[3] = baseIdx + 2;
}

std::vector<double> Interpolator::solveCubicSegmentRoots(double a0,
                                                         double a1,
                                                         double a2,
                                                         double a3) {
    std::vector<double> roots;
    constexpr double eps = 1.0e-12;

    if (std::fabs(a3) < eps) {
        if (std::fabs(a2) < eps) {
            if (std::fabs(a1) < eps) {
                return roots;
            }
            addRootIfValid(roots, -a0 / a1);
            return roots;
        }

        const double disc = a1 * a1 - 4.0 * a2 * a0;
        if (disc < 0.0) {
            return roots;
        }
        const double sq = std::sqrt(std::max(0.0, disc));
        addRootIfValid(roots, (-a1 + sq) / (2.0 * a2));
        addRootIfValid(roots, (-a1 - sq) / (2.0 * a2));
        return roots;
    }

    const double A = a2 / a3;
    const double B = a1 / a3;
    const double C = a0 / a3;

    const double p = B - A * A / 3.0;
    const double q = 2.0 * A * A * A / 27.0 - A * B / 3.0 + C;
    const double disc = 0.25 * q * q + (p * p * p) / 27.0;

    if (disc > eps) {
        const double sq = std::sqrt(disc);
        const double u = std::cbrt(-0.5 * q + sq);
        const double v = std::cbrt(-0.5 * q - sq);
        const double y = u + v;
        addRootIfValid(roots, y - A / 3.0);
    } else if (std::fabs(disc) <= eps) {
        const double u = std::cbrt(-0.5 * q);
        addRootIfValid(roots, 2.0 * u - A / 3.0);
        addRootIfValid(roots, -u - A / 3.0);
    } else {
        const double phi = std::acos(std::max(-1.0, std::min(1.0, (-0.5 * q) / std::sqrt(-(p * p * p) / 27.0))));
        const double m = 2.0 * std::sqrt(-p / 3.0);
        addRootIfValid(roots, m * std::cos(phi / 3.0) - A / 3.0);
        addRootIfValid(roots, m * std::cos((phi + kTwoPi) / 3.0) - A / 3.0);
        addRootIfValid(roots, m * std::cos((phi + 2.0 * kTwoPi) / 3.0) - A / 3.0);
    }

    std::sort(roots.begin(), roots.end());
    return roots;
}
