#ifndef CODEX_CXX2_INTERPOLATOR_H
#define CODEX_CXX2_INTERPOLATOR_H

#include <vector>

class NaturalCubicSpline {
public:
    void build(const std::vector<double>& x, const std::vector<double>& y);
    double eval(double x) const;
    double evalSegment(int segment, double t) const;

    const std::vector<double>& knots() const { return x_; }
    const std::vector<double>& values() const { return y_; }
    const std::vector<double>& b() const { return b_; }
    const std::vector<double>& c() const { return c_; }
    const std::vector<double>& d() const { return d_; }

private:
    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> b_;
    std::vector<double> c_;
    std::vector<double> d_;
};

class Interpolator {
public:
    static int clampInt(int value, int lo, int hi);
    static int lowerBracket(const std::vector<double>& xp, double x);

    static double linear1D(const std::vector<double>& xp,
                           const std::vector<double>& fp,
                           double x);

    static double linear1DStride(const std::vector<double>& xp,
                                 const double* fp,
                                 int n,
                                 int stride,
                                 double x);

    static double bilinear(const std::vector<double>& xcoords,
                           const std::vector<double>& ycoords,
                           const std::vector<double>& data,
                           int nx,
                           int ny,
                           double x,
                           double y);

    static double spline2D(const std::vector<double>& xcoords,
                           const std::vector<double>& ycoords,
                           const std::vector<double>& data,
                           int nx,
                           int ny,
                           double x,
                           double y);

    static double trilinear(const std::vector<double>& xcoords,
                            const std::vector<double>& ycoords,
                            const std::vector<double>& zcoords,
                            const std::vector<double>& data,
                            int nx,
                            int ny,
                            int nz,
                            double x,
                            double y,
                            double z);

    static double lagrange4(double x,
                            double x0,
                            double x1,
                            double x2,
                            double x3,
                            double y0,
                            double y1,
                            double y2,
                            double y3);

    static void cubicStencilCentered(int baseIdx, int nNodes, int outIdx[4]);

    static std::vector<double> solveCubicSegmentRoots(double a0,
                                                      double a1,
                                                      double a2,
                                                      double a3);
};

#endif
