#ifndef CODEX_CXX2_APARDATA_H
#define CODEX_CXX2_APARDATA_H

#include <string>
#include <vector>

class AparData {
public:
    void load(const std::string& aparPath);

    int nx = 0;
    int ny = 0;
    int nz = 0;
    int zperiod = 1;
    int nzG = 0;

    int nx_cfr = 0;
    int ny_cfr = 0;

    int ixsep1 = 0;
    int ixsep2 = 0;
    int ixsep = 0;
    int nypf1 = 0;
    int nypf2 = 0;

    int divertor = 0;  // 0=circ,1=single,2=double

    int jyomp = 0;  // 0-based

    double zmin = 0.0;
    double zmax = 0.0;
    double dz_torus = 0.0;
    double xMin = 0.0;
    double xMax = 0.0;

    // 1D coordinates in MATLAB 1-based index space.
    std::vector<double> xiarray;
    std::vector<double> yiarray;
    std::vector<double> ziarray;
    std::vector<double> zarray;
    std::vector<double> xarray;
    std::vector<double> yiarray_cfr;
    std::vector<double> xiarray_cfr;

    // Raw arrays: internal layout is row-major [ix,iy] or [ix,iy,iz].
    std::vector<double> psixy;      // nx * ny
    std::vector<double> dxdy;       // nx * ny * nzG
    std::vector<double> dzdy;       // nx * ny * nzG
    std::vector<double> dxdy_p1;    // nx * nzG
    std::vector<double> dzdy_p1;    // nx * nzG
    std::vector<double> dxdy_m1;    // nx * nzG
    std::vector<double> dzdy_m1;    // nx * nzG
    std::vector<double> shiftAngle; // nx
    std::vector<double> zShift;     // nx * ny
    std::vector<double> rxy;        // nx * ny
    std::vector<double> zxy;        // nx * ny
    std::vector<double> rxy_cfr;    // nx_cfr * ny_cfr
    std::vector<double> zxy_cfr;    // nx_cfr * ny_cfr
    std::vector<double> zShift_cfr; // nx_cfr * ny_cfr

    std::vector<double> theta;      // ny
    std::vector<double> theta_cfr;  // ny_cfr

    int idx2(int ix, int iy) const { return ix * ny + iy; }
    int idx3(int ix, int iy, int iz) const { return (ix * ny + iy) * nzG + iz; }
    int idx_xz(int ix, int iz) const { return ix * nzG + iz; }
    int idx_cfr(int ix, int iy) const { return ix * ny_cfr + iy; }

    double wrapZ(double z) const;

private:
    void computeDerivedGeometry();
    void computeThetaProfiles();
};

#endif
