#ifndef FIELD_DATA_H
#define FIELD_DATA_H

#include <string>
#include <vector>

class FieldData {
public:
    FieldData();

    // Dimensions and configuration.
    int nx;
    int ny;
    int nz;
    int nzG;
    int zperiod;
    int divertor;
    int ixsep;   // 1-based
    int nypf1;   // 1-based
    int nypf2;   // 1-based
    int jyomp;   // 0-based

    // Scalars.
    float dy0;
    float dz;
    float zmin;
    float zmax;
    float dz_torus;
    float xMin;
    float xMax;
    std::string apar_gridfile_attr;

    // 2D arrays (nx * ny), flat index = ix * ny + iy.
    std::vector<float> rxy;
    std::vector<float> zxy;
    std::vector<float> psixy;
    std::vector<float> zShift;
    std::vector<float> hthe;
    std::vector<float> bxy;
    std::vector<float> btxy;
    std::vector<float> bpxy;
    std::vector<float> sinty;
    std::vector<float> bxcvx;
    std::vector<float> bxcvy;
    std::vector<float> bxcvz;
    std::vector<float> jpar0;
    std::vector<float> nu;

    // 1D coordinate arrays.
    std::vector<float> sa;       // ShiftAngle (nx)
    std::vector<float> xarray;   // psixy[:, jyomp] (nx)
    std::vector<float> xiarray;  // 1..nx
    std::vector<float> yiarray;  // 1..ny
    std::vector<float> ziarray;  // 1..nzG+1
    std::vector<float> zarray;   // 0..2pi

    // 3D arrays (nx * ny * nzG), flat index = (ix * ny + iy) * nzG + iz.
    std::vector<float> apar;
    std::vector<float> dapardx;
    std::vector<float> dapardy;
    std::vector<float> dapardz;
    std::vector<float> dxdy;
    std::vector<float> dzdy;

    // Twist-shift arrays (nx * nzG), flat index = ix * nzG + iz.
    std::vector<float> dxdy_p1;
    std::vector<float> dzdy_p1;
    std::vector<float> dxdy_m1;
    std::vector<float> dzdy_m1;

    // Theta + CFR arrays.
    std::vector<float> theta;       // ny
    std::vector<float> theta_cfr;   // ny_cfr
    std::vector<float> xiarray_cfr; // ixsep
    std::vector<float> yiarray_cfr; // ny_cfr
    std::vector<float> rxy_cfr;     // ixsep * ny_cfr
    std::vector<float> zxy_cfr;     // ixsep * ny_cfr
    std::vector<float> zs_cfr;      // ixsep * ny_cfr
    int ny_cfr;

    // I/O + setup.
    void load_grid(const std::string& grid_file);
    void load_apar(const std::string& apar_file);
    void compute_fields();
    void compute_twist_shift();
    void setup_theta_and_cfr();
    void initialize(const std::string& grid_file, const std::string& apar_file);

    // Index helpers.
    inline int idx2(int ix, int iy) const { return ix * ny + iy; }
    inline int idx3(int ix, int iy, int iz) const { return (ix * ny + iy) * nzG + iz; }
    inline int idx_xz(int ix, int iz) const { return ix * nzG + iz; }
    inline int idx_cfr(int ix, int iy) const { return ix * ny_cfr + iy; }

    // Interpolation helpers.
    float wrap_z(float z) const;
    float interp1(const std::vector<float>& xp,
                  const std::vector<float>& fp,
                  float x) const;
    float interp1_stride(const std::vector<float>& xp,
                         const float* fp,
                         int n,
                         float x,
                         int stride) const;
    float interp2_rect(const std::vector<float>& xcoords,
                       const std::vector<float>& ycoords,
                       const std::vector<float>& data,
                       int nx_local,
                       int ny_local,
                       float x,
                       float y) const;

    // x-z interpolation for field arrays.
    float interp_xz_3d_y(const std::vector<float>& data3d,
                         int y0,
                         float x,
                         float z) const;
    float interp_xz_2d(const std::vector<float>& data2d,
                       float x,
                       float z) const;
    float interp_periodic_row_3d(const std::vector<float>& data3d,
                                 int ix,
                                 int iy,
                                 float z) const;

    // Evaluate dx/dy and dz/dy for RK4 stage.
    // stage = 0 (start), 1 (half), 2 (end)
    void evaluate_stage(float x,
                        float z,
                        int y_start_1b,
                        int region,
                        int direction,
                        int stage,
                        float& dxdy_val,
                        float& dzdy_val) const;

private:
    void load_var_2d(int ncid, const char* name, std::vector<float>& out) const;
    void load_var_1d(int ncid, const char* name, int len, std::vector<float>& out) const;
    void load_var_3d_nz(int ncid, const char* name, int nz_file, std::vector<float>& out) const;
    int lower_bracket(const std::vector<float>& xp, float x) const;
};

#endif
