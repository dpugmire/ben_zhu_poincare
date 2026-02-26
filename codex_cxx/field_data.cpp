#include "field_data.h"

#include <netcdf.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {

const float kPi = 3.14159265358979323846f;
const float kMu0 = 4.0f * kPi * 1.0e-7f;

inline int clamp_int_local(int v, int lo, int hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

struct NaturalCubicSpline {
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> b;
    std::vector<float> c;
    std::vector<float> d;
    int n = 0;

    void build(const std::vector<float>& xin, const std::vector<float>& yin) {
        x = xin;
        y = yin;
        n = static_cast<int>(x.size());
        b.assign(n, 0.0f);
        c.assign(n, 0.0f);
        d.assign(n, 0.0f);

        if (n <= 2 || static_cast<int>(y.size()) != n) {
            if (n == 2) {
                float dx = x[1] - x[0];
                if (std::fabs(dx) > 1.0e-20f) {
                    b[0] = (y[1] - y[0]) / dx;
                }
            }
            return;
        }

        const int nm1 = n - 1;
        std::vector<float> h(nm1, 0.0f);
        for (int i = 0; i < nm1; ++i) {
            h[i] = x[i + 1] - x[i];
        }

        std::vector<float> alpha(n, 0.0f);
        for (int i = 1; i < nm1; ++i) {
            if (std::fabs(h[i]) < 1.0e-20f || std::fabs(h[i - 1]) < 1.0e-20f) continue;
            alpha[i] = 3.0f / h[i] * (y[i + 1] - y[i]) -
                       3.0f / h[i - 1] * (y[i] - y[i - 1]);
        }

        std::vector<float> l(n, 0.0f), mu(n, 0.0f), z(n, 0.0f);
        l[0] = 1.0f;
        for (int i = 1; i < nm1; ++i) {
            l[i] = 2.0f * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
            if (std::fabs(l[i]) < 1.0e-20f) l[i] = 1.0e-20f;
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }
        l[nm1] = 1.0f;
        c[nm1] = 0.0f;

        for (int j = nm1 - 1; j >= 0; --j) {
            c[j] = z[j] - mu[j] * c[j + 1];
            if (std::fabs(h[j]) < 1.0e-20f) {
                b[j] = 0.0f;
                d[j] = 0.0f;
            } else {
                b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2.0f * c[j]) / 3.0f;
                d[j] = (c[j + 1] - c[j]) / (3.0f * h[j]);
            }
        }
    }

    float eval(float xv) const {
        if (n <= 0) return 0.0f;
        if (n == 1) return y[0];
        if (xv <= x.front()) return y.front();
        if (xv >= x.back()) return y.back();

        int lo = 0;
        int hi = n - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (x[mid] <= xv) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        float dx = xv - x[lo];
        return y[lo] + b[lo] * dx + c[lo] * dx * dx + d[lo] * dx * dx * dx;
    }
};

float lagrange4(float x,
                float x0, float x1, float x2, float x3,
                float y0, float y1, float y2, float y3) {
    const float xs[4] = {x0, x1, x2, x3};
    const float ys[4] = {y0, y1, y2, y3};
    float out = 0.0f;
    for (int i = 0; i < 4; ++i) {
        float li = 1.0f;
        for (int j = 0; j < 4; ++j) {
            if (i == j) continue;
            float den = xs[i] - xs[j];
            if (std::fabs(den) < 1.0e-20f) {
                return ys[1];
            }
            li *= (x - xs[j]) / den;
        }
        out += ys[i] * li;
    }
    return out;
}

void cubic_stencil_centered(int base_idx, int n_nodes, int out_idx[4]) {
    if (n_nodes < 4) {
        out_idx[0] = 0;
        out_idx[1] = clamp_int_local(base_idx, 0, std::max(0, n_nodes - 1));
        out_idx[2] = out_idx[1];
        out_idx[3] = out_idx[1];
        return;
    }
    if (base_idx <= 1) {
        out_idx[0] = 0;
        out_idx[1] = 1;
        out_idx[2] = 2;
        out_idx[3] = 3;
        return;
    }
    if (base_idx >= n_nodes - 3) {
        out_idx[0] = n_nodes - 4;
        out_idx[1] = n_nodes - 3;
        out_idx[2] = n_nodes - 2;
        out_idx[3] = n_nodes - 1;
        return;
    }
    out_idx[0] = base_idx - 1;
    out_idx[1] = base_idx;
    out_idx[2] = base_idx + 1;
    out_idx[3] = base_idx + 2;
}

void nc_check(int status, const std::string& where) {
    if (status != NC_NOERR) {
        throw std::runtime_error(where + ": " + nc_strerror(status));
    }
}

int get_scalar_int_var(int ncid, const char* name) {
    int varid = -1;
    nc_check(nc_inq_varid(ncid, name, &varid), std::string("nc_inq_varid(") + name + ")");
    int value = 0;
    nc_check(nc_get_var_int(ncid, varid, &value), std::string("nc_get_var_int(") + name + ")");
    return value;
}

bool get_global_int_attr(int ncid, const char* name, int& value) {
    int status = nc_get_att_int(ncid, NC_GLOBAL, name, &value);
    if (status == NC_ENOTATT) {
        return false;
    }
    nc_check(status, std::string("nc_get_att_int(") + name + ")");
    return true;
}

bool get_global_float_attr(int ncid, const char* name, float& value) {
    int status = nc_get_att_float(ncid, NC_GLOBAL, name, &value);
    if (status == NC_ENOTATT) {
        return false;
    }
    nc_check(status, std::string("nc_get_att_float(") + name + ")");
    return true;
}

bool get_global_text_attr(int ncid, const char* name, std::string& value) {
    size_t len = 0;
    int status = nc_inq_attlen(ncid, NC_GLOBAL, name, &len);
    if (status == NC_ENOTATT) {
        return false;
    }
    nc_check(status, std::string("nc_inq_attlen(") + name + ")");
    if (len == 0) {
        value.clear();
        return true;
    }
    std::vector<char> buf(len + 1, '\0');
    nc_check(nc_get_att_text(ncid, NC_GLOBAL, name, buf.data()),
             std::string("nc_get_att_text(") + name + ")");
    value = std::string(buf.data());
    return true;
}

}  // namespace

FieldData::FieldData()
    : nx(0),
      ny(0),
      nz(0),
      nzG(0),
      zperiod(1),
      divertor(0),
      ixsep(0),
      nypf1(0),
      nypf2(0),
      jyomp(0),
      dy0(0.0f),
      dz(0.0f),
      zmin(0.0f),
      zmax(2.0f * kPi),
      dz_torus(0.0f),
      xMin(0.0f),
      xMax(0.0f),
      ny_cfr(0) {}

int FieldData::lower_bracket(const std::vector<float>& xp, float x) const {
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
        int mid = (lo + hi) / 2;
        if (xp[mid] <= x) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return lo;
}

float FieldData::wrap_z(float z) const {
    const float period = zmax;
    if (period <= 0.0f) {
        return z;
    }
    float out = std::fmod(z, period);
    if (out < 0.0f) {
        out += period;
    }
    if (out >= period) {
        out -= period;
    }
    return out;
}

float FieldData::interp1(const std::vector<float>& xp,
                         const std::vector<float>& fp,
                         float x) const {
    const int n = static_cast<int>(xp.size());
    if (n == 0 || static_cast<int>(fp.size()) != n) {
        return 0.0f;
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

    int i0 = lower_bracket(xp, x);
    int i1 = i0 + 1;
    float denom = xp[i1] - xp[i0];
    if (std::fabs(denom) < 1.0e-20f) {
        return fp[i0];
    }
    float t = (x - xp[i0]) / denom;
    return fp[i0] + t * (fp[i1] - fp[i0]);
}

float FieldData::interp1_stride(const std::vector<float>& xp,
                                const float* fp,
                                int n,
                                float x,
                                int stride) const {
    if (n <= 0 || fp == 0) {
        return 0.0f;
    }
    if (n == 1) {
        return fp[0];
    }
    if (x <= xp[0]) {
        return fp[0];
    }
    if (x >= xp[n - 1]) {
        return fp[(n - 1) * stride];
    }

    int lo = 0;
    int hi = n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (xp[mid] <= x) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    float denom = xp[hi] - xp[lo];
    if (std::fabs(denom) < 1.0e-20f) {
        return fp[lo * stride];
    }
    float t = (x - xp[lo]) / denom;
    float v0 = fp[lo * stride];
    float v1 = fp[hi * stride];
    return v0 + t * (v1 - v0);
}

float FieldData::interp2_rect(const std::vector<float>& xcoords,
                              const std::vector<float>& ycoords,
                              const std::vector<float>& data,
                              int nx_local,
                              int ny_local,
                              float x,
                              float y) const {
    if (nx_local < 2 || ny_local < 2) {
        return 0.0f;
    }
    if (static_cast<int>(xcoords.size()) != nx_local ||
        static_cast<int>(ycoords.size()) != ny_local ||
        static_cast<int>(data.size()) != nx_local * ny_local) {
        return 0.0f;
    }

    int ix0 = 0;
    if (x <= xcoords[0]) {
        ix0 = 0;
    } else if (x >= xcoords[nx_local - 1]) {
        ix0 = nx_local - 2;
    } else {
        int lo = 0;
        int hi = nx_local - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (xcoords[mid] <= x) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        ix0 = lo;
    }
    int ix1 = ix0 + 1;

    int iy0 = 0;
    if (y <= ycoords[0]) {
        iy0 = 0;
    } else if (y >= ycoords[ny_local - 1]) {
        iy0 = ny_local - 2;
    } else {
        int lo = 0;
        int hi = ny_local - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (ycoords[mid] <= y) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        iy0 = lo;
    }
    int iy1 = iy0 + 1;

    float dx_loc = xcoords[ix1] - xcoords[ix0];
    float dy_loc = ycoords[iy1] - ycoords[iy0];
    float tx = 0.0f;
    float ty = 0.0f;
    if (std::fabs(dx_loc) > 1.0e-20f) {
        tx = (x - xcoords[ix0]) / dx_loc;
    }
    if (std::fabs(dy_loc) > 1.0e-20f) {
        ty = (y - ycoords[iy0]) / dy_loc;
    }
    if (tx < 0.0f) tx = 0.0f;
    if (tx > 1.0f) tx = 1.0f;
    if (ty < 0.0f) ty = 0.0f;
    if (ty > 1.0f) ty = 1.0f;

    float v00 = data[ix0 * ny_local + iy0];
    float v10 = data[ix1 * ny_local + iy0];
    float v01 = data[ix0 * ny_local + iy1];
    float v11 = data[ix1 * ny_local + iy1];

    float v0 = v00 + tx * (v10 - v00);
    float v1 = v01 + tx * (v11 - v01);
    return v0 + ty * (v1 - v0);
}

float FieldData::interp2_spline(const std::vector<float>& xcoords,
                                const std::vector<float>& ycoords,
                                const std::vector<float>& data,
                                int nx_local,
                                int ny_local,
                                float x,
                                float y) const {
    if (nx_local < 2 || ny_local < 2) {
        return 0.0f;
    }
    if (static_cast<int>(xcoords.size()) != nx_local ||
        static_cast<int>(ycoords.size()) != ny_local ||
        static_cast<int>(data.size()) != nx_local * ny_local) {
        return 0.0f;
    }

    std::vector<float> along_x(nx_local, 0.0f);
    std::vector<float> row(ny_local, 0.0f);
    NaturalCubicSpline spl;

    for (int ix = 0; ix < nx_local; ++ix) {
        const int base = ix * ny_local;
        for (int iy = 0; iy < ny_local; ++iy) {
            row[iy] = data[base + iy];
        }
        spl.build(ycoords, row);
        along_x[ix] = spl.eval(y);
    }

    spl.build(xcoords, along_x);
    return spl.eval(x);
}

float FieldData::interp_xz_3d_y(const std::vector<float>& data3d,
                                int y0,
                                float x,
                                float z) const {
    if (nx < 2 || nzG < 1) {
        return 0.0f;
    }
    if (y0 < 0) y0 = 0;
    if (y0 >= ny) y0 = ny - 1;

    int ix0 = 0;
    float tx = 0.0f;
    if (x <= xarray.front()) {
        ix0 = 0;
        tx = 0.0f;
    } else if (x >= xarray.back()) {
        ix0 = nx - 2;
        tx = 1.0f;
    } else {
        ix0 = lower_bracket(xarray, x);
        float denom = xarray[ix0 + 1] - xarray[ix0];
        tx = (std::fabs(denom) > 1.0e-20f) ? (x - xarray[ix0]) / denom : 0.0f;
    }
    if (tx < 0.0f) tx = 0.0f;
    if (tx > 1.0f) tx = 1.0f;

    float z_wrapped = wrap_z(z);
    float zf = z_wrapped / dz_torus;
    int iz0 = static_cast<int>(std::floor(zf));
    if (iz0 < 0) iz0 = 0;
    if (iz0 >= nzG) iz0 = nzG - 1;
    int iz1 = (iz0 + 1) % nzG;
    float tz = zf - static_cast<float>(iz0);
    if (tz < 0.0f) tz = 0.0f;
    if (tz > 1.0f) tz = 1.0f;

    int ix1 = ix0 + 1;
    if (ix1 >= nx) ix1 = nx - 1;

    float v00 = data3d[idx3(ix0, y0, iz0)];
    float v01 = data3d[idx3(ix0, y0, iz1)];
    float v10 = data3d[idx3(ix1, y0, iz0)];
    float v11 = data3d[idx3(ix1, y0, iz1)];

    float v0 = v00 + tz * (v01 - v00);
    float v1 = v10 + tz * (v11 - v10);
    return v0 + tx * (v1 - v0);
}

float FieldData::interp_xz_3d_y_spline(const std::vector<float>& data3d,
                                       int y0,
                                       float x,
                                       float z) const {
    if (nx < 4 || nzG < 4) {
        return interp_xz_3d_y(data3d, y0, x, z);
    }
    if (nx < 2 || nzG < 1) {
        return 0.0f;
    }
    y0 = clamp_int_local(y0, 0, ny - 1);
    float xq = x;
    if (xq <= xarray.front()) xq = xarray.front();
    if (xq >= xarray.back()) xq = xarray.back();

    int ix_base = lower_bracket(xarray, xq);
    int ixs[4] = {0, 1, 2, 3};
    cubic_stencil_centered(ix_base, nx, ixs);

    float zq = wrap_z(z);
    int iz_base = static_cast<int>(std::floor(zq / dz_torus));
    iz_base = clamp_int_local(iz_base, 0, nzG - 1);
    int izs[4] = {0, 1, 2, 3};
    cubic_stencil_centered(iz_base, nzG + 1, izs);

    auto sample = [&](int ix, int iz_ext) -> float {
        if (iz_ext <= 0) return data3d[idx3(ix, y0, 0)];
        if (iz_ext >= nzG) return data3d[idx3(ix, y0, nzG - 1)];
        return data3d[idx3(ix, y0, iz_ext)];
    };

    float xvals[4] = {xarray[ixs[0]], xarray[ixs[1]], xarray[ixs[2]], xarray[ixs[3]]};
    float zinterp[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    for (int i = 0; i < 4; ++i) {
        const int ix = ixs[i];
        float zc0 = static_cast<float>(izs[0]) * dz_torus;
        float zc1 = static_cast<float>(izs[1]) * dz_torus;
        float zc2 = static_cast<float>(izs[2]) * dz_torus;
        float zc3 = static_cast<float>(izs[3]) * dz_torus;
        float zv0 = sample(ix, izs[0]);
        float zv1 = sample(ix, izs[1]);
        float zv2 = sample(ix, izs[2]);
        float zv3 = sample(ix, izs[3]);
        zinterp[i] = lagrange4(zq, zc0, zc1, zc2, zc3, zv0, zv1, zv2, zv3);
    }

    return lagrange4(xq,
                     xvals[0], xvals[1], xvals[2], xvals[3],
                     zinterp[0], zinterp[1], zinterp[2], zinterp[3]);
}

float FieldData::interp_xz_2d(const std::vector<float>& data2d,
                              float x,
                              float z) const {
    if (nx < 2 || nzG < 1) {
        return 0.0f;
    }

    int ix0 = 0;
    float tx = 0.0f;
    if (x <= xarray.front()) {
        ix0 = 0;
        tx = 0.0f;
    } else if (x >= xarray.back()) {
        ix0 = nx - 2;
        tx = 1.0f;
    } else {
        ix0 = lower_bracket(xarray, x);
        float denom = xarray[ix0 + 1] - xarray[ix0];
        tx = (std::fabs(denom) > 1.0e-20f) ? (x - xarray[ix0]) / denom : 0.0f;
    }
    if (tx < 0.0f) tx = 0.0f;
    if (tx > 1.0f) tx = 1.0f;

    float z_wrapped = wrap_z(z);
    float zf = z_wrapped / dz_torus;
    int iz0 = static_cast<int>(std::floor(zf));
    if (iz0 < 0) iz0 = 0;
    if (iz0 >= nzG) iz0 = nzG - 1;
    int iz1 = (iz0 + 1) % nzG;
    float tz = zf - static_cast<float>(iz0);
    if (tz < 0.0f) tz = 0.0f;
    if (tz > 1.0f) tz = 1.0f;

    int ix1 = ix0 + 1;
    if (ix1 >= nx) ix1 = nx - 1;

    float v00 = data2d[idx_xz(ix0, iz0)];
    float v01 = data2d[idx_xz(ix0, iz1)];
    float v10 = data2d[idx_xz(ix1, iz0)];
    float v11 = data2d[idx_xz(ix1, iz1)];

    float v0 = v00 + tz * (v01 - v00);
    float v1 = v10 + tz * (v11 - v10);
    return v0 + tx * (v1 - v0);
}

float FieldData::interp_xz_2d_spline(const std::vector<float>& data2d,
                                     float x,
                                     float z) const {
    if (nx < 4 || nzG < 4) {
        return interp_xz_2d(data2d, x, z);
    }
    if (nx < 2 || nzG < 1) {
        return 0.0f;
    }
    float xq = x;
    if (xq <= xarray.front()) xq = xarray.front();
    if (xq >= xarray.back()) xq = xarray.back();

    int ix_base = lower_bracket(xarray, xq);
    int ixs[4] = {0, 1, 2, 3};
    cubic_stencil_centered(ix_base, nx, ixs);

    float zq = wrap_z(z);
    int iz_base = static_cast<int>(std::floor(zq / dz_torus));
    iz_base = clamp_int_local(iz_base, 0, nzG - 1);
    int izs[4] = {0, 1, 2, 3};
    cubic_stencil_centered(iz_base, nzG + 1, izs);

    auto sample = [&](int ix, int iz_ext) -> float {
        if (iz_ext <= 0) return data2d[idx_xz(ix, 0)];
        if (iz_ext >= nzG) return data2d[idx_xz(ix, nzG - 1)];
        return data2d[idx_xz(ix, iz_ext)];
    };

    float xvals[4] = {xarray[ixs[0]], xarray[ixs[1]], xarray[ixs[2]], xarray[ixs[3]]};
    float zinterp[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    for (int i = 0; i < 4; ++i) {
        const int ix = ixs[i];
        float zc0 = static_cast<float>(izs[0]) * dz_torus;
        float zc1 = static_cast<float>(izs[1]) * dz_torus;
        float zc2 = static_cast<float>(izs[2]) * dz_torus;
        float zc3 = static_cast<float>(izs[3]) * dz_torus;
        float zv0 = sample(ix, izs[0]);
        float zv1 = sample(ix, izs[1]);
        float zv2 = sample(ix, izs[2]);
        float zv3 = sample(ix, izs[3]);
        zinterp[i] = lagrange4(zq, zc0, zc1, zc2, zc3, zv0, zv1, zv2, zv3);
    }

    return lagrange4(xq,
                     xvals[0], xvals[1], xvals[2], xvals[3],
                     zinterp[0], zinterp[1], zinterp[2], zinterp[3]);
}

float FieldData::interp_periodic_row_3d(const std::vector<float>& data3d,
                                        int ix,
                                        int iy,
                                        float z) const {
    if (ix < 0) ix = 0;
    if (ix >= nx) ix = nx - 1;
    if (iy < 0) iy = 0;
    if (iy >= ny) iy = ny - 1;

    float z_wrapped = wrap_z(z);
    float zf = z_wrapped / dz_torus;
    int iz0 = static_cast<int>(std::floor(zf));
    if (iz0 < 0) iz0 = 0;
    if (iz0 >= nzG) iz0 = nzG - 1;
    int iz1 = (iz0 + 1) % nzG;
    float t = zf - static_cast<float>(iz0);
    if (t < 0.0f) t = 0.0f;
    if (t > 1.0f) t = 1.0f;

    float v0 = data3d[idx3(ix, iy, iz0)];
    float v1 = data3d[idx3(ix, iy, iz1)];
    return v0 + t * (v1 - v0);
}

float FieldData::interp_periodic_row_3d_spline(const std::vector<float>& data3d,
                                               int ix,
                                               int iy,
                                               float z) const {
    ix = clamp_int_local(ix, 0, nx - 1);
    iy = clamp_int_local(iy, 0, ny - 1);
    if (nzG < 4) {
        return interp_periodic_row_3d(data3d, ix, iy, z);
    }

    float zq = wrap_z(z);
    int iz_base = static_cast<int>(std::floor(zq / dz_torus));
    iz_base = clamp_int_local(iz_base, 0, nzG - 1);
    int izs[4] = {0, 1, 2, 3};
    cubic_stencil_centered(iz_base, nzG + 1, izs);

    auto sample = [&](int iz_ext) -> float {
        if (iz_ext <= 0) return data3d[idx3(ix, iy, 0)];
        if (iz_ext >= nzG) return data3d[idx3(ix, iy, 0)];
        return data3d[idx3(ix, iy, iz_ext)];
    };

    float zc0 = static_cast<float>(izs[0]) * dz_torus;
    float zc1 = static_cast<float>(izs[1]) * dz_torus;
    float zc2 = static_cast<float>(izs[2]) * dz_torus;
    float zc3 = static_cast<float>(izs[3]) * dz_torus;
    float zv0 = sample(izs[0]);
    float zv1 = sample(izs[1]);
    float zv2 = sample(izs[2]);
    float zv3 = sample(izs[3]);
    return lagrange4(zq, zc0, zc1, zc2, zc3, zv0, zv1, zv2, zv3);
}

void FieldData::evaluate_stage(float x,
                               float z,
                               int y_start_1b,
                               int region,
                               int direction,
                               int stage,
                               float& dxdy_val,
                               float& dzdy_val) const {
    if (stage < 0) stage = 0;
    if (stage > 2) stage = 2;

    int yp = y_start_1b - 1;
    if (yp < 0) yp = 0;
    if (yp >= ny) yp = ny - 1;

    bool use_twist = false;
    bool use_plus = false;
    int yn = yp;

    if (direction == 1) {
        if (region == 0 && y_start_1b == nypf2) {
            use_twist = true;
            use_plus = true;
        } else {
            yn = y_start_1b;
            if (yn < 0) yn = 0;
            if (yn >= ny) yn = ny - 1;
        }
    } else if (direction == -1) {
        if (region == 0 && y_start_1b == (nypf1 + 1)) {
            use_twist = true;
            use_plus = false;
        } else {
            yn = y_start_1b - 2;
            if (yn < 0) yn = 0;
            if (yn >= ny) yn = ny - 1;
        }
    } else {
        throw std::runtime_error("direction must be +1 or -1");
    }

    if (stage == 0) {
        dxdy_val = interp_xz_3d_y_spline(dxdy, yp, x, z);
        dzdy_val = interp_xz_3d_y_spline(dzdy, yp, x, z);
        return;
    }

    if (stage == 2) {
        if (use_twist) {
            if (use_plus) {
                dxdy_val = interp_xz_2d_spline(dxdy_p1, x, z);
                dzdy_val = interp_xz_2d_spline(dzdy_p1, x, z);
            } else {
                dxdy_val = interp_xz_2d_spline(dxdy_m1, x, z);
                dzdy_val = interp_xz_2d_spline(dzdy_m1, x, z);
            }
        } else {
            dxdy_val = interp_xz_3d_y_spline(dxdy, yn, x, z);
            dzdy_val = interp_xz_3d_y_spline(dzdy, yn, x, z);
        }
        return;
    }

    float dx_p = interp_xz_3d_y_spline(dxdy, yp, x, z);
    float dz_p = interp_xz_3d_y_spline(dzdy, yp, x, z);
    float dx_n = 0.0f;
    float dz_n = 0.0f;
    if (use_twist) {
        if (use_plus) {
            dx_n = interp_xz_2d_spline(dxdy_p1, x, z);
            dz_n = interp_xz_2d_spline(dzdy_p1, x, z);
        } else {
            dx_n = interp_xz_2d_spline(dxdy_m1, x, z);
            dz_n = interp_xz_2d_spline(dzdy_m1, x, z);
        }
    } else {
        dx_n = interp_xz_3d_y_spline(dxdy, yn, x, z);
        dz_n = interp_xz_3d_y_spline(dzdy, yn, x, z);
    }
    dxdy_val = 0.5f * (dx_p + dx_n);
    dzdy_val = 0.5f * (dz_p + dz_n);
}

void FieldData::load_var_2d(int ncid, const char* name, std::vector<float>& out) const {
    int varid = -1;
    nc_check(nc_inq_varid(ncid, name, &varid), std::string("nc_inq_varid(") + name + ")");

    std::vector<double> temp(static_cast<size_t>(nx) * ny);
    nc_check(nc_get_var_double(ncid, varid, temp.data()),
             std::string("nc_get_var_double(") + name + ")");

    out.resize(temp.size());
    for (size_t i = 0; i < temp.size(); ++i) {
        out[i] = static_cast<float>(temp[i]);
    }
}

void FieldData::load_var_1d(int ncid, const char* name, int len, std::vector<float>& out) const {
    int varid = -1;
    nc_check(nc_inq_varid(ncid, name, &varid), std::string("nc_inq_varid(") + name + ")");

    std::vector<double> temp(len);
    nc_check(nc_get_var_double(ncid, varid, temp.data()),
             std::string("nc_get_var_double(") + name + ")");

    out.resize(len);
    for (int i = 0; i < len; ++i) {
        out[i] = static_cast<float>(temp[i]);
    }
}

void FieldData::load_var_3d_nz(int ncid, const char* name, int nz_file, std::vector<float>& out) const {
    int varid = -1;
    nc_check(nc_inq_varid(ncid, name, &varid), std::string("nc_inq_varid(") + name + ")");

    std::vector<double> temp(static_cast<size_t>(nx) * ny * nz_file);
    nc_check(nc_get_var_double(ncid, varid, temp.data()),
             std::string("nc_get_var_double(") + name + ")");

    // apar NetCDF files are generated from MATLAB/Octave arrays. Convert
    // the incoming 3D payload from MATLAB linearization to C row-major.
    out.resize(temp.size());
    const size_t nxy = static_cast<size_t>(nx) * static_cast<size_t>(ny);
    for (size_t dst = 0; dst < temp.size(); ++dst) {
        int ix = static_cast<int>(dst % static_cast<size_t>(nx));
        size_t tmp = dst / static_cast<size_t>(nx);
        int iy = static_cast<int>(tmp % static_cast<size_t>(ny));
        int iz = static_cast<int>(dst / nxy);

        size_t src = (static_cast<size_t>(ix) * static_cast<size_t>(ny) +
                      static_cast<size_t>(iy)) * static_cast<size_t>(nz_file) +
                     static_cast<size_t>(iz);
        out[dst] = static_cast<float>(temp[src]);
    }
}

void FieldData::load_grid(const std::string& grid_file) {
    std::cout << "Loading grid information from " << grid_file << " ..." << std::endl;

    int ncid = -1;
    nc_check(nc_open(grid_file.c_str(), NC_NOWRITE, &ncid), "nc_open(grid)");

    int rxy_varid = -1;
    nc_check(nc_inq_varid(ncid, "Rxy", &rxy_varid), "nc_inq_varid(Rxy)");

    int ndims = 0;
    nc_check(nc_inq_varndims(ncid, rxy_varid, &ndims), "nc_inq_varndims(Rxy)");
    if (ndims != 2) {
        nc_close(ncid);
        throw std::runtime_error("Rxy is not 2D");
    }

    int dimids[2] = {0, 0};
    nc_check(nc_inq_vardimid(ncid, rxy_varid, dimids), "nc_inq_vardimid(Rxy)");

    size_t dim0 = 0;
    size_t dim1 = 0;
    nc_check(nc_inq_dimlen(ncid, dimids[0], &dim0), "nc_inq_dimlen(Rxy dim0)");
    nc_check(nc_inq_dimlen(ncid, dimids[1], &dim1), "nc_inq_dimlen(Rxy dim1)");
    nx = static_cast<int>(dim0);
    ny = static_cast<int>(dim1);

    load_var_2d(ncid, "Rxy", rxy);
    load_var_2d(ncid, "Zxy", zxy);
    load_var_2d(ncid, "psixy", psixy);
    load_var_2d(ncid, "zShift", zShift);
    load_var_2d(ncid, "hthe", hthe);
    load_var_2d(ncid, "Bxy", bxy);
    load_var_2d(ncid, "Btxy", btxy);
    load_var_2d(ncid, "Bpxy", bpxy);
    load_var_2d(ncid, "sinty", sinty);
    load_var_2d(ncid, "bxcvx", bxcvx);
    load_var_2d(ncid, "bxcvy", bxcvy);
    load_var_2d(ncid, "bxcvz", bxcvz);
    load_var_2d(ncid, "Jpar0", jpar0);

    std::vector<float> dy_raw;
    load_var_2d(ncid, "dy", dy_raw);
    dy0 = dy_raw.empty() ? 0.0f : dy_raw[0];

    load_var_1d(ncid, "ShiftAngle", nx, sa);

    int ixsep1 = get_scalar_int_var(ncid, "ixseps1");
    int ixsep2 = get_scalar_int_var(ncid, "ixseps2");

    if (ixsep2 < nx) {
        divertor = 2;
        ixsep = ixsep1;
        nypf1 = get_scalar_int_var(ncid, "jyseps1_1") + 1;
        nypf2 = get_scalar_int_var(ncid, "jyseps2_2") + 1;
        std::cout << "  Double null configuration" << std::endl;
    } else if (ixsep1 < nx) {
        divertor = 1;
        ixsep = ixsep1;
        nypf1 = get_scalar_int_var(ncid, "jyseps1_1") + 1;
        nypf2 = get_scalar_int_var(ncid, "jyseps2_2") + 1;
        std::cout << "  Single null configuration" << std::endl;
    } else {
        divertor = 0;
        ixsep = nx;
        nypf1 = 0;
        nypf2 = ny;
        std::cout << "  Circular configuration" << std::endl;
    }

    // Outboard midplane index (0-based): argmax R at outer radial boundary.
    jyomp = 0;
    float rmax = rxy[idx2(nx - 1, 0)];
    for (int iy = 1; iy < ny; ++iy) {
        float value = rxy[idx2(nx - 1, iy)];
        if (value > rmax) {
            rmax = value;
            jyomp = iy;
        }
    }

    xiarray.resize(nx);
    yiarray.resize(ny);
    xarray.resize(nx);
    for (int ix = 0; ix < nx; ++ix) {
        xiarray[ix] = static_cast<float>(ix + 1);
        xarray[ix] = psixy[idx2(ix, jyomp)];
    }
    for (int iy = 0; iy < ny; ++iy) {
        yiarray[iy] = static_cast<float>(iy + 1);
    }

    xMin = *std::min_element(xarray.begin(), xarray.end());
    xMax = *std::max_element(xarray.begin(), xarray.end());

    // nu = btxy * hthe / bpxy / rxy
    nu.resize(static_cast<size_t>(nx) * ny);
    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            int i2 = idx2(ix, iy);
            nu[i2] = btxy[i2] * hthe[i2] / bpxy[i2] / rxy[i2];
        }
    }

    nc_check(nc_close(ncid), "nc_close(grid)");
}

void FieldData::load_apar(const std::string& apar_file) {
    std::cout << "Loading apar data from " << apar_file << " ..." << std::endl;

    int ncid = -1;
    nc_check(nc_open(apar_file.c_str(), NC_NOWRITE, &ncid), "nc_open(apar)");

    // Read shape from variable dimensions.
    int apar_varid = -1;
    nc_check(nc_inq_varid(ncid, "apar", &apar_varid), "nc_inq_varid(apar)");
    int ndims = 0;
    nc_check(nc_inq_varndims(ncid, apar_varid, &ndims), "nc_inq_varndims(apar)");
    if (ndims != 3) {
        nc_close(ncid);
        throw std::runtime_error("apar is not 3D");
    }

    int dimids[3] = {0, 0, 0};
    nc_check(nc_inq_vardimid(ncid, apar_varid, dimids), "nc_inq_vardimid(apar)");
    size_t dim_nx = 0;
    size_t dim_ny = 0;
    size_t dim_nz = 0;
    nc_check(nc_inq_dimlen(ncid, dimids[0], &dim_nx), "nc_inq_dimlen(apar nx)");
    nc_check(nc_inq_dimlen(ncid, dimids[1], &dim_ny), "nc_inq_dimlen(apar ny)");
    nc_check(nc_inq_dimlen(ncid, dimids[2], &dim_nz), "nc_inq_dimlen(apar nz)");

    const int nx_file = static_cast<int>(dim_nx);
    const int ny_file = static_cast<int>(dim_ny);
    const int nz_file = static_cast<int>(dim_nz);

    if (nx_file != nx || ny_file != ny) {
        std::ostringstream oss;
        oss << "Dimension mismatch between apar file and grid file. "
            << "apar=(" << nx_file << "," << ny_file << "," << nz_file << "), "
            << "grid=(" << nx << "," << ny << ")";
        nc_close(ncid);
        throw std::runtime_error(oss.str());
    }

    // Optional global attributes.
    int attr_i = 0;
    if (get_global_int_attr(ncid, "zperiod", attr_i)) {
        zperiod = attr_i;
    } else {
        zperiod = 1;
    }
    if (get_global_int_attr(ncid, "divertor", attr_i)) {
        divertor = attr_i;
    }
    if (divertor == 1) {
        if (get_global_int_attr(ncid, "ixsep", attr_i)) ixsep = attr_i;
        if (get_global_int_attr(ncid, "nypf1", attr_i)) nypf1 = attr_i;
        if (get_global_int_attr(ncid, "nypf2", attr_i)) nypf2 = attr_i;
    }

    float attr_f = 0.0f;
    if (get_global_float_attr(ncid, "dy0", attr_f)) {
        dy0 = attr_f;
    }
    if (get_global_float_attr(ncid, "dz", attr_f)) {
        dz = attr_f;
    }

    std::string attr_text;
    if (get_global_text_attr(ncid, "gridfile", attr_text)) {
        apar_gridfile_attr = attr_text;
    } else {
        apar_gridfile_attr.clear();
    }

    nz = nz_file;
    nzG = nz * zperiod;
    if (nzG <= 0) {
        nc_close(ncid);
        throw std::runtime_error("Invalid nzG");
    }

    zmin = 0.0f;
    zmax = 2.0f * kPi;
    dz_torus = (zmax - zmin) / static_cast<float>(nzG);

    ziarray.resize(nzG + 1);
    zarray.resize(nzG + 1);
    for (int i = 0; i <= nzG; ++i) {
        ziarray[i] = static_cast<float>(i + 1);
        zarray[i] = static_cast<float>(i) * dz_torus;
    }

    std::vector<float> apar_file_data;
    std::vector<float> dapardx_file_data;
    std::vector<float> dapardy_file_data;
    std::vector<float> dapardz_file_data;

    load_var_3d_nz(ncid, "apar", nz_file, apar_file_data);
    load_var_3d_nz(ncid, "dapardx", nz_file, dapardx_file_data);
    load_var_3d_nz(ncid, "dapardy", nz_file, dapardy_file_data);
    load_var_3d_nz(ncid, "dapardz", nz_file, dapardz_file_data);

    nc_check(nc_close(ncid), "nc_close(apar)");

    // Fill full torus if needed.
    apar.assign(static_cast<size_t>(nx) * ny * nzG, 0.0f);
    dapardx.assign(static_cast<size_t>(nx) * ny * nzG, 0.0f);
    dapardy.assign(static_cast<size_t>(nx) * ny * nzG, 0.0f);
    dapardz.assign(static_cast<size_t>(nx) * ny * nzG, 0.0f);

    if (nz_file == nzG) {
        apar = apar_file_data;
        dapardx = dapardx_file_data;
        dapardy = dapardy_file_data;
        dapardz = dapardz_file_data;
    } else if (nz_file == nz) {
        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                for (int zp = 0; zp < zperiod; ++zp) {
                    for (int iz = 0; iz < nz; ++iz) {
                        int iz_full = zp * nz + iz;
                        int src_idx = (ix * ny + iy) * nz + iz;
                        int dst_idx = (ix * ny + iy) * nzG + iz_full;
                        apar[dst_idx] = apar_file_data[src_idx];
                        dapardx[dst_idx] = dapardx_file_data[src_idx];
                        dapardy[dst_idx] = dapardy_file_data[src_idx];
                        dapardz[dst_idx] = dapardz_file_data[src_idx];
                    }
                }
            }
        }
    } else {
        std::ostringstream oss;
        oss << "Unexpected apar z dimension: " << nz_file
            << " (expected " << nz << " or " << nzG << ")";
        throw std::runtime_error(oss.str());
    }
}

void FieldData::compute_fields() {
    std::cout << "Computing perturbed field arrays ..." << std::endl;
    dxdy.assign(static_cast<size_t>(nx) * ny * nzG, 0.0f);
    dzdy.assign(static_cast<size_t>(nx) * ny * nzG, 0.0f);

    std::vector<float> A1(static_cast<size_t>(nx) * ny, 0.0f);
    std::vector<float> A2(static_cast<size_t>(nx) * ny, 0.0f);
    std::vector<float> A3(static_cast<size_t>(nx) * ny, 0.0f);
    std::vector<float> JJ(static_cast<size_t>(nx) * ny, 0.0f);
    std::vector<float> b0dgy(static_cast<size_t>(nx) * ny, 0.0f);

    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            int i2 = idx2(ix, iy);
            A1[i2] = rxy[i2] * bpxy[i2] * btxy[i2] / hthe[i2];
            A2[i2] = bxy[i2] * bxy[i2];
            A3[i2] = sinty[i2] * A1[i2];
            JJ[i2] = kMu0 * bpxy[i2] / hthe[i2] / (bxy[i2] * bxy[i2]) * jpar0[i2];
            b0dgy[i2] = bpxy[i2] / hthe[i2];
        }
    }

    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            int i2 = idx2(ix, iy);
            for (int iz = 0; iz < nzG; ++iz) {
                int i3 = idx3(ix, iy, iz);
                float bdgx = (1.0f / bxy[i2]) *
                                 (-A1[i2] * dapardy[i3] - A2[i2] * dapardz[i3]) +
                             apar[i3] * bxcvx[i2];
                float bdgy = (1.0f / bxy[i2]) *
                                 (A1[i2] * dapardx[i3] - A3[i2] * dapardz[i3]) +
                             apar[i3] * (bxcvy[i2] + JJ[i2]);
                float bdgz = (1.0f / bxy[i2]) *
                                 (A2[i2] * dapardx[i3] + A3[i2] * dapardz[i3]) +
                             apar[i3] * bxcvz[i2];

                float denom = b0dgy[i2] + bdgy;
                dxdy[i3] = bdgx / denom;
                dzdy[i3] = bdgz / denom;
            }
        }
    }
}

void FieldData::compute_twist_shift() {
    std::cout << "Computing twist-shift boundary arrays ..." << std::endl;

    dxdy_p1.assign(static_cast<size_t>(nx) * nzG, 0.0f);
    dzdy_p1.assign(static_cast<size_t>(nx) * nzG, 0.0f);
    dxdy_m1.assign(static_cast<size_t>(nx) * nzG, 0.0f);
    dzdy_m1.assign(static_cast<size_t>(nx) * nzG, 0.0f);

    // Match Python behavior: compute only for ix < ixsep.
    for (int ix = 0; ix < ixsep; ++ix) {
        for (int iz = 0; iz < nzG; ++iz) {
            float z_base = zarray[iz];

            float zp = wrap_z(z_base + sa[ix]);
            dxdy_p1[idx_xz(ix, iz)] = interp_periodic_row_3d_spline(dxdy, ix, nypf1, zp);
            dzdy_p1[idx_xz(ix, iz)] = interp_periodic_row_3d_spline(dzdy, ix, nypf1, zp);

            float zm = wrap_z(z_base - sa[ix]);
            dxdy_m1[idx_xz(ix, iz)] = interp_periodic_row_3d_spline(dxdy, ix, nypf2 - 1, zm);
            dzdy_m1[idx_xz(ix, iz)] = interp_periodic_row_3d_spline(dzdy, ix, nypf2 - 1, zm);
        }
    }
}

void FieldData::setup_theta_and_cfr() {
    std::cout << "Setting up theta and CFR geometry ..." << std::endl;

    theta.assign(ny, 0.0f);

    if (divertor == 1) {
        // Single-null theta.
        float core_r_min = rxy[idx2(0, nypf1)];
        float core_r_max = core_r_min;
        float core_z_min = zxy[idx2(0, nypf1)];
        float core_z_max = core_z_min;
        for (int iy = nypf1; iy < ny - nypf1; ++iy) {
            float rv = rxy[idx2(0, iy)];
            float zv = zxy[idx2(0, iy)];
            if (rv < core_r_min) core_r_min = rv;
            if (rv > core_r_max) core_r_max = rv;
            if (zv < core_z_min) core_z_min = zv;
            if (zv > core_z_max) core_z_max = zv;
        }
        float center_x = 0.5f * (core_r_max + core_r_min);
        float center_y = 0.5f * (core_z_max + core_z_min);

        std::vector<int> tmp;
        tmp.reserve(ny);
        for (int iy = 0; iy < nypf1; ++iy) tmp.push_back(iy);
        for (int iy = ny - nypf1 - 1; iy >= nypf1; --iy) tmp.push_back(iy);
        for (int iy = ny - nypf1; iy < ny; ++iy) tmp.push_back(iy);

        std::vector<float> sepx(tmp.size(), 0.0f);
        std::vector<float> sepy(tmp.size(), 0.0f);
        for (size_t i = 0; i < tmp.size(); ++i) {
            int iy = tmp[i];
            sepx[i] = 0.5f * (rxy[idx2(ixsep - 1, iy)] + rxy[idx2(ixsep, iy)]);
            sepy[i] = 0.5f * (zxy[idx2(ixsep - 1, iy)] + zxy[idx2(ixsep, iy)]);
        }

        float xpoint_x = 0.25f * (sepx[nypf1 - 1] + sepx[nypf1] +
                                   sepx[ny - nypf1 - 1] + sepx[ny - nypf1]);
        float xpoint_y = 0.25f * (sepy[nypf1 - 1] + sepy[nypf1] +
                                   sepy[ny - nypf1 - 1] + sepy[ny - nypf1]);

        float ux = center_x - xpoint_x;
        float uy = center_y - xpoint_y;
        for (int iy = 0; iy < ny; ++iy) {
            float vx = center_x - rxy[idx2(0, iy)];
            float vy = center_y - zxy[idx2(0, iy)];
            float cross_norm = std::fabs(ux * vy - uy * vx);
            float dot = ux * vx + uy * vy;
            theta[iy] = std::atan2(cross_norm, dot) / kPi;
        }

        int itheta = 0;
        float thmax = theta[0];
        for (int iy = 1; iy < ny; ++iy) {
            if (theta[iy] > thmax) {
                thmax = theta[iy];
                itheta = iy;
            }
        }
        for (int iy = itheta; iy < ny; ++iy) {
            theta[iy] = 2.0f - theta[iy];
        }

        itheta = 0;
        thmax = theta[0];
        for (int iy = 1; iy < ny; ++iy) {
            if (theta[iy] > thmax) {
                thmax = theta[iy];
                itheta = iy;
            }
        }
        if (itheta != ny - 1) {
            for (int iy = itheta; iy < ny; ++iy) {
                theta[iy] = 4.0f - theta[iy];
            }
        }

        float ref = theta[nypf1];
        for (int iy = 0; iy < ny; ++iy) {
            theta[iy] -= ref;
        }
    } else {
        // Circular (and fallback for unsupported cases).
        float rmin = rxy[idx2(0, 0)];
        float rmax = rmin;
        float zmin_local = zxy[idx2(0, 0)];
        float zmax_local = zmin_local;
        for (int iy = 0; iy < ny; ++iy) {
            float rv = rxy[idx2(0, iy)];
            float zv = zxy[idx2(0, iy)];
            if (rv < rmin) rmin = rv;
            if (rv > rmax) rmax = rv;
            if (zv < zmin_local) zmin_local = zv;
            if (zv > zmax_local) zmax_local = zv;
        }
        float center_x = 0.5f * (rmax + rmin);
        float center_y = 0.5f * (zmax_local + zmin_local);
        float ux = center_x - rxy[idx2(0, 0)];
        float uy = center_y - zxy[idx2(0, 0)];

        for (int iy = 0; iy < ny; ++iy) {
            float vx = center_x - rxy[idx2(0, iy)];
            float vy = center_y - zxy[idx2(0, iy)];
            float cross_norm = std::fabs(ux * vy - uy * vx);
            float dot = ux * vx + uy * vy;
            theta[iy] = std::atan2(cross_norm, dot) / kPi;
        }

        int itheta = 0;
        float thmax = theta[0];
        for (int iy = 1; iy < ny; ++iy) {
            if (theta[iy] > thmax) {
                thmax = theta[iy];
                itheta = iy;
            }
        }
        for (int iy = itheta; iy < ny; ++iy) {
            theta[iy] = 2.0f - theta[iy];
        }

        itheta = 0;
        thmax = theta[0];
        for (int iy = 1; iy < ny; ++iy) {
            if (theta[iy] > thmax) {
                thmax = theta[iy];
                itheta = iy;
            }
        }
        if (itheta != ny - 1) {
            for (int iy = itheta; iy < ny; ++iy) {
                theta[iy] = 4.0f - theta[iy];
            }
        }
    }

    xiarray_cfr.resize(ixsep);
    for (int ix = 0; ix < ixsep; ++ix) {
        xiarray_cfr[ix] = static_cast<float>(ix + 1);
    }

    if (divertor == 0) {
        ny_cfr = ny + 1;
        yiarray_cfr.resize(ny_cfr);
        theta_cfr.assign(ny_cfr, 0.0f);
        for (int iy = 0; iy < ny_cfr; ++iy) {
            yiarray_cfr[iy] = static_cast<float>(iy + 1);
        }
        for (int iy = 0; iy < ny; ++iy) {
            theta_cfr[iy] = theta[iy];
        }
        theta_cfr[ny_cfr - 1] = 2.0f;
    } else {
        ny_cfr = nypf2 - nypf1 + 1;
        yiarray_cfr.resize(ny_cfr);
        theta_cfr.assign(ny_cfr, 0.0f);
        for (int iy = 0; iy < ny_cfr; ++iy) {
            yiarray_cfr[iy] = static_cast<float>(nypf1 + 1 + iy);
            int src_iy = nypf1 + iy;
            if (src_iy >= 0 && src_iy < ny) {
                theta_cfr[iy] = theta[src_iy];
            }
        }
        if (!theta_cfr.empty()) {
            theta_cfr[ny_cfr - 1] = 2.0f;
        }
    }

    rxy_cfr.assign(static_cast<size_t>(ixsep) * ny_cfr, 0.0f);
    zxy_cfr.assign(static_cast<size_t>(ixsep) * ny_cfr, 0.0f);
    zs_cfr.assign(static_cast<size_t>(ixsep) * ny_cfr, 0.0f);

    for (int ix = 0; ix < ixsep; ++ix) {
        // Base slice nypf1 : nypf2 (exclusive) then append wrapped point.
        for (int j = 0; j < ny_cfr - 1; ++j) {
            int iy = nypf1 + j;
            rxy_cfr[idx_cfr(ix, j)] = rxy[idx2(ix, iy)];
            zxy_cfr[idx_cfr(ix, j)] = zxy[idx2(ix, iy)];
            zs_cfr[idx_cfr(ix, j)] = zShift[idx2(ix, iy)];
        }

        rxy_cfr[idx_cfr(ix, ny_cfr - 1)] = rxy_cfr[idx_cfr(ix, 0)];
        zxy_cfr[idx_cfr(ix, ny_cfr - 1)] = zxy_cfr[idx_cfr(ix, 0)];

        float zs_last = 0.5f * (nu[idx2(ix, nypf1)] + nu[idx2(ix, nypf2 - 1)]) * dy0;
        zs_last += zs_cfr[idx_cfr(ix, ny_cfr - 2)];
        zs_cfr[idx_cfr(ix, ny_cfr - 1)] = zs_last;
    }
}

void FieldData::initialize(const std::string& grid_file, const std::string& apar_file) {
    load_grid(grid_file);
    load_apar(apar_file);
    compute_fields();
    compute_twist_shift();
    setup_theta_and_cfr();
}
