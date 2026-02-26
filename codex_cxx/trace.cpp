#include "trace.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "rk4.h"

namespace {

const float kPi = 3.14159265358979323846f;

inline int clamp_int(int v, int lo, int hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

inline size_t traj_index(int row, int col, int nsteps) {
    return static_cast<size_t>(row) * static_cast<size_t>(nsteps) + static_cast<size_t>(col);
}

inline float get_traj(const std::vector<float>& traj, int row, int col, int nsteps) {
    return traj[traj_index(row, col, nsteps)];
}

inline void set_traj(std::vector<float>& traj, int row, int col, int nsteps, float value) {
    traj[traj_index(row, col, nsteps)] = value;
}

float interp_xindex_2d(const FieldData& field,
                       const std::vector<float>& data2d,
                       int yidx,
                       float xind_1b) {
    yidx = clamp_int(yidx, 0, field.ny - 1);
    if (field.nx < 2) {
        return data2d[field.idx2(0, yidx)];
    }

    float xf = xind_1b - 1.0f;
    if (xf <= 0.0f) {
        return data2d[field.idx2(0, yidx)];
    }
    if (xf >= static_cast<float>(field.nx - 1)) {
        return data2d[field.idx2(field.nx - 1, yidx)];
    }

    int ix0 = static_cast<int>(std::floor(xf));
    ix0 = clamp_int(ix0, 0, field.nx - 2);
    int ix1 = ix0 + 1;
    float tx = xf - static_cast<float>(ix0);

    float v0 = data2d[field.idx2(ix0, yidx)];
    float v1 = data2d[field.idx2(ix1, yidx)];
    return v0 + tx * (v1 - v0);
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

    inline float eval_segment(int seg, float t) const {
        return y[seg] + b[seg] * t + c[seg] * t * t + d[seg] * t * t * t;
    }
};

void add_root_if_in_unit_interval(std::vector<float>& roots, float r) {
    const float eps = 1.0e-6f;
    if (r < -eps || r > 1.0f + eps) return;
    if (r < 0.0f) r = 0.0f;
    if (r > 1.0f) r = 1.0f;
    for (size_t i = 0; i < roots.size(); ++i) {
        if (std::fabs(roots[i] - r) < 1.0e-5f) return;
    }
    roots.push_back(r);
}

std::vector<float> solve_segment_roots(float a0, float a1, float a2, float a3) {
    std::vector<float> roots;
    const float eps = 1.0e-12f;

    if (std::fabs(a3) < eps) {
        if (std::fabs(a2) < eps) {
            if (std::fabs(a1) < eps) return roots;
            add_root_if_in_unit_interval(roots, -a0 / a1);
            return roots;
        }
        float disc = a1 * a1 - 4.0f * a2 * a0;
        if (disc < 0.0f) return roots;
        float sq = std::sqrt(std::max(0.0f, disc));
        add_root_if_in_unit_interval(roots, (-a1 + sq) / (2.0f * a2));
        add_root_if_in_unit_interval(roots, (-a1 - sq) / (2.0f * a2));
        return roots;
    }

    float A = a2 / a3;
    float B = a1 / a3;
    float C = a0 / a3;

    float p = B - A * A / 3.0f;
    float q = 2.0f * A * A * A / 27.0f - A * B / 3.0f + C;
    float disc = 0.25f * q * q + (p * p * p) / 27.0f;

    if (disc > eps) {
        float sq = std::sqrt(disc);
        float u = std::cbrt(-0.5f * q + sq);
        float v = std::cbrt(-0.5f * q - sq);
        float y = u + v;
        add_root_if_in_unit_interval(roots, y - A / 3.0f);
    } else if (std::fabs(disc) <= eps) {
        float u = std::cbrt(-0.5f * q);
        add_root_if_in_unit_interval(roots, 2.0f * u - A / 3.0f);
        add_root_if_in_unit_interval(roots, -u - A / 3.0f);
    } else {
        float phi = std::acos(std::max(-1.0f, std::min(1.0f, (-0.5f * q) / std::sqrt(-(p * p * p) / 27.0f))));
        float m = 2.0f * std::sqrt(-p / 3.0f);
        add_root_if_in_unit_interval(roots, m * std::cos(phi / 3.0f) - A / 3.0f);
        add_root_if_in_unit_interval(roots, m * std::cos((phi + 2.0f * kPi) / 3.0f) - A / 3.0f);
        add_root_if_in_unit_interval(roots, m * std::cos((phi + 4.0f * kPi) / 3.0f) - A / 3.0f);
    }

    std::sort(roots.begin(), roots.end());
    return roots;
}

struct CrossingEval {
    float xind;
    float yind;
    float zvalue;
    float ipx;
    float ipy;
    float ipz;
};

CrossingEval evaluate_crossing(const FieldData& field,
                               const std::vector<float>& traj,
                               int nsteps,
                               int tc0,
                               int tc1,
                               int direction,
                               float alpha) {
    if (alpha < 0.0f) alpha = 0.0f;
    if (alpha > 1.0f) alpha = 1.0f;
    const float beta = 1.0f - alpha;

    float xind_tmp = beta * get_traj(traj, 1, tc0, nsteps) + alpha * get_traj(traj, 1, tc1, nsteps);
    float yind_tmp = beta * get_traj(traj, 2, tc0, nsteps) + alpha * get_traj(traj, 2, tc1, nsteps);

    float zvalue = beta * get_traj(traj, 6, tc0, nsteps) + alpha * get_traj(traj, 6, tc1, nsteps);
    if (std::fabs(get_traj(traj, 6, tc0, nsteps) - get_traj(traj, 6, tc1, nsteps)) > 1.0f) {
        float z0w = field.wrap_z(get_traj(traj, 6, tc0, nsteps));
        float z1w = field.wrap_z(get_traj(traj, 6, tc1, nsteps));
        zvalue = beta * z0w + alpha * z1w;
    }

    if (static_cast<int>(std::round(get_traj(traj, 2, tc0, nsteps))) == field.nypf2 &&
        direction == 1 &&
        xind_tmp < static_cast<float>(field.ixsep) + 0.5f) {
        yind_tmp = beta * get_traj(traj, 2, tc0, nsteps) + alpha * static_cast<float>(field.nypf2 + 1);
    } else if (static_cast<int>(std::round(get_traj(traj, 2, tc0, nsteps))) == (field.nypf1 + 1) &&
               direction == -1 &&
               xind_tmp < static_cast<float>(field.ixsep) + 0.5f) {
        yind_tmp = beta * static_cast<float>(field.nypf2 + 1) + alpha * get_traj(traj, 2, tc1, nsteps);
        float shiftangle = field.interp1(field.xiarray, field.sa, xind_tmp);
        zvalue = field.wrap_z(zvalue - shiftangle);
    } else if (tc0 > 0) {
        int y_prev = static_cast<int>(std::round(get_traj(traj, 2, tc0 - 1, nsteps)));
        if (y_prev == field.nypf2 || y_prev == (field.nypf1 + 1)) {
            float z0 = field.interp1(field.ziarray, field.zarray, get_traj(traj, 3, tc0, nsteps));
            float z1 = field.interp1(field.ziarray, field.zarray, get_traj(traj, 3, tc1, nsteps));
            zvalue = beta * z0 + alpha * z1;
        }
    }
    zvalue = field.wrap_z(zvalue);

    float rxyvalue = 0.0f;
    float zxyvalue = 0.0f;
    float zsvalue = 0.0f;
    if (xind_tmp < static_cast<float>(field.ixsep) + 0.5f) {
        rxyvalue = field.interp2_spline(field.xiarray_cfr, field.yiarray_cfr, field.rxy_cfr,
                                        field.ixsep, field.ny_cfr, xind_tmp, yind_tmp);
        zxyvalue = field.interp2_spline(field.xiarray_cfr, field.yiarray_cfr, field.zxy_cfr,
                                        field.ixsep, field.ny_cfr, xind_tmp, yind_tmp);
        zsvalue = field.interp2_spline(field.xiarray_cfr, field.yiarray_cfr, field.zs_cfr,
                                       field.ixsep, field.ny_cfr, xind_tmp, yind_tmp);
    } else {
        rxyvalue = field.interp2_spline(field.xiarray, field.yiarray, field.rxy,
                                        field.nx, field.ny, xind_tmp, yind_tmp);
        zxyvalue = field.interp2_spline(field.xiarray, field.yiarray, field.zxy,
                                        field.nx, field.ny, xind_tmp, yind_tmp);
        zsvalue = field.interp2_spline(field.xiarray, field.yiarray, field.zShift,
                                       field.nx, field.ny, xind_tmp, yind_tmp);
    }

    float ipx3d_tmp = rxyvalue * std::cos(zsvalue);
    float ipy3d_tmp = rxyvalue * std::sin(zsvalue);
    float ipx = ipx3d_tmp * std::cos(zvalue) - ipy3d_tmp * std::sin(zvalue);
    float ipy = ipx3d_tmp * std::sin(zvalue) + ipy3d_tmp * std::cos(zvalue);

    CrossingEval out;
    out.xind = xind_tmp;
    out.yind = yind_tmp;
    out.zvalue = zvalue;
    out.ipx = ipx;
    out.ipy = ipy;
    out.ipz = zxyvalue;
    return out;
}

Point3D reconstruct_trajectory_point(const FieldData& field,
                                     const std::vector<float>& traj,
                                     int nsteps,
                                     int istep) {
    int yidx = static_cast<int>(std::round(get_traj(traj, 2, istep, nsteps))) - 1;
    yidx = clamp_int(yidx, 0, field.ny - 1);
    float xind_step = get_traj(traj, 1, istep, nsteps);

    float rxyvalue = interp_xindex_2d(field, field.rxy, yidx, xind_step);
    float zsvalue = interp_xindex_2d(field, field.zShift, yidx, xind_step);
    float zxyvalue = interp_xindex_2d(field, field.zxy, yidx, xind_step);
    float zvalue = field.interp1(field.ziarray, field.zarray, get_traj(traj, 3, istep, nsteps));

    float x3d_tmp = rxyvalue * std::cos(zsvalue);
    float y3d_tmp = rxyvalue * std::sin(zsvalue);

    Point3D p;
    p.x = x3d_tmp * std::cos(zvalue) - y3d_tmp * std::sin(zvalue);
    p.y = x3d_tmp * std::sin(zvalue) + y3d_tmp * std::cos(zvalue);
    p.z = zxyvalue;
    return p;
}

LineTraceResult trace_single_initial_point(const FieldData& field,
                                           const TraceOptions& options,
                                           const InitialPoint& init_point) {
    const int nsteps = options.nturns * field.ny;
    const int np_max = 1250;

    float xind = init_point.xind;
    float x_start = init_point.x_start;
    int y_start = init_point.y_start;
    float z_start = init_point.z_start;

    std::vector<float> traj(static_cast<size_t>(7) * nsteps, 0.0f);

    int it = 0;
    int iturn = 1;

    int region = 1;
    if (xind < static_cast<float>(field.ixsep) + 0.5f) {
        region = 0;
        if (y_start < field.nypf1 + 1 || y_start > field.nypf2) {
            region = 2;
        }
    }

    float zind = field.interp1(field.zarray, field.ziarray, z_start);

    std::cout << "  line " << init_point.iline << " started at indices ("
              << xind << "," << y_start << "," << zind << ")" << std::endl;

    if (field.divertor == 1) {
        if (y_start == field.ny && options.direction == 1) {
            std::cout << "    line starts on divertor" << std::endl;
            region = 14;
            set_traj(traj, 5, it, nsteps, 0.0f);
        } else if (y_start == 1 && options.direction == -1) {
            std::cout << "    line starts on divertor" << std::endl;
            region = 13;
            set_traj(traj, 5, it, nsteps, 0.0f);
        }
    }

    set_traj(traj, 0, it, nsteps, 1.0f);
    set_traj(traj, 1, it, nsteps, xind);
    set_traj(traj, 2, it, nsteps, static_cast<float>(y_start));
    set_traj(traj, 3, it, nsteps, zind);
    set_traj(traj, 4, it, nsteps, static_cast<float>(region));
    set_traj(traj, 6, it, nsteps, z_start);

    while (region < 10 && iturn < options.nturns) {
        if (iturn % 50 == 1) {
            std::cout << "    line " << init_point.iline << ", turn "
                      << iturn << "/" << options.nturns << std::endl;
        }

        for (int iy = 0; iy < field.ny - 1; ++iy) {
            float x_end = x_start;
            float z_end = z_start;
            int y_end = y_start;

            if (region == 0 && y_start > field.nypf1 && y_start < field.nypf2 + 1) {
                rk4_step(field, x_start, y_start, z_start, region, options.direction, x_end, z_end);
                y_end = (options.direction == 1) ? (y_start + 1) : (y_start - 1);

                set_traj(traj, 6, it + 1, nsteps, z_end);

                if (x_end > field.xMax) {
                    std::cout << "    line " << init_point.iline << " reaches outer boundary" << std::endl;
                    region = 12;
                } else if (x_end < field.xMin) {
                    std::cout << "    line " << init_point.iline << " reaches inner boundary" << std::endl;
                    region = 11;
                } else {
                    xind = field.interp1(field.xarray, field.xiarray, x_end);
                    if (xind > static_cast<float>(field.ixsep) + 0.5f) {
                        region = 1;
                        std::cout << "    line " << init_point.iline << " enters SOL" << std::endl;
                    }
                }

                if (options.direction == 1 && y_start == field.nypf2 && region == 0) {
                    float shiftangle = field.interp1(field.xiarray, field.sa, xind);
                    z_end += shiftangle;
                    y_end = field.nypf1 + 1;
                }
                if (options.direction == -1 && y_start == field.nypf1 + 1 && region == 0) {
                    float shiftangle = field.interp1(field.xiarray, field.sa, xind);
                    z_end -= shiftangle;
                    y_end = field.nypf2;
                }

                z_end = field.wrap_z(z_end);
                zind = field.interp1(field.zarray, field.ziarray, z_end);

                ++it;
                if (it >= nsteps) {
                    region = 15;
                    break;
                }

                int ixh = static_cast<int>(std::round(xind));
                ixh = clamp_int(ixh, 1, field.nx);
                int iyh = clamp_int(y_end, 1, field.ny);

                set_traj(traj, 0, it, nsteps, static_cast<float>(iturn));
                set_traj(traj, 1, it, nsteps, xind);
                set_traj(traj, 2, it, nsteps, static_cast<float>(y_end));
                set_traj(traj, 3, it, nsteps, zind);
                set_traj(traj, 4, it, nsteps, static_cast<float>(region));
                set_traj(traj, 5, it, nsteps, field.hthe[field.idx2(ixh - 1, iyh - 1)]);

                x_start = x_end;
                y_start = y_end;
                z_start = z_end;
            } else if (region == 1 || region == 2) {
                rk4_step(field, x_start, y_start, z_start, region, options.direction, x_end, z_end);
                y_end = (options.direction == 1) ? (y_start + 1) : (y_start - 1);

                set_traj(traj, 6, it + 1, nsteps, z_end);

                if (options.direction == 1 && y_start == field.nypf1 && region == 2) {
                    y_end = field.nypf2 + 1;
                } else if (options.direction == -1 && y_start == field.nypf2 + 1 && region == 2) {
                    y_end = field.nypf1;
                }

                if (x_end > field.xMax) {
                    std::cout << "    line " << init_point.iline << " reaches outer boundary" << std::endl;
                    region = 12;
                } else if (x_end < field.xMin) {
                    std::cout << "    line " << init_point.iline << " reaches inner boundary" << std::endl;
                    region = 11;
                } else {
                    xind = field.interp1(field.xarray, field.xiarray, x_end);
                    if (xind < static_cast<float>(field.ixsep) + 0.5f &&
                        y_end > field.nypf1 && y_end < field.nypf2 + 1) {
                        region = 0;
                    } else if (xind < static_cast<float>(field.ixsep) + 0.5f &&
                               (y_end > field.nypf2 - 1 || y_end < field.nypf1)) {
                        region = 2;
                    }
                }

                if (options.direction == 1 && y_end == field.ny) {
                    region = 14;
                } else if (options.direction == -1 && y_end == 1) {
                    region = 13;
                }

                z_end = field.wrap_z(z_end);
                zind = field.interp1(field.zarray, field.ziarray, z_end);

                ++it;
                if (it >= nsteps) {
                    region = 15;
                    break;
                }

                int ixh = static_cast<int>(std::round(xind));
                ixh = clamp_int(ixh, 1, field.nx);
                int iyh = clamp_int(y_start, 1, field.ny);

                set_traj(traj, 0, it, nsteps, static_cast<float>(iturn));
                set_traj(traj, 1, it, nsteps, xind);
                set_traj(traj, 2, it, nsteps, static_cast<float>(y_end));
                set_traj(traj, 3, it, nsteps, zind);
                set_traj(traj, 4, it, nsteps, static_cast<float>(region));
                set_traj(traj, 5, it, nsteps, field.hthe[field.idx2(ixh - 1, iyh - 1)]);

                x_start = x_end;
                y_start = y_end;
                z_start = z_end;
            }

            if (region >= 10 || it + 1 >= nsteps) {
                break;
            }
        }

        ++iturn;
    }

    int itmax = it + 1;
    if (itmax < 1) {
        itmax = 1;
    }
    if (itmax > nsteps) {
        itmax = nsteps;
    }

    LineTraceResult result;
    result.iline = init_point.iline;
    result.end_region = region;
    result.connection_length = 0.0f;
    result.trajectory_xyz.reserve(itmax);
    result.puncture_steps.reserve(np_max);
    result.puncture_xyz.reserve(np_max);
    result.puncture_theta_psi.reserve(np_max);

    for (int istep = 0; istep < itmax; ++istep) {
        result.trajectory_xyz.push_back(reconstruct_trajectory_point(field, traj, nsteps, istep));
    }

    if (itmax > 1) {
        std::vector<float> fit_x(static_cast<size_t>(itmax), 0.0f);
        std::vector<float> fit_it(static_cast<size_t>(itmax), 0.0f);
        for (int istep = 0; istep < itmax; ++istep) {
            fit_x[static_cast<size_t>(istep)] = result.trajectory_xyz[static_cast<size_t>(istep)].x;
            fit_it[static_cast<size_t>(istep)] = static_cast<float>(istep + 1);
        }

        NaturalCubicSpline fit_spline;
        fit_spline.build(fit_it, fit_x);

        float last_fit_root = -1.0e30f;
        const float endpoint_eps = 1.0e-6f;
        const float dedup_eps = 1.0e-5f;

        for (int seg = 0; seg < itmax - 1 && static_cast<int>(result.puncture_xyz.size()) < np_max; ++seg) {
            std::vector<float> roots = solve_segment_roots(fit_spline.y[seg],
                                                           fit_spline.b[seg],
                                                           fit_spline.c[seg],
                                                           fit_spline.d[seg]);

            if (roots.empty()) {
                float x0 = fit_x[static_cast<size_t>(seg)];
                float x1 = fit_x[static_cast<size_t>(seg + 1)];
                if (x0 * x1 <= 0.0f) {
                    float alpha = 0.5f;
                    float denom = x1 - x0;
                    if (std::fabs(denom) > 1.0e-20f) {
                        alpha = -x0 / denom;
                    }
                    if (alpha < 0.0f) alpha = 0.0f;
                    if (alpha > 1.0f) alpha = 1.0f;
                    roots.push_back(alpha);
                }
            }

            for (size_t ir = 0; ir < roots.size() && static_cast<int>(result.puncture_xyz.size()) < np_max; ++ir) {
                float alpha = roots[ir];
                if (alpha <= endpoint_eps && seg > 0) {
                    continue;
                }
                if (alpha >= 1.0f - endpoint_eps && seg < itmax - 2) {
                    continue;
                }

                float fit_root = static_cast<float>(seg + 1) + alpha;
                if (std::fabs(fit_root - last_fit_root) < dedup_eps) {
                    continue;
                }
                last_fit_root = fit_root;

                int tc0 = seg;
                int tc1 = seg + 1;
                CrossingEval cross =
                    evaluate_crossing(field, traj, nsteps, tc0, tc1, options.direction, alpha);

                if (cross.ipy > 0.0f) {
                    Point3D puncture_xyz;
                    puncture_xyz.x = cross.ipx;
                    puncture_xyz.y = cross.ipy;
                    puncture_xyz.z = cross.ipz;

                    Point2D theta_psi;
                    theta_psi.x = field.interp1(field.yiarray_cfr, field.theta_cfr, cross.yind);
                    theta_psi.y = field.interp1(field.xiarray, field.xarray, cross.xind);

                    int step_idx = static_cast<int>(std::floor(fit_root));
                    if (step_idx < 1) step_idx = 1;
                    if (step_idx >= itmax) step_idx = itmax - 1;

                    result.puncture_steps.push_back(step_idx);
                    result.puncture_xyz.push_back(puncture_xyz);
                    result.puncture_theta_psi.push_back(theta_psi);
                }
            }
        }
    }

    float lc = 0.0f;
    for (int istep = 0; istep < itmax; ++istep) {
        lc += get_traj(traj, 5, istep, nsteps);
    }
    result.connection_length = lc;

    std::cout << "    line " << init_point.iline << ": "
              << result.puncture_xyz.size()
              << " punctures, connection length=" << result.connection_length
              << ", end region=" << result.end_region << std::endl;

    return result;
}

void write_trace_outputs(const std::vector<LineTraceResult>& results) {
    std::ofstream ip_xyz("ip_xyz.txt");
    std::ofstream ip_thetapsi("ip_thetapsi.txt");
    std::ofstream traj_xyz("traj_xyz.txt");
    if (!ip_xyz || !ip_thetapsi || !traj_xyz) {
        throw std::runtime_error("failed to open one or more output files");
    }

    ip_xyz << "iline it ipx ipy ipz\n";
    ip_thetapsi << "iline it theta psi\n";
    traj_xyz << "iline it x y z\n";

    ip_xyz << std::setprecision(16);
    ip_thetapsi << std::setprecision(16);
    traj_xyz << std::setprecision(16);

    for (size_t ir = 0; ir < results.size(); ++ir) {
        const LineTraceResult& result = results[ir];

        for (size_t istep = 0; istep < result.trajectory_xyz.size(); ++istep) {
            const Point3D& p = result.trajectory_xyz[istep];
            traj_xyz << result.iline << " " << (istep + 1) << " "
                     << p.x << " " << p.y << " " << p.z << "\n";
        }

        size_t np = result.puncture_xyz.size();
        np = std::min(np, result.puncture_theta_psi.size());
        np = std::min(np, result.puncture_steps.size());
        for (size_t ip = 0; ip < np; ++ip) {
            const Point3D& pxyz = result.puncture_xyz[ip];
            const Point2D& tpsi = result.puncture_theta_psi[ip];
            int step_idx = result.puncture_steps[ip];

            ip_xyz << result.iline << " " << step_idx << " "
                   << pxyz.x << " " << pxyz.y << " " << pxyz.z << "\n";

            float theta_rad = tpsi.x * kPi;
            ip_thetapsi << result.iline << " " << step_idx << " "
                        << theta_rad << " " << tpsi.y << "\n";
        }
    }
}

}  // namespace

std::vector<InitialPoint> build_initial_points(const FieldData& field, const TraceOptions& options) {
    std::vector<float> lines = options.lines;
    if (lines.empty()) {
        lines.reserve(options.nlines);
        for (int i = 1; i <= options.nlines; ++i) {
            lines.push_back(static_cast<float>(i));
        }
    }

    std::vector<InitialPoint> initial_points;
    initial_points.reserve(lines.size());

    const int default_y_start = field.jyomp + 1;
    const float default_z_start = field.zarray.empty() ? 0.0f : field.zarray[0];

    for (size_t i = 0; i < lines.size(); ++i) {
        float iline = lines[i];
        if (iline < 1.0f || iline > static_cast<float>(field.nx)) {
            std::cout << "  Skipping invalid line index " << iline << std::endl;
            continue;
        }

        InitialPoint point;
        point.iline = iline;
        point.xind = iline;
        point.x_start = field.interp1(field.xiarray, field.xarray, point.xind);
        point.y_start = default_y_start;
        point.z_start = default_z_start;
        initial_points.push_back(point);
    }

    return initial_points;
}

std::vector<LineTraceResult> trace_initial_points(const FieldData& field,
                                                  const TraceOptions& options,
                                                  const std::vector<InitialPoint>& initial_points) {
    std::vector<LineTraceResult> results;
    results.reserve(initial_points.size());

    for (size_t i = 0; i < initial_points.size(); ++i) {
        results.push_back(trace_single_initial_point(field, options, initial_points[i]));
    }

    return results;
}

void trace_field_lines(const FieldData& field, const TraceOptions& options) {
    if (options.direction != 1 && options.direction != -1) {
        throw std::runtime_error("direction must be 1 or -1");
    }
    if (options.nturns <= 0 || options.nlines <= 0) {
        throw std::runtime_error("nturns and nlines must be > 0");
    }

    std::cout << "Starting field-line tracing ..." << std::endl;

    std::vector<InitialPoint> initial_points = build_initial_points(field, options);
    std::vector<LineTraceResult> results = trace_initial_points(field, options, initial_points);

    write_trace_outputs(results);

    std::cout << "Field-line tracing complete." << std::endl;
}
