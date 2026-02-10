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

}  // namespace

void trace_field_lines(const FieldData& field, const TraceOptions& options) {
    if (options.direction != 1 && options.direction != -1) {
        throw std::runtime_error("direction must be 1 or -1");
    }
    if (options.nturns <= 0 || options.nlines <= 0) {
        throw std::runtime_error("nturns and nlines must be > 0");
    }

    std::vector<int> lines = options.lines;
    if (lines.empty()) {
        lines.reserve(options.nlines);
        for (int i = 1; i <= options.nlines; ++i) {
            lines.push_back(i);
        }
    }

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

    const int nsteps = options.nturns * field.ny;
    const int np_max = 1250;

    std::cout << "Starting field-line tracing ..." << std::endl;

    for (size_t iline_i = 0; iline_i < lines.size(); ++iline_i) {
        const int iline = lines[iline_i];
        if (iline < 1 || iline > field.nx) {
            std::cout << "  Skipping invalid line index " << iline << std::endl;
            continue;
        }

        float xind = static_cast<float>(iline);
        float xStart = field.psixy[field.idx2(iline - 1, field.jyomp)];
        int yyy = field.jyomp + 1;
        int yStart = yyy;
        int zzz = 1;
        float zStart = field.zarray[0];

        std::vector<float> traj(static_cast<size_t>(7) * nsteps, 0.0f);
        std::vector<float> fl_x3d(nsteps, 0.0f);
        std::vector<float> fl_y3d(nsteps, 0.0f);
        std::vector<float> fl_z3d(nsteps, 0.0f);

        std::vector<float> px(np_max, 0.0f);
        std::vector<float> py(np_max, 0.0f);
        std::vector<float> pz(np_max, 0.0f);
        std::vector<float> ptheta(np_max, 0.0f);
        std::vector<float> ppsi(np_max, 0.0f);

        int it = 0;
        int iturn = 1;

        int region = 1;
        if (xind < static_cast<float>(field.ixsep) + 0.5f) {
            region = 0;
            if (yStart < field.nypf1 + 1 || yStart > field.nypf2) {
                region = 2;
            }
        }

        float zind = field.interp1(field.zarray, field.ziarray, zStart);

        std::cout << "  line " << iline << " started at indices ("
                  << xind << "," << yStart << "," << zind << ")" << std::endl;

        if (field.divertor == 1) {
            if (yStart == field.ny && options.direction == 1) {
                std::cout << "    line starts on divertor" << std::endl;
                region = 14;
                set_traj(traj, 5, it, nsteps, 0.0f);
            } else if (yStart == 1 && options.direction == -1) {
                std::cout << "    line starts on divertor" << std::endl;
                region = 13;
                set_traj(traj, 5, it, nsteps, 0.0f);
            }
        }

        set_traj(traj, 0, it, nsteps, 1.0f);
        set_traj(traj, 1, it, nsteps, xind);
        set_traj(traj, 2, it, nsteps, static_cast<float>(yStart));
        set_traj(traj, 3, it, nsteps, zind);
        set_traj(traj, 4, it, nsteps, static_cast<float>(region));
        set_traj(traj, 6, it, nsteps, zStart);

        while (region < 10 && iturn < options.nturns) {
            if (iturn % 50 == 1) {
                std::cout << "    line " << iline << ", turn "
                          << iturn << "/" << options.nturns << std::endl;
            }

            for (int iy = 0; iy < field.ny - 1; ++iy) {
                float xEnd = xStart;
                float zEnd = zStart;
                int yEnd = yStart;

                if (region == 0 && yStart > field.nypf1 && yStart < field.nypf2 + 1) {
                    rk4_step(field, xStart, yStart, zStart, region, options.direction, xEnd, zEnd);
                    yEnd = (options.direction == 1) ? (yStart + 1) : (yStart - 1);

                    set_traj(traj, 6, it + 1, nsteps, zEnd);

                    if (xEnd > field.xMax) {
                        std::cout << "    line " << iline << " reaches outer boundary" << std::endl;
                        region = 12;
                    } else if (xEnd < field.xMin) {
                        std::cout << "    line " << iline << " reaches inner boundary" << std::endl;
                        region = 11;
                    } else {
                        xind = field.interp1(field.xarray, field.xiarray, xEnd);
                        if (xind > static_cast<float>(field.ixsep) + 0.5f) {
                            region = 1;
                            std::cout << "    line " << iline << " enters SOL" << std::endl;
                        }
                    }

                    if (options.direction == 1 && yStart == field.nypf2 && region == 0) {
                        float shiftangle = field.interp1(field.xiarray, field.sa, xind);
                        zEnd += shiftangle;
                        yEnd = field.nypf1 + 1;
                    }
                    if (options.direction == -1 && yStart == field.nypf1 + 1 && region == 0) {
                        float shiftangle = field.interp1(field.xiarray, field.sa, xind);
                        zEnd -= shiftangle;
                        yEnd = field.nypf2;
                    }

                    zEnd = field.wrap_z(zEnd);
                    zind = field.interp1(field.zarray, field.ziarray, zEnd);

                    ++it;
                    if (it >= nsteps) {
                        region = 15;
                        break;
                    }

                    int ixh = static_cast<int>(std::round(xind));
                    ixh = clamp_int(ixh, 1, field.nx);
                    int iyh = clamp_int(yEnd, 1, field.ny);

                    set_traj(traj, 0, it, nsteps, static_cast<float>(iturn));
                    set_traj(traj, 1, it, nsteps, xind);
                    set_traj(traj, 2, it, nsteps, static_cast<float>(yEnd));
                    set_traj(traj, 3, it, nsteps, zind);
                    set_traj(traj, 4, it, nsteps, static_cast<float>(region));
                    set_traj(traj, 5, it, nsteps, field.hthe[field.idx2(ixh - 1, iyh - 1)]);

                    xStart = xEnd;
                    yStart = yEnd;
                    zStart = zEnd;
                } else if (region == 1 || region == 2) {
                    rk4_step(field, xStart, yStart, zStart, region, options.direction, xEnd, zEnd);
                    yEnd = (options.direction == 1) ? (yStart + 1) : (yStart - 1);

                    set_traj(traj, 6, it + 1, nsteps, zEnd);

                    if (options.direction == 1 && yStart == field.nypf1 && region == 2) {
                        yEnd = field.nypf2 + 1;
                    } else if (options.direction == -1 && yStart == field.nypf2 + 1 && region == 2) {
                        yEnd = field.nypf1;
                    }

                    if (xEnd > field.xMax) {
                        std::cout << "    line " << iline << " reaches outer boundary" << std::endl;
                        region = 12;
                    } else if (xEnd < field.xMin) {
                        std::cout << "    line " << iline << " reaches inner boundary" << std::endl;
                        region = 11;
                    } else {
                        xind = field.interp1(field.xarray, field.xiarray, xEnd);
                        if (xind < static_cast<float>(field.ixsep) + 0.5f &&
                            yEnd > field.nypf1 && yEnd < field.nypf2 + 1) {
                            region = 0;
                        } else if (xind < static_cast<float>(field.ixsep) + 0.5f &&
                                   (yEnd > field.nypf2 - 1 || yEnd < field.nypf1)) {
                            region = 2;
                        }
                    }

                    if (options.direction == 1 && yEnd == field.ny) {
                        region = 14;
                    } else if (options.direction == -1 && yEnd == 1) {
                        region = 13;
                    }

                    zEnd = field.wrap_z(zEnd);
                    zind = field.interp1(field.zarray, field.ziarray, zEnd);

                    ++it;
                    if (it >= nsteps) {
                        region = 15;
                        break;
                    }

                    int ixh = static_cast<int>(std::round(xind));
                    ixh = clamp_int(ixh, 1, field.nx);
                    int iyh = clamp_int(yStart, 1, field.ny);

                    set_traj(traj, 0, it, nsteps, static_cast<float>(iturn));
                    set_traj(traj, 1, it, nsteps, xind);
                    set_traj(traj, 2, it, nsteps, static_cast<float>(yEnd));
                    set_traj(traj, 3, it, nsteps, zind);
                    set_traj(traj, 4, it, nsteps, static_cast<float>(region));
                    set_traj(traj, 5, it, nsteps, field.hthe[field.idx2(ixh - 1, iyh - 1)]);

                    xStart = xEnd;
                    yStart = yEnd;
                    zStart = zEnd;
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

        for (int istep = 0; istep < itmax; ++istep) {
            int yidx = static_cast<int>(std::round(get_traj(traj, 2, istep, nsteps))) - 1;
            yidx = clamp_int(yidx, 0, field.ny - 1);
            float xind_step = get_traj(traj, 1, istep, nsteps);

            float rxyvalue = interp_xindex_2d(field, field.rxy, yidx, xind_step);
            float zsvalue = interp_xindex_2d(field, field.zShift, yidx, xind_step);
            float zxyvalue = interp_xindex_2d(field, field.zxy, yidx, xind_step);
            float zvalue = field.interp1(field.ziarray, field.zarray, get_traj(traj, 3, istep, nsteps));

            float x3d_tmp = rxyvalue * std::cos(zsvalue);
            float y3d_tmp = rxyvalue * std::sin(zsvalue);
            fl_x3d[istep] = x3d_tmp * std::cos(zvalue) - y3d_tmp * std::sin(zvalue);
            fl_y3d[istep] = x3d_tmp * std::sin(zvalue) + y3d_tmp * std::cos(zvalue);
            fl_z3d[istep] = zxyvalue;

            traj_xyz << iline << " " << (istep + 1) << " "
                     << fl_x3d[istep] << " " << fl_y3d[istep] << " " << fl_z3d[istep] << "\n";
        }

        int ip = 0;
        if (itmax > 1) {
            for (int istep = 1; istep < itmax && ip < np_max; ++istep) {
                float x0 = fl_x3d[istep - 1];
                float x1 = fl_x3d[istep];
                float prod = x0 * x1;
                if (prod > 0.0f) {
                    continue;
                }
                float denom = x1 - x0;
                if (std::fabs(denom) < 1.0e-20f) {
                    continue;
                }

                float a = -x0 / denom;
                if (a < 0.0f || a > 1.0f) {
                    continue;
                }
                float b = 1.0f - a;

                int tc0 = istep - 1;
                int tc1 = istep;

                float xind_tmp = b * get_traj(traj, 1, tc0, nsteps) + a * get_traj(traj, 1, tc1, nsteps);
                float yind_tmp = b * get_traj(traj, 2, tc0, nsteps) + a * get_traj(traj, 2, tc1, nsteps);

                float zvalue = b * get_traj(traj, 6, tc0, nsteps) + a * get_traj(traj, 6, tc1, nsteps);
                if (std::fabs(get_traj(traj, 6, tc0, nsteps) - get_traj(traj, 6, tc1, nsteps)) > 1.0f) {
                    float z0w = field.wrap_z(get_traj(traj, 6, tc0, nsteps));
                    float z1w = field.wrap_z(get_traj(traj, 6, tc1, nsteps));
                    zvalue = b * z0w + a * z1w;
                }

                if (static_cast<int>(std::round(get_traj(traj, 2, tc0, nsteps))) == field.nypf2 &&
                    options.direction == 1 &&
                    xind_tmp < static_cast<float>(field.ixsep) + 0.5f) {
                    yind_tmp = b * get_traj(traj, 2, tc0, nsteps) + a * static_cast<float>(field.nypf2 + 1);
                } else if (static_cast<int>(std::round(get_traj(traj, 2, tc0, nsteps))) == (field.nypf1 + 1) &&
                           options.direction == -1 &&
                           xind_tmp < static_cast<float>(field.ixsep) + 0.5f) {
                    yind_tmp = b * static_cast<float>(field.nypf2 + 1) + a * get_traj(traj, 2, tc1, nsteps);
                    float shiftangle = field.interp1(field.xiarray, field.sa, xind_tmp);
                    zvalue = field.wrap_z(zvalue - shiftangle);
                } else if (tc0 > 0) {
                    int y_prev = static_cast<int>(std::round(get_traj(traj, 2, tc0 - 1, nsteps)));
                    if (y_prev == field.nypf2 || y_prev == (field.nypf1 + 1)) {
                        float z0 = field.interp1(field.ziarray, field.zarray, get_traj(traj, 3, tc0, nsteps));
                        float z1v = field.interp1(field.ziarray, field.zarray, get_traj(traj, 3, tc1, nsteps));
                        zvalue = b * z0 + a * z1v;
                    }
                }
                zvalue = field.wrap_z(zvalue);

                float rxyvalue = 0.0f;
                float zxyvalue = 0.0f;
                float zsvalue = 0.0f;
                if (xind_tmp < static_cast<float>(field.ixsep) + 0.5f) {
                    rxyvalue = field.interp2_rect(field.xiarray_cfr, field.yiarray_cfr, field.rxy_cfr,
                                                  field.ixsep, field.ny_cfr, xind_tmp, yind_tmp);
                    zxyvalue = field.interp2_rect(field.xiarray_cfr, field.yiarray_cfr, field.zxy_cfr,
                                                  field.ixsep, field.ny_cfr, xind_tmp, yind_tmp);
                    zsvalue = field.interp2_rect(field.xiarray_cfr, field.yiarray_cfr, field.zs_cfr,
                                                 field.ixsep, field.ny_cfr, xind_tmp, yind_tmp);
                } else {
                    rxyvalue = field.interp2_rect(field.xiarray, field.yiarray, field.rxy,
                                                  field.nx, field.ny, xind_tmp, yind_tmp);
                    zxyvalue = field.interp2_rect(field.xiarray, field.yiarray, field.zxy,
                                                  field.nx, field.ny, xind_tmp, yind_tmp);
                    zsvalue = field.interp2_rect(field.xiarray, field.yiarray, field.zShift,
                                                 field.nx, field.ny, xind_tmp, yind_tmp);
                }

                float ipx3d_tmp = rxyvalue * std::cos(zsvalue);
                float ipy3d_tmp = rxyvalue * std::sin(zsvalue);
                float ipx = ipx3d_tmp * std::cos(zvalue) - ipy3d_tmp * std::sin(zvalue);
                float ipy = ipx3d_tmp * std::sin(zvalue) + ipy3d_tmp * std::cos(zvalue);
                float ipz = zxyvalue;

                if (ipy > 0.0f) {
                    px[ip] = ipx;
                    py[ip] = ipy;
                    pz[ip] = ipz;
                    ip_xyz << iline << " " << istep << " " << ipx << " " << ipy << " " << ipz << "\n";

                    ptheta[ip] = field.interp1(field.yiarray_cfr, field.theta_cfr, yind_tmp);
                    ppsi[ip] = field.interp1(field.xiarray, field.xarray, xind_tmp);
                    float theta_rad = ptheta[ip] * kPi;
                    ip_thetapsi << iline << " " << istep << " " << theta_rad << " " << ppsi[ip] << "\n";
                    ++ip;
                }
            }
        }

        float lc = 0.0f;
        for (int istep = 0; istep < itmax; ++istep) {
            lc += get_traj(traj, 5, istep, nsteps);
        }

        std::cout << "    line " << iline << ": " << ip
                  << " punctures, connection length=" << lc
                  << ", end region=" << region << std::endl;

        (void)yyy;
        (void)zzz;
        (void)px;
        (void)py;
        (void)pz;
        (void)ptheta;
        (void)ppsi;
    }

    std::cout << "Field-line tracing complete." << std::endl;
}
