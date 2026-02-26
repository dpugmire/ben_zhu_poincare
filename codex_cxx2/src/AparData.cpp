#include "AparData.h"

#include <netcdf.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kPi = 3.1415926535897932384626433832795;

void ncCheck(int status, const std::string& where) {
    if (status != NC_NOERR) {
        throw std::runtime_error(where + ": " + nc_strerror(status));
    }
}

int requireVarId(int ncid, const std::string& name) {
    int varid = -1;
    ncCheck(nc_inq_varid(ncid, name.c_str(), &varid), "nc_inq_varid(" + name + ")");
    return varid;
}

size_t requireDimSize(int ncid, const std::string& name) {
    int dimid = -1;
    ncCheck(nc_inq_dimid(ncid, name.c_str(), &dimid), "nc_inq_dimid(" + name + ")");
    size_t dimLen = 0;
    ncCheck(nc_inq_dimlen(ncid, dimid, &dimLen), "nc_inq_dimlen(" + name + ")");
    return dimLen;
}

int readScalarIntVar(int ncid, const std::string& name) {
    const int varid = requireVarId(ncid, name);
    int value = 0;
    ncCheck(nc_get_var_int(ncid, varid, &value), "nc_get_var_int(" + name + ")");
    return value;
}

std::vector<double> readVarDouble(int ncid, const std::string& name) {
    const int varid = requireVarId(ncid, name);

    int ndims = 0;
    ncCheck(nc_inq_varndims(ncid, varid, &ndims), "nc_inq_varndims(" + name + ")");

    std::vector<int> dimids(static_cast<size_t>(ndims), 0);
    ncCheck(nc_inq_vardimid(ncid, varid, dimids.data()), "nc_inq_vardimid(" + name + ")");

    size_t total = 1;
    for (int dimid : dimids) {
        size_t len = 0;
        ncCheck(nc_inq_dimlen(ncid, dimid, &len), "nc_inq_dimlen(" + name + ")");
        total *= len;
    }

    std::vector<double> data(total, 0.0);
    ncCheck(nc_get_var_double(ncid, varid, data.data()), "nc_get_var_double(" + name + ")");
    return data;
}

std::vector<std::string> getVarDimNames(int ncid, const std::string& name) {
    const int varid = requireVarId(ncid, name);

    int ndims = 0;
    ncCheck(nc_inq_varndims(ncid, varid, &ndims), "nc_inq_varndims(" + name + ")");

    std::vector<int> dimids(static_cast<size_t>(ndims), 0);
    ncCheck(nc_inq_vardimid(ncid, varid, dimids.data()), "nc_inq_vardimid(" + name + ")");

    std::vector<std::string> out;
    out.reserve(static_cast<size_t>(ndims));

    for (int dimid : dimids) {
        char dimName[NC_MAX_NAME + 1] = {0};
        ncCheck(nc_inq_dimname(ncid, dimid, dimName), "nc_inq_dimname(" + name + ")");
        out.emplace_back(dimName);
    }

    return out;
}

std::vector<double> convert2DToIxIy(int ncid,
                                    const std::string& varName,
                                    int nx,
                                    int ny) {
    const std::vector<std::string> dims = getVarDimNames(ncid, varName);
    if (dims.size() != 2) {
        throw std::runtime_error(varName + " must be 2D");
    }

    const std::vector<double> raw = readVarDouble(ncid, varName);
    std::vector<double> out(static_cast<size_t>(nx) * ny, 0.0);

    if (dims[0] == "ny" && dims[1] == "nx") {
        for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
                out[static_cast<size_t>(ix) * ny + iy] = raw[static_cast<size_t>(iy) * nx + ix];
            }
        }
        return out;
    }

    if (dims[0] == "nx" && dims[1] == "ny") {
        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                out[static_cast<size_t>(ix) * ny + iy] = raw[static_cast<size_t>(ix) * ny + iy];
            }
        }
        return out;
    }

    // Fallback assumes [ny,nx].
    for (int iy = 0; iy < ny; ++iy) {
        for (int ix = 0; ix < nx; ++ix) {
            out[static_cast<size_t>(ix) * ny + iy] = raw[static_cast<size_t>(iy) * nx + ix];
        }
    }

    return out;
}

std::vector<double> convert2DToIxIyLocal(int ncid,
                                         const std::string& varName,
                                         int nx,
                                         int ny,
                                         const std::string& dimX,
                                         const std::string& dimY) {
    const std::vector<std::string> dims = getVarDimNames(ncid, varName);
    if (dims.size() != 2) {
        throw std::runtime_error(varName + " must be 2D");
    }

    const std::vector<double> raw = readVarDouble(ncid, varName);
    std::vector<double> out(static_cast<size_t>(nx) * ny, 0.0);

    if (dims[0] == dimY && dims[1] == dimX) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
                out[static_cast<size_t>(ix) * ny + iy] = raw[static_cast<size_t>(iy) * nx + ix];
            }
        }
        return out;
    }

    if (dims[0] == dimX && dims[1] == dimY) {
        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                out[static_cast<size_t>(ix) * ny + iy] = raw[static_cast<size_t>(ix) * ny + iy];
            }
        }
        return out;
    }

    for (int iy = 0; iy < ny; ++iy) {
        for (int ix = 0; ix < nx; ++ix) {
            out[static_cast<size_t>(ix) * ny + iy] = raw[static_cast<size_t>(iy) * nx + ix];
        }
    }

    return out;
}

std::vector<double> convert3DToIxIyIzReplicated(int ncid,
                                                const std::string& varName,
                                                int nx,
                                                int ny,
                                                int nz,
                                                int zperiod) {
    const std::vector<std::string> dims = getVarDimNames(ncid, varName);
    if (dims.size() != 3) {
        throw std::runtime_error(varName + " must be 3D");
    }

    const std::vector<double> raw = readVarDouble(ncid, varName);

    const int nzG = nz * zperiod;
    std::vector<double> out(static_cast<size_t>(nx) * ny * nzG, 0.0);

    auto setOut = [&out, ny, nzG](int ix, int iy, int iz, double value) {
        out[(static_cast<size_t>(ix) * ny + iy) * nzG + iz] = value;
    };

    if (dims[0] == "nz" && dims[1] == "ny" && dims[2] == "nx") {
        for (int iz = 0; iz < nz; ++iz) {
            for (int iy = 0; iy < ny; ++iy) {
                for (int ix = 0; ix < nx; ++ix) {
                    const double value = raw[(static_cast<size_t>(iz) * ny + iy) * nx + ix];
                    for (int zp = 0; zp < zperiod; ++zp) {
                        setOut(ix, iy, zp * nz + iz, value);
                    }
                }
            }
        }
        return out;
    }

    if (dims[0] == "nx" && dims[1] == "ny" && dims[2] == "nz") {
        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                for (int iz = 0; iz < nz; ++iz) {
                    const double value = raw[(static_cast<size_t>(ix) * ny + iy) * nz + iz];
                    for (int zp = 0; zp < zperiod; ++zp) {
                        setOut(ix, iy, zp * nz + iz, value);
                    }
                }
            }
        }
        return out;
    }

    for (int iz = 0; iz < nz; ++iz) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
                const double value = raw[(static_cast<size_t>(iz) * ny + iy) * nx + ix];
                for (int zp = 0; zp < zperiod; ++zp) {
                    setOut(ix, iy, zp * nz + iz, value);
                }
            }
        }
    }

    return out;
}

std::vector<double> convert2DXZReplicated(int ncid,
                                          const std::string& varName,
                                          int nx,
                                          int nz,
                                          int zperiod) {
    const std::vector<std::string> dims = getVarDimNames(ncid, varName);
    if (dims.size() != 2) {
        throw std::runtime_error(varName + " must be 2D");
    }

    const std::vector<double> raw = readVarDouble(ncid, varName);

    const int nzG = nz * zperiod;
    std::vector<double> out(static_cast<size_t>(nx) * nzG, 0.0);

    auto setOut = [&out, nzG](int ix, int iz, double value) {
        out[static_cast<size_t>(ix) * nzG + iz] = value;
    };

    if (dims[0] == "nz" && dims[1] == "nx") {
        for (int iz = 0; iz < nz; ++iz) {
            for (int ix = 0; ix < nx; ++ix) {
                const double value = raw[static_cast<size_t>(iz) * nx + ix];
                for (int zp = 0; zp < zperiod; ++zp) {
                    setOut(ix, zp * nz + iz, value);
                }
            }
        }
        return out;
    }

    if (dims[0] == "nx" && dims[1] == "nz") {
        for (int ix = 0; ix < nx; ++ix) {
            for (int iz = 0; iz < nz; ++iz) {
                const double value = raw[static_cast<size_t>(ix) * nz + iz];
                for (int zp = 0; zp < zperiod; ++zp) {
                    setOut(ix, zp * nz + iz, value);
                }
            }
        }
        return out;
    }

    for (int iz = 0; iz < nz; ++iz) {
        for (int ix = 0; ix < nx; ++ix) {
            const double value = raw[static_cast<size_t>(iz) * nx + ix];
            for (int zp = 0; zp < zperiod; ++zp) {
                setOut(ix, zp * nz + iz, value);
            }
        }
    }

    return out;
}

}  // namespace

void AparData::load(const std::string& aparPath) {
    int ncid = -1;
    ncCheck(nc_open(aparPath.c_str(), NC_NOWRITE, &ncid), "nc_open(" + aparPath + ")");

    try {
        nx = static_cast<int>(requireDimSize(ncid, "nx"));
        ny = static_cast<int>(requireDimSize(ncid, "ny"));
        nz = static_cast<int>(requireDimSize(ncid, "nz"));
        nx_cfr = static_cast<int>(requireDimSize(ncid, "nx_cfr"));
        ny_cfr = static_cast<int>(requireDimSize(ncid, "ny_cfr"));

        zperiod = readScalarIntVar(ncid, "zperiod");

        if (nx <= 0 || ny <= 0 || nz <= 0 || zperiod <= 0) {
            throw std::runtime_error("Invalid dimensions in apar NetCDF file");
        }

        nzG = nz * zperiod;

        ixsep1 = readScalarIntVar(ncid, "ixsep1");
        ixsep2 = readScalarIntVar(ncid, "ixsep2");
        nypf1 = readScalarIntVar(ncid, "nypf1");
        nypf2 = readScalarIntVar(ncid, "nypf2");

        psixy = convert2DToIxIy(ncid, "psixy", nx, ny);
        dxdy = convert3DToIxIyIzReplicated(ncid, "dxdy", nx, ny, nz, zperiod);
        dzdy = convert3DToIxIyIzReplicated(ncid, "dzdy", nx, ny, nz, zperiod);
        dxdy_p1 = convert2DXZReplicated(ncid, "dxdy_p1", nx, nz, zperiod);
        dzdy_p1 = convert2DXZReplicated(ncid, "dzdy_p1", nx, nz, zperiod);
        dxdy_m1 = convert2DXZReplicated(ncid, "dxdy_m1", nx, nz, zperiod);
        dzdy_m1 = convert2DXZReplicated(ncid, "dzdy_m1", nx, nz, zperiod);

        shiftAngle = readVarDouble(ncid, "shiftAngle");
        if (static_cast<int>(shiftAngle.size()) != nx) {
            throw std::runtime_error("shiftAngle size mismatch");
        }

        zShift = convert2DToIxIy(ncid, "zShift", nx, ny);
        rxy = convert2DToIxIy(ncid, "rxy", nx, ny);
        zxy = convert2DToIxIy(ncid, "zxy", nx, ny);

        rxy_cfr = convert2DToIxIyLocal(ncid, "rxy_cfr", nx_cfr, ny_cfr, "nx_cfr", "ny_cfr");
        zxy_cfr = convert2DToIxIyLocal(ncid, "zxy_cfr", nx_cfr, ny_cfr, "nx_cfr", "ny_cfr");
        zShift_cfr = convert2DToIxIyLocal(ncid, "zShift_cfr", nx_cfr, ny_cfr, "nx_cfr", "ny_cfr");
    } catch (...) {
        nc_close(ncid);
        throw;
    }

    ncCheck(nc_close(ncid), "nc_close(" + aparPath + ")");

    computeDerivedGeometry();
    computeThetaProfiles();
}

void AparData::computeDerivedGeometry() {
    if (ixsep2 < nx) {
        divertor = 2;
        ixsep = ixsep1;
    } else if (ixsep1 < nx) {
        divertor = 1;
        ixsep = ixsep1;
    } else {
        divertor = 0;
        ixsep = nx;
        nypf1 = 0;
        nypf2 = ny;
    }

    zmin = 0.0;
    zmax = 2.0 * kPi;
    dz_torus = (zmax - zmin) / static_cast<double>(nzG);

    xiarray.resize(nx);
    yiarray.resize(ny);
    ziarray.resize(nzG + 1);
    zarray.resize(nzG + 1);

    for (int ix = 0; ix < nx; ++ix) {
        xiarray[ix] = static_cast<double>(ix + 1);
    }
    for (int iy = 0; iy < ny; ++iy) {
        yiarray[iy] = static_cast<double>(iy + 1);
    }
    for (int i = 0; i <= nzG; ++i) {
        ziarray[i] = static_cast<double>(i + 1);
        zarray[i] = static_cast<double>(i) * dz_torus;
    }

    jyomp = 0;
    double rmax = rxy[idx2(nx - 1, 0)];
    for (int iy = 1; iy < ny; ++iy) {
        const double value = rxy[idx2(nx - 1, iy)];
        if (value > rmax) {
            rmax = value;
            jyomp = iy;
        }
    }

    xarray.resize(nx);
    for (int ix = 0; ix < nx; ++ix) {
        xarray[ix] = psixy[idx2(ix, jyomp)];
    }

    const auto minmaxX = std::minmax_element(xarray.begin(), xarray.end());
    xMin = *minmaxX.first;
    xMax = *minmaxX.second;

    xiarray_cfr.resize(nx_cfr);
    for (int ix = 0; ix < nx_cfr; ++ix) {
        xiarray_cfr[ix] = static_cast<double>(ix + 1);
    }

    yiarray_cfr.resize(ny_cfr);
    if (divertor == 0) {
        for (int iy = 0; iy < ny_cfr; ++iy) {
            yiarray_cfr[iy] = static_cast<double>(iy + 1);
        }
    } else {
        for (int iy = 0; iy < ny_cfr; ++iy) {
            yiarray_cfr[iy] = static_cast<double>(nypf1 + 1 + iy);
        }
    }
}

void AparData::computeThetaProfiles() {
    theta.assign(ny, 0.0);

    if (divertor == 1) {
        double coreRmin = rxy[idx2(0, nypf1)];
        double coreRmax = coreRmin;
        double coreZmin = zxy[idx2(0, nypf1)];
        double coreZmax = coreZmin;

        for (int iy = nypf1; iy < ny - nypf1; ++iy) {
            const double rv = rxy[idx2(0, iy)];
            const double zv = zxy[idx2(0, iy)];
            coreRmin = std::min(coreRmin, rv);
            coreRmax = std::max(coreRmax, rv);
            coreZmin = std::min(coreZmin, zv);
            coreZmax = std::max(coreZmax, zv);
        }

        const double centerX = 0.5 * (coreRmax + coreRmin);
        const double centerY = 0.5 * (coreZmax + coreZmin);

        std::vector<int> tmp;
        tmp.reserve(ny);
        for (int iy = 0; iy < nypf1; ++iy) {
            tmp.push_back(iy);
        }
        for (int iy = ny - nypf1 - 1; iy >= nypf1; --iy) {
            tmp.push_back(iy);
        }
        for (int iy = ny - nypf1; iy < ny; ++iy) {
            tmp.push_back(iy);
        }

        std::vector<double> sepx(tmp.size(), 0.0);
        std::vector<double> sepy(tmp.size(), 0.0);
        for (size_t i = 0; i < tmp.size(); ++i) {
            const int iy = tmp[i];
            sepx[i] = 0.5 * (rxy[idx2(ixsep - 1, iy)] + rxy[idx2(ixsep, iy)]);
            sepy[i] = 0.5 * (zxy[idx2(ixsep - 1, iy)] + zxy[idx2(ixsep, iy)]);
        }

        const double xpointX = 0.25 * (sepx[nypf1 - 1] + sepx[nypf1] +
                                       sepx[ny - nypf1 - 1] + sepx[ny - nypf1]);
        const double xpointY = 0.25 * (sepy[nypf1 - 1] + sepy[nypf1] +
                                       sepy[ny - nypf1 - 1] + sepy[ny - nypf1]);

        const double ux = centerX - xpointX;
        const double uy = centerY - xpointY;

        for (int iy = 0; iy < ny; ++iy) {
            const double vx = centerX - rxy[idx2(0, iy)];
            const double vy = centerY - zxy[idx2(0, iy)];
            const double cross = std::fabs(ux * vy - uy * vx);
            const double dot = ux * vx + uy * vy;
            theta[iy] = std::atan2(cross, dot) / kPi;
        }

        int itheta = 0;
        double thmax = theta[0];
        for (int iy = 1; iy < ny; ++iy) {
            if (theta[iy] > thmax) {
                thmax = theta[iy];
                itheta = iy;
            }
        }
        for (int iy = itheta; iy < ny; ++iy) {
            theta[iy] = 2.0 - theta[iy];
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
                theta[iy] = 4.0 - theta[iy];
            }
        }

        const double ref = theta[nypf1];
        for (double& v : theta) {
            v -= ref;
        }
    } else {
        double rmin = rxy[idx2(0, 0)];
        double rmax = rmin;
        double zminLocal = zxy[idx2(0, 0)];
        double zmaxLocal = zminLocal;

        for (int iy = 0; iy < ny; ++iy) {
            const double rv = rxy[idx2(0, iy)];
            const double zv = zxy[idx2(0, iy)];
            rmin = std::min(rmin, rv);
            rmax = std::max(rmax, rv);
            zminLocal = std::min(zminLocal, zv);
            zmaxLocal = std::max(zmaxLocal, zv);
        }

        const double centerX = 0.5 * (rmax + rmin);
        const double centerY = 0.5 * (zmaxLocal + zminLocal);
        const double ux = centerX - rxy[idx2(0, 0)];
        const double uy = centerY - zxy[idx2(0, 0)];

        for (int iy = 0; iy < ny; ++iy) {
            const double vx = centerX - rxy[idx2(0, iy)];
            const double vy = centerY - zxy[idx2(0, iy)];
            const double cross = std::fabs(ux * vy - uy * vx);
            const double dot = ux * vx + uy * vy;
            theta[iy] = std::atan2(cross, dot) / kPi;
        }

        int itheta = 0;
        double thmax = theta[0];
        for (int iy = 1; iy < ny; ++iy) {
            if (theta[iy] > thmax) {
                thmax = theta[iy];
                itheta = iy;
            }
        }
        for (int iy = itheta; iy < ny; ++iy) {
            theta[iy] = 2.0 - theta[iy];
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
                theta[iy] = 4.0 - theta[iy];
            }
        }
    }

    theta_cfr.assign(ny_cfr, 0.0);
    if (divertor == 0) {
        for (int iy = 0; iy < ny_cfr - 1 && iy < ny; ++iy) {
            theta_cfr[iy] = theta[iy];
        }
    } else {
        for (int iy = 0; iy < ny_cfr; ++iy) {
            const int src = nypf1 + iy;
            if (src >= 0 && src < ny) {
                theta_cfr[iy] = theta[src];
            }
        }
    }

    if (!theta_cfr.empty()) {
        theta_cfr.back() = 2.0;
    }
}

double AparData::wrapZ(double z) const {
    const double period = zmax;
    if (period <= 0.0) {
        return z;
    }
    double out = std::fmod(z, period);
    if (out < 0.0) {
        out += period;
    }
    if (out >= period) {
        out -= period;
    }
    return out;
}
