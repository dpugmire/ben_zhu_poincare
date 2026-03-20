#ifndef CODEX_CXX2_VISKORESTRACEWORKLET_H
#define CODEX_CXX2_VISKORESTRACEWORKLET_H

#include "AparData.h"
#include "Types.h"

#if defined(CODEX_USE_VISKORES)

#include <cmath>

#include <viskores/cont/ArrayHandle.h>
#include <viskores/cont/ExecutionObjectBase.h>
#include <viskores/worklet/WorkletMapField.h>

class ViskoresAparFieldExecutionObject
{
public:
  using FloatType = viskores::FloatDefault;
  using FloatArray = viskores::cont::ArrayHandle<FloatType>;
  using FloatPortal = typename FloatArray::ReadPortalType;

  ViskoresAparFieldExecutionObject(const FloatPortal& xiarrayIn,
                                   const FloatPortal& xarrayIn,
                                   const FloatPortal& ziarrayIn,
                                   const FloatPortal& zarrayIn,
                                   const FloatPortal& dxdyIn,
                                   const FloatPortal& dzdyIn,
                                   const FloatPortal& dxdyP1In,
                                   const FloatPortal& dzdyP1In,
                                   const FloatPortal& dxdyM1In,
                                   const FloatPortal& dzdyM1In,
                                   const FloatPortal& shiftAngleIn,
                                   int nxIn,
                                   int nyIn,
                                   int nzGIn,
                                   int ixsepIn,
                                   int nypf1In,
                                   int nypf2In,
                                   FloatType zmaxIn,
                                   FloatType dzTorusIn,
                                   FloatType xMinIn,
                                   FloatType xMaxIn)
    : xiarray(xiarrayIn)
    , xarray(xarrayIn)
    , ziarray(ziarrayIn)
    , zarray(zarrayIn)
    , dxdy(dxdyIn)
    , dzdy(dzdyIn)
    , dxdy_p1(dxdyP1In)
    , dzdy_p1(dzdyP1In)
    , dxdy_m1(dxdyM1In)
    , dzdy_m1(dzdyM1In)
    , shiftAngle(shiftAngleIn)
    , nx(nxIn)
    , ny(nyIn)
    , nzG(nzGIn)
    , ixsep(ixsepIn)
    , nypf1(nypf1In)
    , nypf2(nypf2In)
    , zmax(zmaxIn)
    , dz_torus(dzTorusIn)
    , xMin(xMinIn)
    , xMax(xMaxIn)
  {
  }

  VISKORES_EXEC int clampInt(int value, int lo, int hi) const
  {
    if (value < lo)
    {
      return lo;
    }
    if (value > hi)
    {
      return hi;
    }
    return value;
  }

  template <typename PortalType>
  VISKORES_EXEC int lowerBracket(const PortalType& xp, FloatType x) const
  {
    const int n = static_cast<int>(xp.GetNumberOfValues());
    if (n < 2)
    {
      return 0;
    }
    if (x <= xp.Get(0))
    {
      return 0;
    }
    if (x >= xp.Get(n - 1))
    {
      return n - 2;
    }

    int lo = 0;
    int hi = n - 1;
    while (hi - lo > 1)
    {
      const int mid = lo + (hi - lo) / 2;
      if (xp.Get(mid) <= x)
      {
        lo = mid;
      }
      else
      {
        hi = mid;
      }
    }
    return lo;
  }

  template <typename XPortalType, typename FPortalType>
  VISKORES_EXEC FloatType interp1(const XPortalType& xp, const FPortalType& fp, FloatType x) const
  {
    const int n = static_cast<int>(xp.GetNumberOfValues());
    if (n <= 0 || fp.GetNumberOfValues() != xp.GetNumberOfValues())
    {
      return 0.0;
    }
    if (n == 1)
    {
      return fp.Get(0);
    }
    if (x <= xp.Get(0))
    {
      return fp.Get(0);
    }
    if (x >= xp.Get(n - 1))
    {
      return fp.Get(n - 1);
    }

    const int i0 = lowerBracket(xp, x);
    const int i1 = i0 + 1;
    const FloatType denom = xp.Get(i1) - xp.Get(i0);
    if (std::fabs(denom) < 1.0e-20)
    {
      return fp.Get(i0);
    }
    const FloatType t = (x - xp.Get(i0)) / denom;
    return fp.Get(i0) + t * (fp.Get(i1) - fp.Get(i0));
  }

  VISKORES_EXEC FloatType wrapZ(FloatType z) const
  {
    if (zmax <= 0.0)
    {
      return z;
    }
    FloatType out = std::fmod(z, zmax);
    if (out < 0.0)
    {
      out += zmax;
    }
    if (out >= zmax)
    {
      out -= zmax;
    }
    return out;
  }

  VISKORES_EXEC viskores::Id idx3(viskores::Id ix, viskores::Id iy, viskores::Id iz) const
  {
    return (static_cast<viskores::Id>(ix) * static_cast<viskores::Id>(ny) + static_cast<viskores::Id>(iy)) * static_cast<viskores::Id>(nzG) +
      static_cast<viskores::Id>(iz);
  }

  VISKORES_EXEC viskores::Id idxXZ(viskores::Id ix, viskores::Id iz) const
  {
    return static_cast<viskores::Id>(ix) * static_cast<viskores::Id>(nzG) + static_cast<viskores::Id>(iz);
  }

  VISKORES_EXEC void cubicStencilCentered(viskores::Id baseIdx, viskores::Id nNodes, viskores::Id outIdx[4]) const
  {
    if (nNodes < 4)
    {
      outIdx[0] = 0;
      outIdx[1] = clampInt(baseIdx, 0, (nNodes > 0) ? (nNodes - 1) : 0);
      outIdx[2] = outIdx[1];
      outIdx[3] = outIdx[1];
      return;
    }

    if (baseIdx <= 1)
    {
      outIdx[0] = 0;
      outIdx[1] = 1;
      outIdx[2] = 2;
      outIdx[3] = 3;
      return;
    }

    if (baseIdx >= nNodes - 3)
    {
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

  VISKORES_EXEC FloatType
  lagrange4(FloatType x, FloatType x0, FloatType x1, FloatType x2, FloatType x3, FloatType y0, FloatType y1, FloatType y2, FloatType y3) const
  {
    const FloatType xs[4] = { x0, x1, x2, x3 };
    const FloatType ys[4] = { y0, y1, y2, y3 };

    FloatType out = 0.0;
    for (int i = 0; i < 4; ++i)
    {
      FloatType li = 1.0;
      for (int j = 0; j < 4; ++j)
      {
        if (i == j)
        {
          continue;
        }
        const FloatType den = xs[i] - xs[j];
        if (std::fabs(den) < 1.0e-20)
        {
          return ys[1];
        }
        li *= (x - xs[j]) / den;
      }
      out += ys[i] * li;
    }
    return out;
  }

  VISKORES_EXEC FloatType sampleClampedRow3D(const FloatPortal& data3d, int ix, int iy, int izExt) const
  {
    if (izExt <= 0)
    {
      return data3d.Get(idx3(ix, iy, 0));
    }
    if (izExt >= nzG)
    {
      return data3d.Get(idx3(ix, iy, nzG - 1));
    }
    return data3d.Get(idx3(ix, iy, izExt));
  }

  VISKORES_EXEC FloatType sampleClampedXZ(const FloatPortal& data2d, int ix, int izExt) const
  {
    if (izExt <= 0)
    {
      return data2d.Get(idxXZ(ix, 0));
    }
    if (izExt >= nzG)
    {
      return data2d.Get(idxXZ(ix, nzG - 1));
    }
    return data2d.Get(idxXZ(ix, izExt));
  }

  VISKORES_EXEC FloatType interpXZ3DAtY(const FloatPortal& data3d, int y0, FloatType x, FloatType z) const
  {
    if (nx < 2 || nzG < 1)
    {
      return 0.0;
    }

    y0 = clampInt(y0, 0, ny - 1);

    int ix0 = 0;
    FloatType tx = 0.0;
    if (x <= xarray.Get(0))
    {
      ix0 = 0;
      tx = 0.0;
    }
    else if (x >= xarray.Get(nx - 1))
    {
      ix0 = nx - 2;
      tx = 1.0;
    }
    else
    {
      ix0 = lowerBracket(xarray, x);
      const FloatType denom = xarray.Get(ix0 + 1) - xarray.Get(ix0);
      tx = (std::fabs(denom) > 1.0e-20) ? ((x - xarray.Get(ix0)) / denom) : 0.0;
    }
    if (tx < 0.0)
    {
      tx = 0.0;
    }
    if (tx > 1.0)
    {
      tx = 1.0;
    }

    const FloatType zWrapped = wrapZ(z);
    const FloatType zf = zWrapped / dz_torus;
    int iz0 = clampInt(static_cast<int>(std::floor(zf)), 0, nzG - 1);
    const int iz1 = (iz0 + 1) % nzG;
    FloatType tz = zf - static_cast<FloatType>(iz0);
    if (tz < 0.0)
    {
      tz = 0.0;
    }
    if (tz > 1.0)
    {
      tz = 1.0;
    }

    int ix1 = ix0 + 1;
    if (ix1 >= nx)
    {
      ix1 = nx - 1;
    }

    const FloatType v00 = data3d.Get(idx3(ix0, y0, iz0));
    const FloatType v01 = data3d.Get(idx3(ix0, y0, iz1));
    const FloatType v10 = data3d.Get(idx3(ix1, y0, iz0));
    const FloatType v11 = data3d.Get(idx3(ix1, y0, iz1));

    const FloatType v0 = v00 + tz * (v01 - v00);
    const FloatType v1 = v10 + tz * (v11 - v10);
    return v0 + tx * (v1 - v0);
  }

  VISKORES_EXEC FloatType interpXZ3DAtYSpline(const FloatPortal& data3d, int y0, FloatType x, FloatType z) const
  {
    if (nx < 4 || nzG < 4)
    {
      return interpXZ3DAtY(data3d, y0, x, z);
    }

    y0 = clampInt(y0, 0, ny - 1);

    FloatType xq = x;
    if (xq < xarray.Get(0))
    {
      xq = xarray.Get(0);
    }
    if (xq > xarray.Get(nx - 1))
    {
      xq = xarray.Get(nx - 1);
    }

    const int ixBase = lowerBracket(xarray, xq);
    viskores::Id ixs[4] = { 0, 1, 2, 3 };
    cubicStencilCentered(ixBase, nx, ixs);

    const FloatType zq = wrapZ(z);
    int izBase = clampInt(static_cast<int>(std::floor(zq / dz_torus)), 0, nzG - 1);
    viskores::Id izs[4] = { 0, 1, 2, 3 };
    cubicStencilCentered(izBase, nzG + 1, izs);

    const FloatType xvals[4] = { xarray.Get(ixs[0]), xarray.Get(ixs[1]), xarray.Get(ixs[2]), xarray.Get(ixs[3]) };

    FloatType zinterp[4] = { 0.0, 0.0, 0.0, 0.0 };
    for (int i = 0; i < 4; ++i)
    {
      const int ix = ixs[i];

      const FloatType zc0 = static_cast<FloatType>(izs[0]) * dz_torus;
      const FloatType zc1 = static_cast<FloatType>(izs[1]) * dz_torus;
      const FloatType zc2 = static_cast<FloatType>(izs[2]) * dz_torus;
      const FloatType zc3 = static_cast<FloatType>(izs[3]) * dz_torus;

      const FloatType zv0 = sampleClampedRow3D(data3d, ix, y0, izs[0]);
      const FloatType zv1 = sampleClampedRow3D(data3d, ix, y0, izs[1]);
      const FloatType zv2 = sampleClampedRow3D(data3d, ix, y0, izs[2]);
      const FloatType zv3 = sampleClampedRow3D(data3d, ix, y0, izs[3]);

      zinterp[i] = lagrange4(zq, zc0, zc1, zc2, zc3, zv0, zv1, zv2, zv3);
    }

    return lagrange4(xq, xvals[0], xvals[1], xvals[2], xvals[3], zinterp[0], zinterp[1], zinterp[2], zinterp[3]);
  }

  VISKORES_EXEC FloatType interpXZ2D(const FloatPortal& data2d, FloatType x, FloatType z) const
  {
    if (nx < 2 || nzG < 1)
    {
      return 0.0;
    }

    int ix0 = 0;
    FloatType tx = 0.0;
    if (x <= xarray.Get(0))
    {
      ix0 = 0;
      tx = 0.0;
    }
    else if (x >= xarray.Get(nx - 1))
    {
      ix0 = nx - 2;
      tx = 1.0;
    }
    else
    {
      ix0 = lowerBracket(xarray, x);
      const FloatType denom = xarray.Get(ix0 + 1) - xarray.Get(ix0);
      tx = (std::fabs(denom) > 1.0e-20) ? ((x - xarray.Get(ix0)) / denom) : 0.0;
    }
    if (tx < 0.0)
    {
      tx = 0.0;
    }
    if (tx > 1.0)
    {
      tx = 1.0;
    }

    const FloatType zWrapped = wrapZ(z);
    const FloatType zf = zWrapped / dz_torus;
    int iz0 = clampInt(static_cast<int>(std::floor(zf)), 0, nzG - 1);
    const int iz1 = (iz0 + 1) % nzG;
    FloatType tz = zf - static_cast<FloatType>(iz0);
    if (tz < 0.0)
    {
      tz = 0.0;
    }
    if (tz > 1.0)
    {
      tz = 1.0;
    }

    int ix1 = ix0 + 1;
    if (ix1 >= nx)
    {
      ix1 = nx - 1;
    }

    const FloatType v00 = data2d.Get(idxXZ(ix0, iz0));
    const FloatType v01 = data2d.Get(idxXZ(ix0, iz1));
    const FloatType v10 = data2d.Get(idxXZ(ix1, iz0));
    const FloatType v11 = data2d.Get(idxXZ(ix1, iz1));

    const FloatType v0 = v00 + tz * (v01 - v00);
    const FloatType v1 = v10 + tz * (v11 - v10);
    return v0 + tx * (v1 - v0);
  }

  VISKORES_EXEC FloatType interpXZ2DSpline(const FloatPortal& data2d, FloatType x, FloatType z) const
  {
    if (nx < 4 || nzG < 4)
    {
      return interpXZ2D(data2d, x, z);
    }

    FloatType xq = x;
    if (xq < xarray.Get(0))
    {
      xq = xarray.Get(0);
    }
    if (xq > xarray.Get(nx - 1))
    {
      xq = xarray.Get(nx - 1);
    }

    const int ixBase = lowerBracket(xarray, xq);
    viskores::Id ixs[4] = { 0, 1, 2, 3 };
    cubicStencilCentered(ixBase, nx, ixs);

    const FloatType zq = wrapZ(z);
    int izBase = clampInt(static_cast<int>(std::floor(zq / dz_torus)), 0, nzG - 1);
    viskores::Id izs[4] = { 0, 1, 2, 3 };
    cubicStencilCentered(izBase, nzG + 1, izs);

    const FloatType xvals[4] = { xarray.Get(ixs[0]), xarray.Get(ixs[1]), xarray.Get(ixs[2]), xarray.Get(ixs[3]) };

    FloatType zinterp[4] = { 0.0, 0.0, 0.0, 0.0 };
    for (int i = 0; i < 4; ++i)
    {
      const int ix = ixs[i];

      const FloatType zc0 = static_cast<FloatType>(izs[0]) * dz_torus;
      const FloatType zc1 = static_cast<FloatType>(izs[1]) * dz_torus;
      const FloatType zc2 = static_cast<FloatType>(izs[2]) * dz_torus;
      const FloatType zc3 = static_cast<FloatType>(izs[3]) * dz_torus;

      const FloatType zv0 = sampleClampedXZ(data2d, ix, izs[0]);
      const FloatType zv1 = sampleClampedXZ(data2d, ix, izs[1]);
      const FloatType zv2 = sampleClampedXZ(data2d, ix, izs[2]);
      const FloatType zv3 = sampleClampedXZ(data2d, ix, izs[3]);

      zinterp[i] = lagrange4(zq, zc0, zc1, zc2, zc3, zv0, zv1, zv2, zv3);
    }

    return lagrange4(xq, xvals[0], xvals[1], xvals[2], xvals[3], zinterp[0], zinterp[1], zinterp[2], zinterp[3]);
  }

  VISKORES_EXEC void evaluateStage(const XZPoint& point, int yStart1b, int region, int direction, int stage, XZDeriv& deriv) const
  {
    stage = clampInt(stage, 0, 2);

    int yp = clampInt(yStart1b - 1, 0, ny - 1);

    bool useTwist = false;
    bool usePlus = false;
    int yn = yp;

    if (direction == 1)
    {
      if (region == 0 && yStart1b == nypf2)
      {
        useTwist = true;
        usePlus = true;
      }
      else
      {
        yn = clampInt(yStart1b, 0, ny - 1);
      }
    }
    else if (direction == -1)
    {
      if (region == 0 && yStart1b == (nypf1 + 1))
      {
        useTwist = true;
        usePlus = false;
      }
      else
      {
        yn = clampInt(yStart1b - 2, 0, ny - 1);
      }
    }
    else
    {
      deriv.dxdy = 0.0;
      deriv.dzdy = 0.0;
      return;
    }

    if (stage == 0)
    {
      deriv.dxdy = interpXZ3DAtYSpline(dxdy, yp, point.x, point.z);
      deriv.dzdy = interpXZ3DAtYSpline(dzdy, yp, point.x, point.z);
      return;
    }

    if (stage == 2)
    {
      if (useTwist)
      {
        if (usePlus)
        {
          deriv.dxdy = interpXZ2DSpline(dxdy_p1, point.x, point.z);
          deriv.dzdy = interpXZ2DSpline(dzdy_p1, point.x, point.z);
        }
        else
        {
          deriv.dxdy = interpXZ2DSpline(dxdy_m1, point.x, point.z);
          deriv.dzdy = interpXZ2DSpline(dzdy_m1, point.x, point.z);
        }
      }
      else
      {
        deriv.dxdy = interpXZ3DAtYSpline(dxdy, yn, point.x, point.z);
        deriv.dzdy = interpXZ3DAtYSpline(dzdy, yn, point.x, point.z);
      }
      return;
    }

    const FloatType dxP = interpXZ3DAtYSpline(dxdy, yp, point.x, point.z);
    const FloatType dzP = interpXZ3DAtYSpline(dzdy, yp, point.x, point.z);

    FloatType dxN = 0.0;
    FloatType dzN = 0.0;
    if (useTwist)
    {
      if (usePlus)
      {
        dxN = interpXZ2DSpline(dxdy_p1, point.x, point.z);
        dzN = interpXZ2DSpline(dzdy_p1, point.x, point.z);
      }
      else
      {
        dxN = interpXZ2DSpline(dxdy_m1, point.x, point.z);
        dzN = interpXZ2DSpline(dzdy_m1, point.x, point.z);
      }
    }
    else
    {
      dxN = interpXZ3DAtYSpline(dxdy, yn, point.x, point.z);
      dzN = interpXZ3DAtYSpline(dzdy, yn, point.x, point.z);
    }

    deriv.dxdy = 0.5 * (dxP + dxN);
    deriv.dzdy = 0.5 * (dzP + dzN);
  }

  FloatPortal xiarray;
  FloatPortal xarray;
  FloatPortal ziarray;
  FloatPortal zarray;
  FloatPortal dxdy;
  FloatPortal dzdy;
  FloatPortal dxdy_p1;
  FloatPortal dzdy_p1;
  FloatPortal dxdy_m1;
  FloatPortal dzdy_m1;
  FloatPortal shiftAngle;

  int nx = 0;
  int ny = 0;
  int nzG = 0;
  int ixsep = 0;
  int nypf1 = 0;
  int nypf2 = 0;
  FloatType zmax = 0.0;
  FloatType dz_torus = 0.0;
  FloatType xMin = 0.0;
  FloatType xMax = 0.0;
};

class ViskoresAparField : public viskores::cont::ExecutionObjectBase
{
public:
  using FloatType = viskores::FloatDefault;
  using FloatArray = viskores::cont::ArrayHandle<FloatType>;

  explicit ViskoresAparField(const AparData& data)
    : xiarray_(makeFloatArray(data.xiarray))
    , xarray_(makeFloatArray(data.xarray))
    , ziarray_(makeFloatArray(data.ziarray))
    , zarray_(makeFloatArray(data.zarray))
    , dxdy_(makeFloatArray(data.dxdy))
    , dzdy_(makeFloatArray(data.dzdy))
    , dxdyP1_(makeFloatArray(data.dxdy_p1))
    , dzdyP1_(makeFloatArray(data.dzdy_p1))
    , dxdyM1_(makeFloatArray(data.dxdy_m1))
    , dzdyM1_(makeFloatArray(data.dzdy_m1))
    , shiftAngle_(makeFloatArray(data.shiftAngle))
    , nx_(data.nx)
    , ny_(data.ny)
    , nzG_(data.nzG)
    , ixsep_(data.ixsep)
    , nypf1_(data.nypf1)
    , nypf2_(data.nypf2)
    , zmax_(data.zmax)
    , dzTorus_(data.dz_torus)
    , xMin_(data.xMin)
    , xMax_(data.xMax)
  {
  }

  VISKORES_CONT ViskoresAparFieldExecutionObject PrepareForExecution(viskores::cont::DeviceAdapterId device, viskores::cont::Token& token) const
  {
    return ViskoresAparFieldExecutionObject(xiarray_.PrepareForInput(device, token),
                                            xarray_.PrepareForInput(device, token),
                                            ziarray_.PrepareForInput(device, token),
                                            zarray_.PrepareForInput(device, token),
                                            dxdy_.PrepareForInput(device, token),
                                            dzdy_.PrepareForInput(device, token),
                                            dxdyP1_.PrepareForInput(device, token),
                                            dzdyP1_.PrepareForInput(device, token),
                                            dxdyM1_.PrepareForInput(device, token),
                                            dzdyM1_.PrepareForInput(device, token),
                                            shiftAngle_.PrepareForInput(device, token),
                                            nx_,
                                            ny_,
                                            nzG_,
                                            ixsep_,
                                            nypf1_,
                                            nypf2_,
                                            zmax_,
                                            dzTorus_,
                                            xMin_,
                                            xMax_);
  }

private:
  template <typename ValueType>
  static FloatArray makeFloatArray(const std::vector<ValueType>& values)
  {
    FloatArray out;
    out.Allocate(static_cast<viskores::Id>(values.size()));
    auto portal = out.WritePortal();
    for (std::size_t i = 0; i < values.size(); ++i)
    {
      portal.Set(static_cast<viskores::Id>(i), static_cast<FloatType>(values[i]));
    }
    return out;
  }

  FloatArray xiarray_;
  FloatArray xarray_;
  FloatArray ziarray_;
  FloatArray zarray_;
  FloatArray dxdy_;
  FloatArray dzdy_;
  FloatArray dxdyP1_;
  FloatArray dzdyP1_;
  FloatArray dxdyM1_;
  FloatArray dzdyM1_;
  FloatArray shiftAngle_;

  int nx_ = 0;
  int ny_ = 0;
  int nzG_ = 0;
  int ixsep_ = 0;
  int nypf1_ = 0;
  int nypf2_ = 0;
  FloatType zmax_ = 0.0;
  FloatType dzTorus_ = 0.0;
  FloatType xMin_ = 0.0;
  FloatType xMax_ = 0.0;
};

class ViskoresTraceStatesWorklet : public viskores::worklet::WorkletMapField
{
public:
  using FloatType = viskores::FloatDefault;

  ViskoresTraceStatesWorklet(viskores::Id maxStatesPerSeed, int direction)
    : maxStatesPerSeed_(maxStatesPerSeed)
    , direction_(direction)
  {
  }

  using ControlSignature =
    void(FieldIn seeds, ExecObject field, WholeArrayOut stateCounts, WholeArrayOut endRegions, WholeArrayOut statusCodes, WholeArrayOut states);
  using ExecutionSignature = void(InputIndex, _1, _2, _3, _4, _5, _6);
  using InputDomain = _1;

  template <typename FieldExecType, typename IntPortal1, typename IntPortal2, typename IntPortal3, typename StatePortalType>
  VISKORES_EXEC void operator()(const viskores::Id& idx,
                                const Point3D& seedInd,
                                const FieldExecType& field,
                                IntPortal1& stateCounts,
                                IntPortal2& endRegions,
                                IntPortal3& statusCodes,
                                StatePortalType& statesOut) const
  {
    const viskores::Id stateBase = idx * static_cast<viskores::Id>(maxStatesPerSeed_);

    stateCounts.Set(idx, 0);
    endRegions.Set(idx, 0);
    statusCodes.Set(idx, static_cast<viskores::Id>(TraceStatus::Ok));

    if (maxStatesPerSeed_ <= 0)
    {
      statusCodes.Set(idx, static_cast<viskores::Id>(TraceStatus::InvalidConfiguration));
      endRegions.Set(idx, 98);
      return;
    }

    const FloatType xindSeed = seedInd.x;
    if (xindSeed < 1.0 || xindSeed > static_cast<FloatType>(field.nx))
    {
      statusCodes.Set(idx, static_cast<viskores::Id>(TraceStatus::InvalidSeed));
      endRegions.Set(idx, 99);
      return;
    }

    FloatType xind = xindSeed;
    XZPoint current;
    current.x = field.interp1(field.xiarray, field.xarray, xind);

    int yStart = static_cast<int>(std::llround(seedInd.y));
    if (yStart < 1)
    {
      yStart = 1;
    }
    else if (yStart > field.ny)
    {
      yStart = field.ny;
    }

    const FloatType zind0 = seedInd.z;
    current.z = field.interp1(field.ziarray, field.zarray, zind0);

    int region = 1;
    if (xind < static_cast<FloatType>(field.ixsep) + 0.5)
    {
      region = 0;
      if (yStart < field.nypf1 + 1 || yStart > field.nypf2)
      {
        region = 2;
      }
    }

    viskores::Id localStateCount = 0;

    TrajectoryState initial;
    initial.turn = 1;
    initial.ind.x = xind;
    initial.ind.y = static_cast<FloatType>(yStart);
    initial.ind.z = zind0;
    initial.region = region;
    initial.segmentLength = 0.0;
    initial.rawZ = current.z;

    statesOut.Set(stateBase + static_cast<viskores::Id>(localStateCount), initial);
    ++localStateCount;

    int iturn = 1;
    while (region < 10 && localStateCount < maxStatesPerSeed_)
    {
      for (int iy = 0; iy < field.ny - 1; ++iy)
      {
        if (region >= 10 || localStateCount >= maxStatesPerSeed_)
        {
          break;
        }

        XZPoint next = current;
        int yEnd = yStart;

        if (region == 0 && yStart > field.nypf1 && yStart < field.nypf2 + 1)
        {
          rk4Step(field, current, yStart, region, direction_, next);
          const FloatType rawZEnd = next.z;

          yEnd = (direction_ == 1) ? (yStart + 1) : (yStart - 1);

          if (next.x > field.xMax)
          {
            region = 12;
          }
          else if (next.x < field.xMin)
          {
            region = 11;
          }
          else
          {
            xind = field.interp1(field.xarray, field.xiarray, next.x);
            if (xind > static_cast<FloatType>(field.ixsep) + 0.5)
            {
              region = 1;
            }
          }

          if (direction_ == 1 && yStart == field.nypf2 && region == 0)
          {
            const FloatType shift = field.interp1(field.xiarray, field.shiftAngle, xind);
            next.z += shift;
            yEnd = field.nypf1 + 1;
          }
          if (direction_ == -1 && yStart == (field.nypf1 + 1) && region == 0)
          {
            const FloatType shift = field.interp1(field.xiarray, field.shiftAngle, xind);
            next.z -= shift;
            yEnd = field.nypf2;
          }

          next.z = field.wrapZ(next.z);
          const FloatType zind = field.interp1(field.zarray, field.ziarray, next.z);

          TrajectoryState step;
          step.turn = iturn;
          step.ind.x = xind;
          step.ind.y = static_cast<FloatType>(yEnd);
          step.ind.z = zind;
          step.region = region;
          step.segmentLength = 0.0;
          step.rawZ = rawZEnd;

          statesOut.Set(stateBase + static_cast<viskores::Id>(localStateCount), step);
          ++localStateCount;

          current = next;
          yStart = yEnd;
        }
        else if (region == 1 || region == 2)
        {
          rk4Step(field, current, yStart, region, direction_, next);
          const FloatType rawZEnd = next.z;

          yEnd = (direction_ == 1) ? (yStart + 1) : (yStart - 1);

          if (direction_ == 1 && yStart == field.nypf1 && region == 2)
          {
            yEnd = field.nypf2 + 1;
          }
          else if (direction_ == -1 && yStart == field.nypf2 + 1 && region == 2)
          {
            yEnd = field.nypf1;
          }

          if (next.x > field.xMax)
          {
            region = 12;
          }
          else if (next.x < field.xMin)
          {
            region = 11;
          }
          else
          {
            xind = field.interp1(field.xarray, field.xiarray, next.x);
            if (xind < static_cast<FloatType>(field.ixsep) + 0.5 && yEnd > field.nypf1 && yEnd < field.nypf2 + 1)
            {
              region = 0;
            }
            else if (xind < static_cast<FloatType>(field.ixsep) + 0.5 && (yEnd > field.nypf2 - 1 || yEnd < field.nypf1))
            {
              region = 2;
            }
          }

          if (direction_ == 1 && yEnd == field.ny)
          {
            region = 14;
          }
          else if (direction_ == -1 && yEnd == 1)
          {
            region = 13;
          }

          next.z = field.wrapZ(next.z);
          const FloatType zind = field.interp1(field.zarray, field.ziarray, next.z);

          TrajectoryState step;
          step.turn = iturn;
          step.ind.x = xind;
          step.ind.y = static_cast<FloatType>(yEnd);
          step.ind.z = zind;
          step.region = region;
          step.segmentLength = 0.0;
          step.rawZ = rawZEnd;

          statesOut.Set(stateBase + static_cast<viskores::Id>(localStateCount), step);
          ++localStateCount;

          current = next;
          yStart = yEnd;
        }
        else
        {
          break;
        }
      }

      ++iturn;
    }

    stateCounts.Set(idx, localStateCount);
    endRegions.Set(idx, region);
    if (region < 10)
    {
      statusCodes.Set(idx, static_cast<viskores::Id>(TraceStatus::MaxStepLimitReached));
    }
  }

private:
  template <typename FieldExecType>
  VISKORES_EXEC void rk4Step(const FieldExecType& field, const XZPoint& start, int yStart, int region, int direction, XZPoint& end) const
  {
    constexpr FloatType h = 1.0;
    const FloatType hh = 0.5 * h;
    const FloatType h6 = h / 6.0;

    XZDeriv k1;
    XZDeriv k2;
    XZDeriv k3;
    XZDeriv k4;

    field.evaluateStage(start, yStart, region, direction, 0, k1);
    const XZPoint p1{ start.x + direction * hh * k1.dxdy, start.z + direction * hh * k1.dzdy };

    field.evaluateStage(p1, yStart, region, direction, 1, k2);
    const XZPoint p2{ start.x + direction * hh * k2.dxdy, start.z + direction * hh * k2.dzdy };

    field.evaluateStage(p2, yStart, region, direction, 1, k3);
    const XZPoint p3{ start.x + direction * k3.dxdy, start.z + direction * k3.dzdy };

    field.evaluateStage(p3, yStart, region, direction, 2, k4);

    end.x = start.x + direction * h6 * (k1.dxdy + 2.0 * k2.dxdy + 2.0 * k3.dxdy + k4.dxdy);
    end.z = start.z + direction * h6 * (k1.dzdy + 2.0 * k2.dzdy + 2.0 * k3.dzdy + k4.dzdy);
  }

  viskores::Id maxStatesPerSeed_ = 1;
  int direction_ = 1;
};

#endif

#endif
