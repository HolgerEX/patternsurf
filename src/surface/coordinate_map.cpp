#include "pattern/surface/coordinate_map.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

namespace pattern {
namespace {

Diagnostic diagnostic(ErrorCode code, std::string message)
{
    return {code, Severity::error, std::move(message)};
}

bool finite(ParameterPoint point)
{
    return std::isfinite(point.u) && std::isfinite(point.v);
}

bool finite(CoordinatePoint2 point)
{
    return std::isfinite(point.h) && std::isfinite(point.k);
}

}  // namespace

CoordinateMapCapabilities IdentityCoordinateMap::capabilities() const noexcept
{
    return {false, true, true, true};
}

Result<CoordinatePoint2> IdentityCoordinateMap::evaluate(
    ParameterPoint point) const
{
    if (!finite(point)) {
        return Result<CoordinatePoint2>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Identity map input must be finite."));
    }
    return Result<CoordinatePoint2>::success({point.u, point.v});
}

Result<Jacobian2> IdentityCoordinateMap::jacobian(ParameterPoint point) const
{
    if (!finite(point)) {
        return Result<Jacobian2>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Identity map input must be finite."));
    }
    return Result<Jacobian2>::success({1.0, 0.0, 0.0, 1.0});
}

Result<ParameterPoint> IdentityCoordinateMap::inverse(
    CoordinatePoint2 point) const
{
    if (!finite(point)) {
        return Result<ParameterPoint>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Identity map inverse input must be finite."));
    }
    return Result<ParameterPoint>::success({point.h, point.k});
}

AffineCoordinateMap::AffineCoordinateMap(
    Jacobian2 jacobian,
    double offset_h,
    double offset_k,
    double matrix_scale,
    double normalized_determinant,
    bool invertible)
    : jacobian_(jacobian),
      offset_h_(offset_h),
      offset_k_(offset_k),
      matrix_scale_(matrix_scale),
      normalized_determinant_(normalized_determinant),
      invertible_(invertible)
{
}

Result<AffineCoordinateMap> AffineCoordinateMap::create(
    double dh_du,
    double dh_dv,
    double dk_du,
    double dk_dv,
    double offset_h,
    double offset_k)
{
    if (!std::isfinite(dh_du) || !std::isfinite(dh_dv)
        || !std::isfinite(dk_du) || !std::isfinite(dk_dv)
        || !std::isfinite(offset_h) || !std::isfinite(offset_k)) {
        return Result<AffineCoordinateMap>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Affine coordinate-map coefficients and offsets must be finite."));
    }

    const double matrix_scale = std::max({
        std::abs(dh_du),
        std::abs(dh_dv),
        std::abs(dk_du),
        std::abs(dk_dv)
    });
    double normalized_determinant = 0.0;
    if (matrix_scale > 0.0) {
        normalized_determinant =
            (dh_du / matrix_scale) * (dk_dv / matrix_scale)
            - (dh_dv / matrix_scale) * (dk_du / matrix_scale);
    }
    const bool invertible =
        std::abs(normalized_determinant)
        > std::numeric_limits<double>::epsilon() * 32.0;

    return Result<AffineCoordinateMap>::success(AffineCoordinateMap(
        {dh_du, dh_dv, dk_du, dk_dv},
        offset_h,
        offset_k,
        matrix_scale,
        normalized_determinant,
        invertible));
}

CoordinateMapCapabilities AffineCoordinateMap::capabilities() const noexcept
{
    return {false, true, invertible_, invertible_};
}

Result<CoordinatePoint2> AffineCoordinateMap::evaluate(
    ParameterPoint point) const
{
    if (!finite(point)) {
        return Result<CoordinatePoint2>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Affine map input must be finite."));
    }

    const CoordinatePoint2 mapped{
        jacobian_.dh_du * point.u + jacobian_.dh_dv * point.v + offset_h_,
        jacobian_.dk_du * point.u + jacobian_.dk_dv * point.v + offset_k_
    };
    if (!finite(mapped)) {
        return Result<CoordinatePoint2>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Affine map produced a non-finite coordinate."));
    }
    return Result<CoordinatePoint2>::success(mapped);
}

Result<Jacobian2> AffineCoordinateMap::jacobian(ParameterPoint point) const
{
    if (!finite(point)) {
        return Result<Jacobian2>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Affine map input must be finite."));
    }
    return Result<Jacobian2>::success(jacobian_);
}

Result<ParameterPoint> AffineCoordinateMap::inverse(
    CoordinatePoint2 point) const
{
    if (!finite(point)) {
        return Result<ParameterPoint>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Affine map inverse input must be finite."));
    }
    if (!invertible_) {
        return Result<ParameterPoint>::failure(diagnostic(
            ErrorCode::unsupported_operation,
            "The affine map is singular or numerically non-invertible."));
    }

    const double h = (point.h - offset_h_) / matrix_scale_;
    const double k = (point.k - offset_k_) / matrix_scale_;
    const double a = jacobian_.dh_du / matrix_scale_;
    const double b = jacobian_.dh_dv / matrix_scale_;
    const double c = jacobian_.dk_du / matrix_scale_;
    const double d = jacobian_.dk_dv / matrix_scale_;
    const double u = (d * h - b * k) / normalized_determinant_;
    const double v = (a * k - c * h) / normalized_determinant_;
    if (!std::isfinite(u) || !std::isfinite(v)) {
        return Result<ParameterPoint>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Affine inverse produced a non-finite parameter point."));
    }

    return Result<ParameterPoint>::success({u, v});
}

}  // namespace pattern
