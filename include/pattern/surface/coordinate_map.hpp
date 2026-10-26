#pragma once

#include "pattern/core/result.hpp"
#include "pattern/surface/domain2.hpp"

namespace pattern {

struct CoordinatePoint2 {
    double h;
    double k;
};

struct Vector2 {
    double x;
    double y;
};

struct Jacobian2 {
    double dh_du;
    double dh_dv;
    double dk_du;
    double dk_dv;
};

struct CoordinateMapCapabilities {
    bool time_dependent;
    bool has_jacobian;
    bool has_inverse;
    bool globally_invertible;
};

class CoordinateMap2D {
public:
    virtual ~CoordinateMap2D() = default;

    [[nodiscard]] virtual CoordinateMapCapabilities capabilities() const noexcept = 0;
    [[nodiscard]] virtual Result<CoordinatePoint2> evaluate(
        ParameterPoint point) const = 0;
    [[nodiscard]] virtual Result<Jacobian2> jacobian(
        ParameterPoint point) const = 0;
    [[nodiscard]] virtual Result<ParameterPoint> inverse(
        CoordinatePoint2 point) const = 0;
};

class IdentityCoordinateMap final : public CoordinateMap2D {
public:
    [[nodiscard]] CoordinateMapCapabilities capabilities() const noexcept override;
    [[nodiscard]] Result<CoordinatePoint2> evaluate(
        ParameterPoint point) const override;
    [[nodiscard]] Result<Jacobian2> jacobian(
        ParameterPoint point) const override;
    [[nodiscard]] Result<ParameterPoint> inverse(
        CoordinatePoint2 point) const override;
};

class AffineCoordinateMap final : public CoordinateMap2D {
public:
    [[nodiscard]] static Result<AffineCoordinateMap> create(
        double dh_du,
        double dh_dv,
        double dk_du,
        double dk_dv,
        double offset_h,
        double offset_k);

    [[nodiscard]] CoordinateMapCapabilities capabilities() const noexcept override;
    [[nodiscard]] Result<CoordinatePoint2> evaluate(
        ParameterPoint point) const override;
    [[nodiscard]] Result<Jacobian2> jacobian(
        ParameterPoint point) const override;
    [[nodiscard]] Result<ParameterPoint> inverse(
        CoordinatePoint2 point) const override;

private:
    AffineCoordinateMap(
        Jacobian2 jacobian,
        double offset_h,
        double offset_k,
        double matrix_scale,
        double normalized_determinant,
        bool invertible);

    Jacobian2 jacobian_;
    double offset_h_;
    double offset_k_;
    double matrix_scale_;
    double normalized_determinant_;
    bool invertible_;
};

}  // namespace pattern
