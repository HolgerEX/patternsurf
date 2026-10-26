#include "pattern/surface/cylinder_patch.hpp"

#include <cmath>
#include <utility>

namespace pattern {

CylinderPatch::CylinderPatch(double radius, Domain2 domain)
    : radius_(radius),
      domain_(std::move(domain))
{
}

Result<CylinderPatch> CylinderPatch::create(double radius, Domain2 domain)
{
    if (!std::isfinite(radius) || radius <= 0.0) {
        return Result<CylinderPatch>::failure({
            ErrorCode::invalid_argument,
            Severity::error,
            "Cylinder radius must be finite and positive."
        });
    }

    return Result<CylinderPatch>::success(
        CylinderPatch(radius, std::move(domain)));
}

double CylinderPatch::radius() const noexcept
{
    return radius_;
}

const Domain2& CylinderPatch::parameter_domain() const noexcept
{
    return domain_;
}

Result<SurfaceSample> CylinderPatch::evaluate(ParameterPoint point) const
{
    if (!domain_.contains(point)) {
        return Result<SurfaceSample>::failure({
            ErrorCode::outside_domain,
            Severity::error,
            "Point is outside the cylinder patch parameter domain."
        });
    }

    const double cosine = std::cos(point.u);
    const double sine = std::sin(point.u);

    return Result<SurfaceSample>::success({
        {radius_ * cosine, radius_ * sine, point.v},
        {-radius_ * sine, radius_ * cosine, 0.0},
        {0.0, 0.0, 1.0}
    });
}

}  // namespace pattern
