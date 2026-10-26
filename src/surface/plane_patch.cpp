#include "pattern/surface/plane_patch.hpp"

#include <utility>

namespace pattern {

PlanePatch::PlanePatch(Domain2 domain)
    : domain_(std::move(domain))
{
}

const Domain2& PlanePatch::parameter_domain() const noexcept
{
    return domain_;
}

Result<SurfaceSample> PlanePatch::evaluate(ParameterPoint point) const
{
    if (!domain_.contains(point)) {
        return Result<SurfaceSample>::failure({
            ErrorCode::outside_domain,
            Severity::error,
            "Point is outside the plane patch parameter domain."
        });
    }

    return Result<SurfaceSample>::success({
        {point.u, point.v, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}
    });
}

}  // namespace pattern
