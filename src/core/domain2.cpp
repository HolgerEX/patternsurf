#include "pattern/surface/domain2.hpp"

#include <cmath>
#include <string>
#include <utility>

namespace pattern {
namespace {

Diagnostic invalid_domain(std::string message)
{
    return {
        ErrorCode::invalid_argument,
        Severity::error,
        std::move(message)
    };
}

}  // namespace

Domain2::Domain2(
    double min_u,
    double min_v,
    double max_u,
    double max_v)
    : min_u_(min_u),
      min_v_(min_v),
      max_u_(max_u),
      max_v_(max_v)
{
}

Result<Domain2> Domain2::rectangle(
    double min_u,
    double min_v,
    double max_u,
    double max_v)
{
    if (!std::isfinite(min_u) || !std::isfinite(min_v)
        || !std::isfinite(max_u) || !std::isfinite(max_v)) {
        return Result<Domain2>::failure(
            invalid_domain("Rectangle bounds must be finite."));
    }

    if (!(min_u < max_u) || !(min_v < max_v)) {
        return Result<Domain2>::failure(
            invalid_domain("Rectangle bounds must have positive extent."));
    }

    const double width = max_u - min_u;
    const double height = max_v - min_v;
    const double area = width * height;
    if (!std::isfinite(width) || !std::isfinite(height)
        || !std::isfinite(area) || area <= 0.0) {
        return Result<Domain2>::failure(
            invalid_domain("Rectangle extents and area must be finite and positive."));
    }

    return Result<Domain2>::success(
        Domain2(min_u, min_v, max_u, max_v));
}

bool Domain2::contains(ParameterPoint point) const noexcept
{
    return std::isfinite(point.u)
        && std::isfinite(point.v)
        && point.u >= min_u_
        && point.u <= max_u_
        && point.v >= min_v_
        && point.v <= max_v_;
}

double Domain2::min_u() const noexcept
{
    return min_u_;
}

double Domain2::min_v() const noexcept
{
    return min_v_;
}

double Domain2::max_u() const noexcept
{
    return max_u_;
}

double Domain2::max_v() const noexcept
{
    return max_v_;
}

double Domain2::width() const noexcept
{
    return max_u_ - min_u_;
}

double Domain2::height() const noexcept
{
    return max_v_ - min_v_;
}

}  // namespace pattern
