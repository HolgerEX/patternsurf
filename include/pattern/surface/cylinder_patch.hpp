#pragma once

#include "pattern/surface/surface_patch.hpp"

namespace pattern {

class CylinderPatch final : public SurfacePatch {
public:
    [[nodiscard]] static Result<CylinderPatch> create(
        double radius,
        Domain2 domain);

    [[nodiscard]] double radius() const noexcept;
    [[nodiscard]] const Domain2& parameter_domain() const noexcept override;
    [[nodiscard]] Result<SurfaceSample> evaluate(
        ParameterPoint point) const override;

private:
    CylinderPatch(double radius, Domain2 domain);

    double radius_;
    Domain2 domain_;
};

}  // namespace pattern
