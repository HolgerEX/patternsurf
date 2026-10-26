#pragma once

#include "pattern/surface/surface_patch.hpp"

namespace pattern {

class PlanePatch final : public SurfacePatch {
public:
    explicit PlanePatch(Domain2 domain);

    [[nodiscard]] const Domain2& parameter_domain() const noexcept override;
    [[nodiscard]] Result<SurfaceSample> evaluate(
        ParameterPoint point) const override;

private:
    Domain2 domain_;
};

}  // namespace pattern
