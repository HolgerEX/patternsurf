#pragma once

#include "pattern/core/result.hpp"
#include "pattern/core/vector3.hpp"
#include "pattern/surface/domain2.hpp"

namespace pattern {

struct SurfaceSample {
    Vector3 position;
    Vector3 du;
    Vector3 dv;
};

class SurfacePatch {
public:
    virtual ~SurfacePatch() = default;

    [[nodiscard]] virtual const Domain2& parameter_domain() const noexcept = 0;
    [[nodiscard]] virtual Result<SurfaceSample> evaluate(
        ParameterPoint point) const = 0;
};

}  // namespace pattern
