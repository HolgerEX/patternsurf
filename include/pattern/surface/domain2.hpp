#pragma once

#include "pattern/core/result.hpp"

namespace pattern {

struct ParameterPoint {
    double u;
    double v;
};

class Domain2 {
public:
    [[nodiscard]] static Result<Domain2> rectangle(
        double min_u,
        double min_v,
        double max_u,
        double max_v);

    [[nodiscard]] bool contains(ParameterPoint point) const noexcept;

    [[nodiscard]] double min_u() const noexcept;
    [[nodiscard]] double min_v() const noexcept;
    [[nodiscard]] double max_u() const noexcept;
    [[nodiscard]] double max_v() const noexcept;
    [[nodiscard]] double width() const noexcept;
    [[nodiscard]] double height() const noexcept;

private:
    Domain2(double min_u, double min_v, double max_u, double max_v);

    double min_u_;
    double min_v_;
    double max_u_;
    double max_v_;
};

}  // namespace pattern
