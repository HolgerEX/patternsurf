#pragma once

#include "pattern/core/result.hpp"
#include "pattern/partition/partition.hpp"

#include <cstdint>
#include <string_view>
#include <vector>

namespace pattern {

struct Polyline {
    std::vector<ParameterPoint> points;
};

struct CurveSet {
    std::vector<Polyline> curves;
};

struct VectorFieldConfig {
    double seed_spacing{0.1};
    double step_size{0.02};
    double maximum_length{2.0};
    double base_angle{};
    double warp_amplitude{};
    double warp_frequency{1.0};
};

struct GeneratorDescriptor {
    std::string_view type_id;
    std::uint32_t schema_version;
    std::string_view display_name;
    bool deterministic;
};

class VectorFieldCurveGenerator {
public:
    [[nodiscard]] static Result<VectorFieldCurveGenerator> create(
        VectorFieldConfig config);

    [[nodiscard]] static constexpr GeneratorDescriptor descriptor() noexcept
    {
        return {
            "pattern.vector_field.curves",
            1,
            "Vector-field curves",
            true
        };
    }

    [[nodiscard]] Result<CurveSet> generate(
        const Partition& partition,
        RegionId region) const;

private:
    explicit VectorFieldCurveGenerator(VectorFieldConfig config);

    VectorFieldConfig config_;
};

}  // namespace pattern
