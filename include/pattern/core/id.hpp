#pragma once

#include <compare>
#include <cstdint>

namespace pattern {

template <typename Tag>
struct Id {
    std::uint64_t value{};

    explicit constexpr Id(std::uint64_t id_value)
        : value(id_value)
    {
    }

    constexpr auto operator<=>(const Id&) const = default;
};

struct RegionIdTag;
struct VertexIdTag;
struct BoundaryIdTag;

using RegionId = Id<RegionIdTag>;
using VertexId = Id<VertexIdTag>;
using BoundaryId = Id<BoundaryIdTag>;

}  // namespace pattern
