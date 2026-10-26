#pragma once

#include "pattern/core/id.hpp"
#include "pattern/core/result.hpp"
#include "pattern/surface/domain2.hpp"

#include <map>
#include <optional>
#include <vector>

namespace pattern {

struct BoundaryUse {
    BoundaryId boundary;
    bool forward;
};

struct Boundary {
    BoundaryId id;
    VertexId from;
    VertexId to;
    std::optional<RegionId> left_region;
    std::optional<RegionId> right_region;
};

struct Region {
    RegionId id;
    std::vector<BoundaryUse> boundary;
    bool background;
};

class Partition {
public:
    [[nodiscard]] static Partition create(Domain2 domain);

    [[nodiscard]] const Domain2& domain() const noexcept;
    [[nodiscard]] const std::map<VertexId, ParameterPoint>& vertices() const noexcept;
    [[nodiscard]] const std::map<BoundaryId, Boundary>& boundaries() const noexcept;
    [[nodiscard]] const std::map<RegionId, Region>& regions() const noexcept;
    [[nodiscard]] const Region* find_region(RegionId id) const noexcept;
    [[nodiscard]] Result<std::vector<RegionId>> adjacent_regions(RegionId id) const;

    [[nodiscard]] std::vector<Diagnostic> validate() const;

    [[nodiscard]] Result<VertexId> move_vertex(
        VertexId vertex,
        ParameterPoint position);
    [[nodiscard]] Result<VertexId> insert_boundary_vertex(
        BoundaryId boundary,
        ParameterPoint position);
    [[nodiscard]] Result<RegionId> split_region(
        RegionId region,
        VertexId first,
        VertexId second);
    [[nodiscard]] Result<RegionId> merge_regions(
        RegionId retained,
        RegionId removed);
    [[nodiscard]] Result<RegionId> delete_region(RegionId region);

private:
    explicit Partition(Domain2 domain);

    [[nodiscard]] Result<VertexId> move_vertex_in_place(
        VertexId vertex,
        ParameterPoint position);
    [[nodiscard]] Result<VertexId> insert_boundary_vertex_in_place(
        BoundaryId boundary,
        ParameterPoint position);
    [[nodiscard]] Result<RegionId> split_region_in_place(
        RegionId region,
        VertexId first,
        VertexId second);
    [[nodiscard]] Result<RegionId> merge_regions_in_place(
        RegionId retained,
        RegionId removed);
    [[nodiscard]] Result<RegionId> delete_region_in_place(RegionId region);

    Domain2 domain_;
    std::map<VertexId, ParameterPoint> vertices_;
    std::map<BoundaryId, Boundary> boundaries_;
    std::map<RegionId, Region> regions_;
    std::uint64_t next_region_id_{1};
    std::uint64_t next_vertex_id_{1};
    std::uint64_t next_boundary_id_{1};
};

}  // namespace pattern
