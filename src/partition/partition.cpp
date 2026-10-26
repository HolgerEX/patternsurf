#include "pattern/partition/partition.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <numeric>
#include <set>
#include <string>
#include <utility>

namespace pattern {
namespace {

Diagnostic diagnostic(ErrorCode code, std::string message)
{
    return {code, Severity::error, std::move(message)};
}

struct GeometryTolerance {
    double linear;
    double cross;
    double area;
};

GeometryTolerance geometry_tolerance(const Domain2& domain)
{
    const double extent = std::max(domain.width(), domain.height());
    const double coordinate_scale = std::max({
        std::abs(domain.min_u()),
        std::abs(domain.min_v()),
        std::abs(domain.max_u()),
        std::abs(domain.max_v())
    });
    const double linear = std::max(
        extent * 1e-12,
        coordinate_scale * std::numeric_limits<double>::epsilon() * 8.0);
    const double domain_area = domain.width() * domain.height();
    return {
        linear,
        linear * extent,
        std::max(domain_area * 1e-12, linear * extent * 2.0)
    };
}

double cross(ParameterPoint a, ParameterPoint b, ParameterPoint c)
{
    return (b.u - a.u) * (c.v - a.v) - (b.v - a.v) * (c.u - a.u);
}

double distance(ParameterPoint a, ParameterPoint b)
{
    return std::hypot(a.u - b.u, a.v - b.v);
}

VertexId use_start_vertex(const Boundary& boundary, const BoundaryUse& use)
{
    return use.forward ? boundary.from : boundary.to;
}

VertexId use_end_vertex(const Boundary& boundary, const BoundaryUse& use)
{
    return use.forward ? boundary.to : boundary.from;
}

std::vector<ParameterPoint> region_points(
    const Region& region,
    const std::map<BoundaryId, Boundary>& boundaries,
    const std::map<VertexId, ParameterPoint>& vertices)
{
    std::vector<ParameterPoint> points;
    points.reserve(region.boundary.size());
    for (const BoundaryUse& use : region.boundary) {
        const Boundary& boundary = boundaries.at(use.boundary);
        points.push_back(vertices.at(use_start_vertex(boundary, use)));
    }
    return points;
}

double signed_area(const std::vector<ParameterPoint>& points)
{
    if (points.size() < 3) {
        return 0.0;
    }

    const ParameterPoint origin = points.front();
    double twice_area = 0.0;
    for (std::size_t index = 1; index + 1 < points.size(); ++index) {
        twice_area += cross(origin, points[index], points[index + 1]);
    }
    return 0.5 * twice_area;
}

bool point_on_segment(
    ParameterPoint point,
    ParameterPoint first,
    ParameterPoint second,
    double linear_tolerance,
    double cross_tolerance)
{
    if (std::abs(cross(first, second, point)) > cross_tolerance) {
        return false;
    }

    return point.u >= std::min(first.u, second.u) - linear_tolerance
        && point.u <= std::max(first.u, second.u) + linear_tolerance
        && point.v >= std::min(first.v, second.v) - linear_tolerance
        && point.v <= std::max(first.v, second.v) + linear_tolerance;
}

int orientation(
    ParameterPoint first,
    ParameterPoint second,
    ParameterPoint third,
    double cross_tolerance)
{
    const double value = cross(first, second, third);
    if (value > cross_tolerance) {
        return 1;
    }
    if (value < -cross_tolerance) {
        return -1;
    }
    return 0;
}

bool segments_intersect(
    ParameterPoint a,
    ParameterPoint b,
    ParameterPoint c,
    ParameterPoint d,
    double linear_tolerance,
    double cross_tolerance)
{
    const int abc = orientation(a, b, c, cross_tolerance);
    const int abd = orientation(a, b, d, cross_tolerance);
    const int cda = orientation(c, d, a, cross_tolerance);
    const int cdb = orientation(c, d, b, cross_tolerance);

    if (abc * abd < 0 && cda * cdb < 0) {
        return true;
    }

    return (abc == 0 && point_on_segment(c, a, b, linear_tolerance, cross_tolerance))
        || (abd == 0 && point_on_segment(d, a, b, linear_tolerance, cross_tolerance))
        || (cda == 0 && point_on_segment(a, c, d, linear_tolerance, cross_tolerance))
        || (cdb == 0 && point_on_segment(b, c, d, linear_tolerance, cross_tolerance));
}

bool point_in_polygon(
    ParameterPoint point,
    const std::vector<ParameterPoint>& polygon,
    double linear_tolerance,
    double cross_tolerance)
{
    bool inside = false;
    for (std::size_t index = 0; index < polygon.size(); ++index) {
        const ParameterPoint first = polygon[index];
        const ParameterPoint second = polygon[(index + 1) % polygon.size()];
        if (point_on_segment(point, first, second, linear_tolerance, cross_tolerance)) {
            return false;
        }

        const bool straddles = (first.v > point.v) != (second.v > point.v);
        if (straddles) {
            const double crossing_u = first.u
                + (point.v - first.v) * (second.u - first.u) / (second.v - first.v);
            if (point.u < crossing_u) {
                inside = !inside;
            }
        }
    }
    return inside;
}

bool lies_on_domain_boundary(
    ParameterPoint first,
    ParameterPoint second,
    const Domain2& domain,
    double tolerance)
{
    return (std::abs(first.u - domain.min_u()) <= tolerance
            && std::abs(second.u - domain.min_u()) <= tolerance)
        || (std::abs(first.u - domain.max_u()) <= tolerance
            && std::abs(second.u - domain.max_u()) <= tolerance)
        || (std::abs(first.v - domain.min_v()) <= tolerance
            && std::abs(second.v - domain.min_v()) <= tolerance)
        || (std::abs(first.v - domain.max_v()) <= tolerance
            && std::abs(second.v - domain.max_v()) <= tolerance);
}

}  // namespace

Partition::Partition(Domain2 domain)
    : domain_(std::move(domain))
{
}

Partition Partition::create(Domain2 domain)
{
    Partition partition(std::move(domain));
    const RegionId region_id{partition.next_region_id_++};

    const VertexId bottom_left{partition.next_vertex_id_++};
    const VertexId bottom_right{partition.next_vertex_id_++};
    const VertexId top_right{partition.next_vertex_id_++};
    const VertexId top_left{partition.next_vertex_id_++};
    partition.vertices_.emplace(
        bottom_left,
        ParameterPoint{partition.domain_.min_u(), partition.domain_.min_v()});
    partition.vertices_.emplace(
        bottom_right,
        ParameterPoint{partition.domain_.max_u(), partition.domain_.min_v()});
    partition.vertices_.emplace(
        top_right,
        ParameterPoint{partition.domain_.max_u(), partition.domain_.max_v()});
    partition.vertices_.emplace(
        top_left,
        ParameterPoint{partition.domain_.min_u(), partition.domain_.max_v()});

    const std::vector<VertexId> cycle{
        bottom_left, bottom_right, top_right, top_left
    };
    std::vector<BoundaryUse> uses;
    uses.reserve(cycle.size());
    for (std::size_t index = 0; index < cycle.size(); ++index) {
        const BoundaryId boundary_id{partition.next_boundary_id_++};
        partition.boundaries_.emplace(
            boundary_id,
            Boundary{
                boundary_id,
                cycle[index],
                cycle[(index + 1) % cycle.size()],
                region_id,
                std::nullopt
            });
        uses.push_back({boundary_id, true});
    }

    partition.regions_.emplace(region_id, Region{region_id, std::move(uses), false});
    return partition;
}

const Domain2& Partition::domain() const noexcept
{
    return domain_;
}

const std::map<VertexId, ParameterPoint>& Partition::vertices() const noexcept
{
    return vertices_;
}

const std::map<BoundaryId, Boundary>& Partition::boundaries() const noexcept
{
    return boundaries_;
}

const std::map<RegionId, Region>& Partition::regions() const noexcept
{
    return regions_;
}

const Region* Partition::find_region(RegionId id) const noexcept
{
    const auto found = regions_.find(id);
    return found == regions_.end() ? nullptr : &found->second;
}

Result<std::vector<RegionId>> Partition::adjacent_regions(RegionId id) const
{
    if (!regions_.contains(id)) {
        return Result<std::vector<RegionId>>::failure(diagnostic(
            ErrorCode::not_found,
            "Cannot query adjacency for a missing region."));
    }

    std::set<RegionId> adjacent;
    for (const auto& [boundary_id, boundary] : boundaries_) {
        (void)boundary_id;
        if (boundary.left_region == id && boundary.right_region) {
            adjacent.insert(*boundary.right_region);
        } else if (boundary.right_region == id && boundary.left_region) {
            adjacent.insert(*boundary.left_region);
        }
    }

    return Result<std::vector<RegionId>>::success(
        std::vector<RegionId>(adjacent.begin(), adjacent.end()));
}

std::vector<Diagnostic> Partition::validate() const
{
    std::vector<Diagnostic> issues;
    if (regions_.empty()) {
        issues.push_back(diagnostic(
            ErrorCode::invalid_topology,
            "A full-coverage partition must contain at least one region."));
        return issues;
    }

    const GeometryTolerance tolerance = geometry_tolerance(domain_);

    for (const auto& [id, point] : vertices_) {
        if (id.value == 0 || !std::isfinite(point.u) || !std::isfinite(point.v)
            || !domain_.contains(point)) {
            issues.push_back(diagnostic(
                ErrorCode::invalid_topology,
                "A vertex has an invalid ID or lies outside the parameter domain."));
            return issues;
        }
    }

    for (auto first = vertices_.begin(); first != vertices_.end(); ++first) {
        for (auto second = std::next(first); second != vertices_.end(); ++second) {
            if (distance(first->second, second->second) <= tolerance.linear) {
                issues.push_back(diagnostic(
                    ErrorCode::invalid_topology,
                    "Distinct vertices occupy the same parameter-space position."));
                return issues;
            }
        }
    }

    std::map<BoundaryId, std::size_t> forward_uses;
    std::map<BoundaryId, std::size_t> reverse_uses;
    for (const auto& [region_id, region] : regions_) {
        if (region_id.value == 0 || region.id != region_id
            || region.boundary.size() < 3) {
            issues.push_back(diagnostic(
                ErrorCode::invalid_topology,
                "A region has an invalid ID or fewer than three boundary edges."));
            return issues;
        }

        std::set<BoundaryId> unique_boundaries;
        for (std::size_t index = 0; index < region.boundary.size(); ++index) {
            const BoundaryUse& use = region.boundary[index];
            const auto boundary_it = boundaries_.find(use.boundary);
            if (boundary_it == boundaries_.end()
                || !unique_boundaries.insert(use.boundary).second) {
                issues.push_back(diagnostic(
                    ErrorCode::invalid_topology,
                    "A region references a missing or repeated boundary."));
                return issues;
            }

            const Boundary& boundary = boundary_it->second;
            const VertexId end = use_end_vertex(boundary, use);
            const BoundaryUse& next_use =
                region.boundary[(index + 1) % region.boundary.size()];
            const auto next_boundary = boundaries_.find(next_use.boundary);
            if (next_boundary == boundaries_.end()
                || end != use_start_vertex(next_boundary->second, next_use)) {
                issues.push_back(diagnostic(
                    ErrorCode::invalid_topology,
                    "A region boundary cycle is not connected."));
                return issues;
            }

            const std::optional<RegionId>& expected_side =
                use.forward ? boundary.left_region : boundary.right_region;
            if (!expected_side || *expected_side != region_id) {
                issues.push_back(diagnostic(
                    ErrorCode::invalid_topology,
                    "A boundary use disagrees with its incident region."));
                return issues;
            }
            if (use.forward) {
                ++forward_uses[use.boundary];
            } else {
                ++reverse_uses[use.boundary];
            }
        }
    }

    for (const auto& [boundary_id, boundary] : boundaries_) {
        if (boundary_id.value == 0 || boundary.id != boundary_id
            || boundary.from == boundary.to
            || !vertices_.contains(boundary.from)
            || !vertices_.contains(boundary.to)
            || (!boundary.left_region && !boundary.right_region)
            || (boundary.left_region && boundary.right_region
                && boundary.left_region == boundary.right_region)) {
            issues.push_back(diagnostic(
                ErrorCode::invalid_topology,
                "A boundary has invalid endpoints or incident regions."));
            return issues;
        }

        const std::size_t expected_forward = boundary.left_region ? 1U : 0U;
        const std::size_t expected_reverse = boundary.right_region ? 1U : 0U;
        if (forward_uses[boundary_id] != expected_forward
            || reverse_uses[boundary_id] != expected_reverse
            || (boundary.left_region && !regions_.contains(*boundary.left_region))
            || (boundary.right_region && !regions_.contains(*boundary.right_region))) {
            issues.push_back(diagnostic(
                ErrorCode::invalid_topology,
                "A boundary's half-edge references are inconsistent."));
            return issues;
        }

        const ParameterPoint from = vertices_.at(boundary.from);
        const ParameterPoint to = vertices_.at(boundary.to);
        if (distance(from, to) <= tolerance.linear) {
            issues.push_back(diagnostic(
                ErrorCode::invalid_topology,
                "A boundary has zero or near-zero length."));
            return issues;
        }

        const bool outer = lies_on_domain_boundary(
            from, to, domain_, tolerance.linear);
        if (outer != (expected_forward + expected_reverse == 1U)) {
            issues.push_back(diagnostic(
                ErrorCode::invalid_topology,
                "Boundary incidence does not match the domain perimeter."));
            return issues;
        }
    }

    double total_area = 0.0;
    std::map<RegionId, std::vector<ParameterPoint>> polygons;
    for (const auto& [region_id, region] : regions_) {
        std::vector<ParameterPoint> points =
            region_points(region, boundaries_, vertices_);
        const double area = signed_area(points);
        if (!std::isfinite(area) || area <= tolerance.area) {
            issues.push_back(diagnostic(
                ErrorCode::invalid_topology,
                "A region has zero area or is not oriented counterclockwise."));
            return issues;
        }
        total_area += area;
        polygons.emplace(region_id, std::move(points));
    }

    for (auto first = boundaries_.begin(); first != boundaries_.end(); ++first) {
        const Boundary& a = first->second;
        const ParameterPoint a0 = vertices_.at(a.from);
        const ParameterPoint a1 = vertices_.at(a.to);
        auto second = std::next(first);
        for (; second != boundaries_.end(); ++second) {
            const Boundary& b = second->second;
            const ParameterPoint b0 = vertices_.at(b.from);
            const ParameterPoint b1 = vertices_.at(b.to);
            if (!segments_intersect(
                    a0, a1, b0, b1, tolerance.linear, tolerance.cross)) {
                continue;
            }

            const bool share_from = a.from == b.from || a.from == b.to;
            const bool share_to = a.to == b.from || a.to == b.to;
            if (share_from != share_to) {
                const VertexId shared = share_from ? a.from : a.to;
                const VertexId a_other = a.from == shared ? a.to : a.from;
                const VertexId b_other = b.from == shared ? b.to : b.from;
                const ParameterPoint origin = vertices_.at(shared);
                const ParameterPoint a_point = vertices_.at(a_other);
                const ParameterPoint b_point = vertices_.at(b_other);
                const double a_length = distance(a_point, origin);
                const double b_length = distance(b_point, origin);
                const double dot =
                    ((a_point.u - origin.u) / a_length)
                        * ((b_point.u - origin.u) / b_length)
                    + ((a_point.v - origin.v) / a_length)
                        * ((b_point.v - origin.v) / b_length);
                if (std::abs(cross(origin, a_point, b_point)) <= tolerance.cross
                    && dot > 0.0) {
                    issues.push_back(diagnostic(
                        ErrorCode::invalid_topology,
                        "Two boundaries overlap beyond a shared vertex."));
                    return issues;
                }
                continue;
            }

            issues.push_back(diagnostic(
                ErrorCode::invalid_topology,
                "Boundaries cross, overlap, or meet without a shared vertex."));
            return issues;
        }
    }

    for (auto first = polygons.begin(); first != polygons.end(); ++first) {
        auto second = std::next(first);
        for (; second != polygons.end(); ++second) {
            const bool first_inside = std::any_of(
                first->second.begin(),
                first->second.end(),
                [&](ParameterPoint point) {
                    return point_in_polygon(
                        point, second->second, tolerance.linear, tolerance.cross);
                });
            const bool second_inside = std::any_of(
                second->second.begin(),
                second->second.end(),
                [&](ParameterPoint point) {
                    return point_in_polygon(
                        point, first->second, tolerance.linear, tolerance.cross);
                });
            if (first_inside || second_inside) {
                issues.push_back(diagnostic(
                    ErrorCode::invalid_topology,
                    "Region interiors overlap."));
                return issues;
            }
        }
    }

    const double domain_area = domain_.width() * domain_.height();
    const double coverage_tolerance =
        std::max(domain_area * 1e-10, tolerance.area * 4.0);
    if (std::abs(total_area - domain_area) > coverage_tolerance) {
        issues.push_back(diagnostic(
            ErrorCode::invalid_topology,
            "Region areas do not cover the complete parameter domain."));
    }

    return issues;
}

Result<VertexId> Partition::move_vertex(
    VertexId vertex,
    ParameterPoint position)
{
    Partition candidate = *this;
    auto result = candidate.move_vertex_in_place(vertex, position);
    if (!result) {
        return Result<VertexId>::failure(result.error());
    }

    const std::vector<Diagnostic> issues = candidate.validate();
    if (!issues.empty()) {
        return Result<VertexId>::failure(issues.front());
    }

    *this = std::move(candidate);
    return Result<VertexId>::success(result.value());
}

Result<VertexId> Partition::move_vertex_in_place(
    VertexId vertex,
    ParameterPoint position)
{
    auto found = vertices_.find(vertex);
    if (found == vertices_.end()) {
        return Result<VertexId>::failure(diagnostic(
            ErrorCode::not_found,
            "Cannot move a missing vertex."));
    }
    if (!domain_.contains(position)) {
        return Result<VertexId>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Moved vertex must remain inside the parameter domain."));
    }

    const bool on_domain_boundary = std::any_of(
        boundaries_.begin(),
        boundaries_.end(),
        [vertex](const auto& item) {
            const Boundary& boundary = item.second;
            const bool incident = boundary.from == vertex || boundary.to == vertex;
            const bool exterior =
                static_cast<bool>(boundary.left_region)
                != static_cast<bool>(boundary.right_region);
            return incident && exterior;
        });
    if (on_domain_boundary) {
        return Result<VertexId>::failure(diagnostic(
            ErrorCode::unsupported_operation,
            "Domain-perimeter vertices are fixed in the initial partition model."));
    }

    found->second = position;
    return Result<VertexId>::success(vertex);
}

Result<VertexId> Partition::insert_boundary_vertex(
    BoundaryId boundary,
    ParameterPoint position)
{
    Partition candidate = *this;
    auto result = candidate.insert_boundary_vertex_in_place(boundary, position);
    if (!result) {
        return Result<VertexId>::failure(result.error());
    }

    const std::vector<Diagnostic> issues = candidate.validate();
    if (!issues.empty()) {
        return Result<VertexId>::failure(issues.front());
    }

    *this = std::move(candidate);
    return Result<VertexId>::success(result.value());
}

Result<VertexId> Partition::insert_boundary_vertex_in_place(
    BoundaryId boundary_id,
    ParameterPoint position)
{
    const auto found = boundaries_.find(boundary_id);
    if (found == boundaries_.end()) {
        return Result<VertexId>::failure(diagnostic(
            ErrorCode::not_found,
            "Cannot insert a vertex on a missing boundary."));
    }
    if (!domain_.contains(position)) {
        return Result<VertexId>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Inserted vertex must lie inside the parameter domain."));
    }

    const GeometryTolerance tolerance = geometry_tolerance(domain_);
    for (const auto& [vertex_id, point] : vertices_) {
        (void)vertex_id;
        if (distance(position, point) <= tolerance.linear) {
            return Result<VertexId>::failure(diagnostic(
                ErrorCode::invalid_argument,
                "Inserted vertex duplicates an existing vertex."));
        }
    }

    const Boundary old_boundary = found->second;
    const ParameterPoint from = vertices_.at(old_boundary.from);
    const ParameterPoint to = vertices_.at(old_boundary.to);
    if (!point_on_segment(
            position, from, to, tolerance.linear, tolerance.cross)) {
        return Result<VertexId>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Inserted vertex must lie on the selected boundary."));
    }

    const double delta_u = to.u - from.u;
    const double delta_v = to.v - from.v;
    const double segment_length = std::hypot(delta_u, delta_v);
    const double fraction =
        ((position.u - from.u) / segment_length) * (delta_u / segment_length)
        + ((position.v - from.v) / segment_length) * (delta_v / segment_length);
    const double fraction_tolerance = std::max(
        tolerance.linear / segment_length,
        std::numeric_limits<double>::epsilon() * 8.0);
    if (fraction <= fraction_tolerance
        || fraction >= 1.0 - fraction_tolerance) {
        return Result<VertexId>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Inserted vertex must be strictly inside the boundary segment."));
    }

    const VertexId inserted{next_vertex_id_++};
    const BoundaryId first_id{next_boundary_id_++};
    const BoundaryId second_id{next_boundary_id_++};
    vertices_.emplace(inserted, position);
    boundaries_.erase(boundary_id);
    boundaries_.emplace(
        first_id,
        Boundary{first_id, old_boundary.from, inserted,
                 old_boundary.left_region, old_boundary.right_region});
    boundaries_.emplace(
        second_id,
        Boundary{second_id, inserted, old_boundary.to,
                 old_boundary.left_region, old_boundary.right_region});

    for (auto& [region_id, region] : regions_) {
        (void)region_id;
        std::vector<BoundaryUse> updated;
        updated.reserve(region.boundary.size() + 1);
        for (const BoundaryUse& use : region.boundary) {
            if (use.boundary != boundary_id) {
                updated.push_back(use);
            } else if (use.forward) {
                updated.push_back({first_id, true});
                updated.push_back({second_id, true});
            } else {
                updated.push_back({second_id, false});
                updated.push_back({first_id, false});
            }
        }
        region.boundary = std::move(updated);
    }

    return Result<VertexId>::success(inserted);
}

Result<RegionId> Partition::split_region(
    RegionId region,
    VertexId first,
    VertexId second)
{
    Partition candidate = *this;
    auto result = candidate.split_region_in_place(region, first, second);
    if (!result) {
        return Result<RegionId>::failure(result.error());
    }

    const std::vector<Diagnostic> issues = candidate.validate();
    if (!issues.empty()) {
        return Result<RegionId>::failure(issues.front());
    }

    *this = std::move(candidate);
    return Result<RegionId>::success(result.value());
}

Result<RegionId> Partition::split_region_in_place(
    RegionId region_id,
    VertexId first,
    VertexId second)
{
    auto region_it = regions_.find(region_id);
    if (region_it == regions_.end()) {
        return Result<RegionId>::failure(diagnostic(
            ErrorCode::not_found,
            "Cannot split a missing region."));
    }
    if (first == second || !vertices_.contains(first) || !vertices_.contains(second)) {
        return Result<RegionId>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "A split requires two distinct existing boundary vertices."));
    }

    Region& region = region_it->second;
    std::size_t first_index = region.boundary.size();
    std::size_t second_index = region.boundary.size();
    for (std::size_t index = 0; index < region.boundary.size(); ++index) {
        const BoundaryUse& use = region.boundary[index];
        const Boundary& boundary = boundaries_.at(use.boundary);
        const VertexId start = use_start_vertex(boundary, use);
        if (start == first) {
            first_index = index;
        }
        if (start == second) {
            second_index = index;
        }
    }
    if (first_index == region.boundary.size()
        || second_index == region.boundary.size()) {
        return Result<RegionId>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Split vertices must both lie on the selected region boundary."));
    }

    const std::size_t count = region.boundary.size();
    if ((first_index + 1) % count == second_index
        || (second_index + 1) % count == first_index) {
        return Result<RegionId>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Adjacent vertices cannot split a region."));
    }

    const ParameterPoint first_point = vertices_.at(first);
    const ParameterPoint second_point = vertices_.at(second);
    const GeometryTolerance tolerance = geometry_tolerance(domain_);
    for (const BoundaryUse& use : region.boundary) {
        const Boundary& boundary = boundaries_.at(use.boundary);
        const bool incident = boundary.from == first || boundary.to == first
            || boundary.from == second || boundary.to == second;
        if (!incident && segments_intersect(
                first_point,
                second_point,
                vertices_.at(boundary.from),
                vertices_.at(boundary.to),
                tolerance.linear,
                tolerance.cross)) {
            return Result<RegionId>::failure(diagnostic(
                ErrorCode::invalid_argument,
                "The split segment crosses the region boundary."));
        }
    }

    const std::vector<ParameterPoint> old_points =
        region_points(region, boundaries_, vertices_);
    const ParameterPoint midpoint{
        std::midpoint(first_point.u, second_point.u),
        std::midpoint(first_point.v, second_point.v)
    };
    if (!point_in_polygon(
            midpoint, old_points, tolerance.linear, tolerance.cross)) {
        return Result<RegionId>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "The split segment must lie inside the region."));
    }

    std::vector<BoundaryUse> first_path;
    for (std::size_t index = first_index; index != second_index;
         index = (index + 1) % count) {
        first_path.push_back(region.boundary[index]);
    }
    std::vector<BoundaryUse> second_path;
    for (std::size_t index = second_index; index != first_index;
         index = (index + 1) % count) {
        second_path.push_back(region.boundary[index]);
    }

    const RegionId new_region{next_region_id_++};
    const BoundaryId dividing_boundary{next_boundary_id_++};
    for (const BoundaryUse& use : second_path) {
        Boundary& boundary = boundaries_.at(use.boundary);
        if (boundary.left_region == region_id) {
            boundary.left_region = new_region;
        } else if (boundary.right_region == region_id) {
            boundary.right_region = new_region;
        }
    }

    boundaries_.emplace(
        dividing_boundary,
        Boundary{
            dividing_boundary,
            first,
            second,
            new_region,
            region_id
        });

    const bool background = region.background;
    region.boundary = std::move(first_path);
    region.boundary.push_back({dividing_boundary, false});
    regions_.emplace(
        new_region,
        Region{
            new_region,
            std::move(second_path),
            background
        });
    regions_.at(new_region).boundary.push_back({dividing_boundary, true});
    return Result<RegionId>::success(new_region);
}

Result<RegionId> Partition::merge_regions(
    RegionId retained,
    RegionId removed)
{
    Partition candidate = *this;
    auto result = candidate.merge_regions_in_place(retained, removed);
    if (!result) {
        return Result<RegionId>::failure(result.error());
    }

    const std::vector<Diagnostic> issues = candidate.validate();
    if (!issues.empty()) {
        return Result<RegionId>::failure(issues.front());
    }

    *this = std::move(candidate);
    return Result<RegionId>::success(result.value());
}

Result<RegionId> Partition::merge_regions_in_place(
    RegionId retained,
    RegionId removed)
{
    if (retained == removed) {
        return Result<RegionId>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "A region cannot be merged with itself."));
    }

    auto retained_it = regions_.find(retained);
    auto removed_it = regions_.find(removed);
    if (retained_it == regions_.end() || removed_it == regions_.end()) {
        return Result<RegionId>::failure(diagnostic(
            ErrorCode::not_found,
            "Both regions must exist before they can be merged."));
    }
    if (retained_it->second.background != removed_it->second.background) {
        return Result<RegionId>::failure(diagnostic(
            ErrorCode::unsupported_operation,
            "Active and background regions cannot be merged together."));
    }

    std::set<BoundaryId> shared;
    for (const auto& [boundary_id, boundary] : boundaries_) {
        if ((boundary.left_region == retained && boundary.right_region == removed)
            || (boundary.left_region == removed && boundary.right_region == retained)) {
            shared.insert(boundary_id);
        }
    }
    if (shared.empty()) {
        return Result<RegionId>::failure(diagnostic(
            ErrorCode::unsupported_operation,
            "Regions must share a boundary before they can be merged."));
    }

    const auto remaining_path = [&](const Region& region)
        -> std::optional<std::vector<BoundaryUse>> {
        const std::size_t count = region.boundary.size();
        std::size_t transitions = 0;
        std::size_t start = count;
        std::size_t shared_count = 0;
        for (std::size_t index = 0; index < count; ++index) {
            const bool current_shared =
                shared.contains(region.boundary[index].boundary);
            const bool previous_shared =
                shared.contains(region.boundary[(index + count - 1) % count].boundary);
            if (current_shared) {
                ++shared_count;
            }
            if (current_shared != previous_shared) {
                ++transitions;
            }
            if (!current_shared && previous_shared) {
                start = index;
            }
        }

        if (shared_count == 0 || transitions != 2 || start == count) {
            return std::nullopt;
        }

        std::vector<BoundaryUse> path;
        path.reserve(count - shared_count);
        for (std::size_t offset = 0; offset < count - shared_count; ++offset) {
            const BoundaryUse& use = region.boundary[(start + offset) % count];
            if (shared.contains(use.boundary)) {
                return std::nullopt;
            }
            path.push_back(use);
        }
        return path;
    };

    std::optional<std::vector<BoundaryUse>> retained_path =
        remaining_path(retained_it->second);
    std::optional<std::vector<BoundaryUse>> removed_path =
        remaining_path(removed_it->second);
    if (!retained_path || !removed_path
        || retained_path->empty() || removed_path->empty()) {
        return Result<RegionId>::failure(diagnostic(
            ErrorCode::unsupported_operation,
            "The shared interface must be one contiguous boundary chain."));
    }

    const BoundaryUse& retained_last = retained_path->back();
    const Boundary& retained_last_boundary =
        boundaries_.at(retained_last.boundary);
    const BoundaryUse& removed_first = removed_path->front();
    const Boundary& removed_first_boundary =
        boundaries_.at(removed_first.boundary);
    if (use_end_vertex(retained_last_boundary, retained_last)
        != use_start_vertex(removed_first_boundary, removed_first)) {
        return Result<RegionId>::failure(diagnostic(
            ErrorCode::invalid_topology,
            "The incident region boundary chains do not meet."));
    }

    std::vector<BoundaryUse> merged = std::move(*retained_path);
    merged.insert(merged.end(), removed_path->begin(), removed_path->end());

    for (auto& [boundary_id, boundary] : boundaries_) {
        if (shared.contains(boundary_id)) {
            continue;
        }
        if (boundary.left_region == removed) {
            boundary.left_region = retained;
        }
        if (boundary.right_region == removed) {
            boundary.right_region = retained;
        }
    }

    for (BoundaryId boundary_id : shared) {
        boundaries_.erase(boundary_id);
    }
    regions_.erase(removed);
    regions_.at(retained).boundary = std::move(merged);
    return Result<RegionId>::success(retained);
}

Result<RegionId> Partition::delete_region(RegionId region)
{
    Partition candidate = *this;
    auto result = candidate.delete_region_in_place(region);
    if (!result) {
        return Result<RegionId>::failure(result.error());
    }

    const std::vector<Diagnostic> issues = candidate.validate();
    if (!issues.empty()) {
        return Result<RegionId>::failure(issues.front());
    }

    *this = std::move(candidate);
    return Result<RegionId>::success(result.value());
}

Result<RegionId> Partition::delete_region_in_place(RegionId region_id)
{
    auto found = regions_.find(region_id);
    if (found == regions_.end()) {
        return Result<RegionId>::failure(diagnostic(
            ErrorCode::not_found,
            "Cannot delete a missing region."));
    }
    if (found->second.background) {
        return Result<RegionId>::failure(diagnostic(
            ErrorCode::unsupported_operation,
            "A background region cannot be deleted while full coverage is required."));
    }

    const RegionId background_id{next_region_id_++};
    Region background{background_id, found->second.boundary, true};
    for (const BoundaryUse& use : background.boundary) {
        Boundary& boundary = boundaries_.at(use.boundary);
        if (boundary.left_region == region_id) {
            boundary.left_region = background_id;
        } else if (boundary.right_region == region_id) {
            boundary.right_region = background_id;
        }
    }

    regions_.erase(found);
    regions_.emplace(background_id, std::move(background));
    return Result<RegionId>::success(background_id);
}

}  // namespace pattern
