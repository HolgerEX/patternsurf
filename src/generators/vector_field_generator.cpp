#include "pattern/generators/vector_field_generator.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <optional>
#include <string>
#include <utility>

namespace pattern {
namespace {

constexpr double max_curve_steps = 100000.0;
constexpr double max_seed_count = 250000.0;
constexpr double max_total_steps = 10000000.0;

Diagnostic diagnostic(ErrorCode code, std::string message)
{
    return {code, Severity::error, std::move(message)};
}

VertexId start_vertex(const Boundary& boundary, const BoundaryUse& use)
{
    return use.forward ? boundary.from : boundary.to;
}

std::vector<ParameterPoint> polygon_for_region(
    const Partition& partition,
    const Region& region)
{
    std::vector<ParameterPoint> polygon;
    polygon.reserve(region.boundary.size());
    for (const BoundaryUse& use : region.boundary) {
        const Boundary& boundary = partition.boundaries().at(use.boundary);
        polygon.push_back(partition.vertices().at(start_vertex(boundary, use)));
    }
    return polygon;
}

double cross(ParameterPoint first, ParameterPoint second)
{
    return first.u * second.v - first.v * second.u;
}

double cross(
    ParameterPoint origin,
    ParameterPoint first,
    ParameterPoint second)
{
    return cross(
        {first.u - origin.u, first.v - origin.v},
        {second.u - origin.u, second.v - origin.v});
}

double distance(ParameterPoint first, ParameterPoint second)
{
    return std::hypot(first.u - second.u, first.v - second.v);
}

bool point_on_segment(
    ParameterPoint point,
    ParameterPoint first,
    ParameterPoint second,
    double linear_tolerance,
    double cross_tolerance)
{
    return std::abs(cross(first, second, point)) <= cross_tolerance
        && point.u >= std::min(first.u, second.u) - linear_tolerance
        && point.u <= std::max(first.u, second.u) + linear_tolerance
        && point.v >= std::min(first.v, second.v) - linear_tolerance
        && point.v <= std::max(first.v, second.v) + linear_tolerance;
}

bool contains(
    ParameterPoint point,
    const std::vector<ParameterPoint>& polygon,
    double linear_tolerance,
    double cross_tolerance)
{
    bool inside = false;
    for (std::size_t index = 0; index < polygon.size(); ++index) {
        const ParameterPoint first = polygon[index];
        const ParameterPoint second = polygon[(index + 1) % polygon.size()];
        if (point_on_segment(
                point, first, second, linear_tolerance, cross_tolerance)) {
            return true;
        }
        if ((first.v > point.v) != (second.v > point.v)) {
            const double intersection_u = first.u
                + (point.v - first.v) * (second.u - first.u)
                    / (second.v - first.v);
            if (point.u < intersection_u) {
                inside = !inside;
            }
        }
    }
    return inside;
}

std::optional<ParameterPoint> first_boundary_intersection(
    ParameterPoint from,
    ParameterPoint to,
    const std::vector<ParameterPoint>& polygon,
    double tolerance)
{
    const ParameterPoint direction{to.u - from.u, to.v - from.v};
    double earliest = std::numeric_limits<double>::infinity();

    for (std::size_t index = 0; index < polygon.size(); ++index) {
        const ParameterPoint edge_start = polygon[index];
        const ParameterPoint edge_end = polygon[(index + 1) % polygon.size()];
        const ParameterPoint edge{
            edge_end.u - edge_start.u,
            edge_end.v - edge_start.v
        };
        const double step_length = distance(from, to);
        const double edge_length = distance(edge_start, edge_end);
        const double denominator = cross(direction, edge);
        const double denominator_tolerance =
            tolerance * std::max(step_length, edge_length);
        if (std::abs(denominator) <= denominator_tolerance) {
            continue;
        }

        const double step_parameter_tolerance = tolerance / step_length;
        const double edge_parameter_tolerance = tolerance / edge_length;
        const ParameterPoint offset{
            edge_start.u - from.u,
            edge_start.v - from.v
        };
        const double along_step = cross(offset, edge) / denominator;
        const double along_edge = cross(offset, direction) / denominator;
        if (along_step < -step_parameter_tolerance
            || along_step > 1.0 + step_parameter_tolerance
            || along_edge < -edge_parameter_tolerance
            || along_edge > 1.0 + edge_parameter_tolerance) {
            continue;
        }

        const double clamped = std::clamp(along_step, 0.0, 1.0);
        earliest = std::min(earliest, clamped);
    }

    if (!std::isfinite(earliest)) {
        return std::nullopt;
    }

    return ParameterPoint{
        from.u + earliest * direction.u,
        from.v + earliest * direction.v
    };
}

Result<Vector2> direction(
    ParameterPoint point,
    const VectorFieldConfig& config,
    const CoordinateMap2D& coordinate_map,
    double sign)
{
    const auto mapped = coordinate_map.evaluate(point);
    if (!mapped) {
        return Result<Vector2>::failure(mapped.error());
    }
    const auto jacobian = coordinate_map.jacobian(point);
    if (!jacobian) {
        return Result<Vector2>::failure(jacobian.error());
    }

    const double phase_h = config.warp_frequency * mapped.value().h;
    const double phase_k = config.warp_frequency * mapped.value().k;
    if (!std::isfinite(phase_h) || !std::isfinite(phase_k)) {
        return Result<Vector2>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Warp frequency overflows in derived coordinates."));
    }
    const double angle = config.base_angle
        + config.warp_amplitude
            * std::sin(phase_h)
            * std::sin(phase_k);
    if (!std::isfinite(angle)) {
        return Result<Vector2>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Vector-field angle is non-finite."));
    }

    const Jacobian2& matrix = jacobian.value();
    const double scale = std::max({
        std::abs(matrix.dh_du),
        std::abs(matrix.dh_dv),
        std::abs(matrix.dk_du),
        std::abs(matrix.dk_dv)
    });
    if (!std::isfinite(scale) || scale == 0.0) {
        return Result<Vector2>::failure(diagnostic(
            ErrorCode::unsupported_operation,
            "Vector field requires a locally invertible coordinate map."));
    }

    const double a = matrix.dh_du / scale;
    const double b = matrix.dh_dv / scale;
    const double c = matrix.dk_du / scale;
    const double d = matrix.dk_dv / scale;
    const double determinant = a * d - b * c;
    if (!std::isfinite(determinant)
        || std::abs(determinant) <= std::numeric_limits<double>::epsilon() * 32.0) {
        return Result<Vector2>::failure(diagnostic(
            ErrorCode::unsupported_operation,
            "Vector field requires a locally invertible coordinate map."));
    }

    const double field_h = std::cos(angle) * sign;
    const double field_k = std::sin(angle) * sign;
    const double vector_u = ((d * field_h - b * field_k) / determinant) / scale;
    const double vector_v = ((a * field_k - c * field_h) / determinant) / scale;
    const double magnitude = std::hypot(vector_u, vector_v);
    if (!std::isfinite(magnitude) || magnitude == 0.0) {
        return Result<Vector2>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Coordinate-map pullback produced an invalid direction."));
    }
    return Result<Vector2>::success({
        vector_u / magnitude,
        vector_v / magnitude
    });
}

Result<ParameterPoint> advance(
    ParameterPoint point,
    double step,
    const VectorFieldConfig& config,
    const CoordinateMap2D& coordinate_map,
    double sign)
{
    const auto initial = direction(point, config, coordinate_map, sign);
    if (!initial) {
        return Result<ParameterPoint>::failure(initial.error());
    }
    const ParameterPoint midpoint{
        point.u + 0.5 * step * initial.value().x,
        point.v + 0.5 * step * initial.value().y
    };
    const auto middle = direction(midpoint, config, coordinate_map, sign);
    if (!middle) {
        return Result<ParameterPoint>::failure(middle.error());
    }
    return Result<ParameterPoint>::success({
        point.u + step * middle.value().x,
        point.v + step * middle.value().y
    });
}

Result<std::vector<ParameterPoint>> trace(
    ParameterPoint seed,
    double sign,
    const VectorFieldConfig& config,
    const CoordinateMap2D& coordinate_map,
    const std::vector<ParameterPoint>& polygon,
    double linear_tolerance,
    double cross_tolerance,
    double branch_length,
    std::size_t step_count)
{
    std::vector<ParameterPoint> points{seed};
    ParameterPoint current = seed;
    double travelled = 0.0;

    for (std::size_t index = 0; index < step_count
         && travelled < branch_length; ++index) {
        const double step = std::min(
            config.step_size,
            branch_length - travelled);
        const auto advanced = advance(
            current, step, config, coordinate_map, sign);
        if (!advanced) {
            return Result<std::vector<ParameterPoint>>::failure(advanced.error());
        }
        const ParameterPoint next = advanced.value();
        if (!std::isfinite(next.u) || !std::isfinite(next.v)) {
            return Result<std::vector<ParameterPoint>>::failure(diagnostic(
                ErrorCode::invalid_argument,
                "The vector field produced a non-finite integration point."));
        }
        const auto intersection = first_boundary_intersection(
            current, next, polygon, linear_tolerance);
        if (intersection) {
            if (distance(current, *intersection) > linear_tolerance) {
                points.push_back(*intersection);
            }
            break;
        }
        if (!contains(next, polygon, linear_tolerance, cross_tolerance)) {
            return Result<std::vector<ParameterPoint>>::failure(diagnostic(
                ErrorCode::invalid_topology,
                "A curve step left the region without intersecting its boundary."));
        }

        travelled += distance(current, next);
        current = next;
        points.push_back(current);
    }

    return Result<std::vector<ParameterPoint>>::success(std::move(points));
}

}  // namespace

VectorFieldCurveGenerator::VectorFieldCurveGenerator(
    VectorFieldConfig config,
    std::shared_ptr<const CoordinateMap2D> coordinate_map)
    : config_(config),
      coordinate_map_(std::move(coordinate_map))
{
}

Result<VectorFieldCurveGenerator> VectorFieldCurveGenerator::create(
    VectorFieldConfig config,
    std::shared_ptr<const CoordinateMap2D> coordinate_map)
{
    if (!std::isfinite(config.seed_spacing) || config.seed_spacing <= 0.0
        || !std::isfinite(config.step_size) || config.step_size <= 0.0
        || !std::isfinite(config.maximum_length) || config.maximum_length <= 0.0
        || !std::isfinite(config.base_angle)
        || !std::isfinite(config.warp_amplitude)
        || !std::isfinite(config.warp_frequency)) {
        return Result<VectorFieldCurveGenerator>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Vector-field parameters must be finite and spacing, step size, and length positive."));
    }

    const double steps = std::ceil(
        (0.5 * config.maximum_length) / config.step_size);
    if (!std::isfinite(steps) || steps > max_curve_steps) {
        return Result<VectorFieldCurveGenerator>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Requested curve length and step size exceed the supported step limit."));
    }

    if (!coordinate_map) {
        coordinate_map = std::make_shared<IdentityCoordinateMap>();
    }
    const CoordinateMapCapabilities capabilities = coordinate_map->capabilities();
    if (!capabilities.has_jacobian || !capabilities.has_inverse) {
        return Result<VectorFieldCurveGenerator>::failure(diagnostic(
            ErrorCode::unsupported_operation,
            "Vector-field curves require a coordinate map with a Jacobian and inverse."));
    }

    return Result<VectorFieldCurveGenerator>::success(
        VectorFieldCurveGenerator(config, std::move(coordinate_map)));
}

Result<CurveSet> VectorFieldCurveGenerator::generate(
    const Partition& partition,
    RegionId region_id) const
{
    const Region* region = partition.find_region(region_id);
    if (region == nullptr) {
        return Result<CurveSet>::failure(diagnostic(
            ErrorCode::not_found,
            "Cannot generate curves for a missing region."));
    }
    if (region->background) {
        return Result<CurveSet>::failure(diagnostic(
            ErrorCode::unsupported_operation,
            "Pattern output cannot be generated for a background region."));
    }

    const std::vector<ParameterPoint> polygon =
        polygon_for_region(partition, *region);
    const Domain2& domain = partition.domain();
    double min_u = domain.max_u();
    double min_v = domain.max_v();
    double max_u = domain.min_u();
    double max_v = domain.min_v();
    for (ParameterPoint point : polygon) {
        min_u = std::min(min_u, point.u);
        min_v = std::min(min_v, point.v);
        max_u = std::max(max_u, point.u);
        max_v = std::max(max_v, point.v);
    }

    const double columns = std::max(
        1.0, std::ceil((max_u - min_u) / config_.seed_spacing));
    const double rows = std::max(
        1.0, std::ceil((max_v - min_v) / config_.seed_spacing));
    if (!std::isfinite(columns) || !std::isfinite(rows)
        || columns > max_seed_count
        || rows > max_seed_count
        || columns * rows > max_seed_count) {
        return Result<CurveSet>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Seed spacing would create more than the supported number of samples."));
    }

    const std::size_t column_count = static_cast<std::size_t>(columns);
    const std::size_t row_count = static_cast<std::size_t>(rows);
    const double cell_width = (max_u - min_u) / columns;
    const double cell_height = (max_v - min_v) / rows;
    const double scale = std::max(domain.width(), domain.height());
    const double linear_tolerance = scale * 1e-12;
    const double cross_tolerance = linear_tolerance * scale;
    if (config_.step_size <= linear_tolerance) {
        return Result<CurveSet>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Integration step size is below the geometric tolerance."));
    }
    const double branch_length = 0.5 * config_.maximum_length;
    const std::size_t step_count = static_cast<std::size_t>(
        std::ceil(branch_length / config_.step_size));
    if (columns * rows * static_cast<double>(step_count) * 2.0
        > max_total_steps) {
        return Result<CurveSet>::failure(diagnostic(
            ErrorCode::invalid_argument,
            "Requested seed density and integration length exceed the work limit."));
    }

    CurveSet output;
    for (std::size_t row = 0; row < row_count; ++row) {
        for (std::size_t column = 0; column < column_count; ++column) {
            const ParameterPoint seed{
                min_u + (static_cast<double>(column) + 0.5) * cell_width,
                min_v + (static_cast<double>(row) + 0.5) * cell_height
            };
            if (!contains(seed, polygon, linear_tolerance, cross_tolerance)) {
                continue;
            }

            auto backward = trace(
                seed,
                -1.0,
                config_,
                *coordinate_map_,
                polygon,
                linear_tolerance,
                cross_tolerance,
                branch_length,
                step_count);
            if (!backward) {
                return Result<CurveSet>::failure(backward.error());
            }
            auto forward = trace(
                seed,
                1.0,
                config_,
                *coordinate_map_,
                polygon,
                linear_tolerance,
                cross_tolerance,
                branch_length,
                step_count);
            if (!forward) {
                return Result<CurveSet>::failure(forward.error());
            }

            Polyline curve;
            curve.points.reserve(
                backward.value().size() + forward.value().size() - 1);
            for (auto point = backward.value().rbegin();
                 point != backward.value().rend(); ++point) {
                curve.points.push_back(*point);
            }
            curve.points.insert(
                curve.points.end(),
                std::next(forward.value().begin()),
                forward.value().end());
            if (curve.points.size() >= 2) {
                output.curves.push_back(std::move(curve));
            }
        }
    }

    return Result<CurveSet>::success(std::move(output));
}

}  // namespace pattern
