#include "pattern/generators/vector_field_generator.hpp"

#include <cstddef>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string_view>

namespace {

int failures = 0;

void expect(bool condition, std::string_view message)
{
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

pattern::Partition make_partition()
{
    const auto domain = pattern::Domain2::rectangle(0.0, 0.0, 1.0, 1.0);
    if (!domain) {
        std::cerr << "Test setup failed: " << domain.error().message << '\n';
        std::exit(EXIT_FAILURE);
    }
    return pattern::Partition::create(domain.value());
}

void test_generator_descriptor_and_parameters()
{
    constexpr auto descriptor = pattern::VectorFieldCurveGenerator::descriptor();
    expect(descriptor.type_id == "pattern.vector_field.curves",
        "generator has a stable serialized type ID");
    expect(descriptor.schema_version == 1,
        "generator exposes a schema version");
    expect(descriptor.deterministic, "generator declares deterministic output");

    pattern::VectorFieldConfig invalid;
    invalid.step_size = 0.0;
    expect(!pattern::VectorFieldCurveGenerator::create(invalid),
        "zero integration step is rejected");

    invalid = {};
    invalid.maximum_length = std::numeric_limits<double>::infinity();
    expect(!pattern::VectorFieldCurveGenerator::create(invalid),
        "non-finite maximum curve length is rejected");
}

void test_curves_are_deterministic_and_clipped()
{
    pattern::VectorFieldConfig config;
    config.seed_spacing = 0.3;
    config.step_size = 0.05;
    config.maximum_length = 2.0;
    config.base_angle = 0.0;
    const auto generator = pattern::VectorFieldCurveGenerator::create(config);
    expect(generator.has_value(), "valid vector-field generator is created");
    if (!generator) {
        return;
    }

    const pattern::Partition partition = make_partition();
    const auto first = generator.value().generate(partition, pattern::RegionId{1});
    const auto second = generator.value().generate(partition, pattern::RegionId{1});
    expect(first.has_value() && second.has_value(),
        "generator produces curves for a valid region");
    if (!first || !second) {
        return;
    }

    expect(!first.value().curves.empty(), "seeds produce at least one curve");
    expect(first.value().curves.size() == second.value().curves.size(),
        "identical inputs produce the same number of curves");
    if (first.value().curves.size() != second.value().curves.size()) {
        return;
    }

    for (std::size_t curve_index = 0;
         curve_index < first.value().curves.size();
         ++curve_index) {
        const auto& curve = first.value().curves[curve_index];
        const auto& repeated = second.value().curves[curve_index];
        expect(curve.points.size() == repeated.points.size(),
            "identical inputs produce the same curve lengths");
        if (curve.points.size() != repeated.points.size()) {
            continue;
        }

        expect(curve.points.size() >= 2, "each output curve has at least two points");
        expect(curve.points.front().u <= 1e-10
                && curve.points.back().u >= 1.0 - 1e-10,
            "constant horizontal flow is clipped to both region boundaries");
        const double total_length =
            curve.points.back().u - curve.points.front().u;
        expect(total_length <= config.maximum_length + 1e-10,
            "combined curve length respects its configured maximum");

        for (std::size_t point_index = 0;
             point_index < curve.points.size();
             ++point_index) {
            const auto point = curve.points[point_index];
            const auto repeated_point = repeated.points[point_index];
            expect(std::abs(point.u - repeated_point.u) <= 1e-14
                    && std::abs(point.v - repeated_point.v) <= 1e-14,
                "identical inputs reproduce curve coordinates");
            expect(point.u >= -1e-10 && point.u <= 1.0 + 1e-10
                    && point.v >= -1e-10 && point.v <= 1.0 + 1e-10,
                "curve samples stay within the region");
        }
    }
}

void test_missing_region_and_seed_limit()
{
    const auto generator = pattern::VectorFieldCurveGenerator::create({});
    expect(generator.has_value(), "default generator configuration is valid");
    if (!generator) {
        return;
    }

    const pattern::Partition partition = make_partition();
    expect(!generator.value().generate(partition, pattern::RegionId{99}),
        "missing region is reported");

    pattern::VectorFieldConfig dense;
    dense.seed_spacing = 1e-4;
    const auto dense_generator = pattern::VectorFieldCurveGenerator::create(dense);
    expect(dense_generator.has_value(), "dense generator parameters are valid");
    if (dense_generator) {
        expect(!dense_generator.value().generate(partition, pattern::RegionId{1}),
            "excessive seed counts are rejected rather than allocated");
    }

    pattern::VectorFieldConfig costly;
    costly.seed_spacing = 0.01;
    costly.step_size = 0.0001;
    const auto costly_generator = pattern::VectorFieldCurveGenerator::create(costly);
    expect(costly_generator.has_value(), "costly but valid parameters are accepted");
    if (costly_generator) {
        expect(!costly_generator.value().generate(partition, pattern::RegionId{1}),
            "excessive total integration work is rejected");
    }
}

void test_curves_clip_to_polygonal_region()
{
    const auto domain = pattern::Domain2::rectangle(0.0, 0.0, 2.0, 1.0);
    if (!domain) {
        expect(false, "polygon clipping test domain is valid");
        return;
    }

    pattern::Partition partition = pattern::Partition::create(domain.value());
    const auto split = partition.split_region(
        pattern::RegionId{1}, pattern::VertexId{1}, pattern::VertexId{3});
    expect(split.has_value(), "polygon clipping setup split succeeds");
    if (!split) {
        return;
    }

    pattern::VectorFieldConfig config;
    config.seed_spacing = 0.25;
    config.step_size = 0.025;
    config.maximum_length = 3.0;
    const auto generator = pattern::VectorFieldCurveGenerator::create(config);
    expect(generator.has_value(), "polygon clipping generator configuration is valid");
    if (!generator) {
        return;
    }

    const auto generated = generator.value().generate(
        partition, pattern::RegionId{1});
    expect(generated.has_value(), "curves are generated on a triangular face");
    if (!generated) {
        return;
    }
    expect(!generated.value().curves.empty(),
        "triangular face contains curve seeds");
    for (const auto& curve : generated.value().curves) {
        for (const auto point : curve.points) {
            expect(point.v <= point.u * 0.5 + 1e-9,
                "curve samples do not leak across the shared diagonal");
        }
    }
}

}  // namespace

int main()
{
    test_generator_descriptor_and_parameters();
    test_curves_are_deterministic_and_clipped();
    test_missing_region_and_seed_limit();
    test_curves_clip_to_polygonal_region();

    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return EXIT_FAILURE;
    }

    std::cout << "All generator tests passed\n";
    return EXIT_SUCCESS;
}
