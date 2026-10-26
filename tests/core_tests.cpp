#include "pattern/surface/cylinder_patch.hpp"
#include "pattern/surface/domain2.hpp"
#include "pattern/surface/plane_patch.hpp"

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

bool close(double lhs, double rhs, double tolerance = 1e-12)
{
    return std::abs(lhs - rhs) <= tolerance;
}

void test_domain_validation()
{
    const auto valid = pattern::Domain2::rectangle(-1.0, 2.0, 3.0, 5.0);
    expect(valid.has_value(), "valid rectangular domain is accepted");
    if (!valid) {
        return;
    }

    expect(valid.value().contains({-1.0, 2.0}), "domain includes its boundary");
    expect(valid.value().contains({1.0, 4.0}), "domain contains interior points");
    expect(!valid.value().contains({3.1, 4.0}), "domain excludes exterior points");
    expect(!valid.value().contains({
        std::numeric_limits<double>::quiet_NaN(), 3.0
    }), "domain excludes non-finite query points");

    const auto zero_extent = pattern::Domain2::rectangle(0.0, 0.0, 0.0, 1.0);
    expect(!zero_extent, "zero-width domain is rejected");

    const auto infinite_bound = pattern::Domain2::rectangle(
        0.0, 0.0, std::numeric_limits<double>::infinity(), 1.0);
    expect(!infinite_bound, "non-finite domain bounds are rejected");

    const auto unrepresentable_area = pattern::Domain2::rectangle(
        0.0, 0.0, 1e200, 1e200);
    expect(!unrepresentable_area, "overflowing domain area is rejected");

    const auto underflowing_area = pattern::Domain2::rectangle(
        0.0, 0.0, 1e-200, 1e-200);
    expect(!underflowing_area, "underflowing domain area is rejected");
}

void test_plane_patch()
{
    auto domain = pattern::Domain2::rectangle(-2.0, -3.0, 2.0, 3.0);
    expect(domain.has_value(), "plane test domain is valid");
    if (!domain) {
        return;
    }

    const pattern::PlanePatch plane(domain.value());
    const auto sample = plane.evaluate({0.25, -0.5});
    expect(sample.has_value(), "plane evaluates points in its domain");
    if (sample) {
        expect(close(sample.value().position.x, 0.25), "plane x position");
        expect(close(sample.value().position.y, -0.5), "plane y position");
        expect(close(sample.value().position.z, 0.0), "plane z position");
        expect(close(sample.value().du.x, 1.0), "plane u tangent");
        expect(close(sample.value().dv.y, 1.0), "plane v tangent");
    }

    expect(!plane.evaluate({3.0, 0.0}), "plane rejects points outside its domain");
}

void test_cylinder_patch()
{
    auto domain = pattern::Domain2::rectangle(0.0, -1.0, 2.0, 1.0);
    expect(domain.has_value(), "cylinder test domain is valid");
    if (!domain) {
        return;
    }

    const auto invalid_radius = pattern::CylinderPatch::create(0.0, domain.value());
    expect(!invalid_radius, "zero cylinder radius is rejected");

    const auto cylinder = pattern::CylinderPatch::create(2.0, domain.value());
    expect(cylinder.has_value(), "positive cylinder radius is accepted");
    if (!cylinder) {
        return;
    }

    const auto sample = cylinder.value().evaluate({0.0, 0.5});
    expect(sample.has_value(), "cylinder evaluates points in its domain");
    if (sample) {
        expect(close(sample.value().position.x, 2.0), "cylinder x position");
        expect(close(sample.value().position.y, 0.0), "cylinder y position");
        expect(close(sample.value().position.z, 0.5), "cylinder z position");
        expect(close(sample.value().du.y, 2.0), "cylinder u tangent");
        expect(close(sample.value().dv.z, 1.0), "cylinder v tangent");
    }

    expect(!cylinder.value().evaluate({3.0, 0.0}),
        "cylinder rejects points outside its domain");
}

}  // namespace

int main()
{
    test_domain_validation();
    test_plane_patch();
    test_cylinder_patch();

    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return EXIT_FAILURE;
    }

    std::cout << "All core tests passed\n";
    return EXIT_SUCCESS;
}
