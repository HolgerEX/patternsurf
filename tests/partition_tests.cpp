#include "pattern/partition/partition.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string_view>
#include <vector>

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
    const auto domain = pattern::Domain2::rectangle(0.0, 0.0, 2.0, 1.0);
    if (!domain) {
        std::cerr << "Test setup failed: " << domain.error().message << '\n';
        std::exit(EXIT_FAILURE);
    }
    return pattern::Partition::create(domain.value());
}

void test_initial_partition()
{
    const pattern::Partition partition = make_partition();
    expect(partition.validate().empty(), "initial rectangular partition is valid");
    expect(partition.regions().size() == 1, "initial partition has one region");
    expect(partition.boundaries().size() == 4, "initial region has four boundaries");
    expect(partition.vertices().size() == 4, "initial region has four vertices");
}

void test_split_and_merge()
{
    pattern::Partition partition = make_partition();
    const pattern::RegionId original{1};
    const auto split = partition.split_region(
        original, pattern::VertexId{1}, pattern::VertexId{3});
    expect(split.has_value(), "diagonal split succeeds");
    if (!split) {
        return;
    }

    expect(partition.validate().empty(), "partition remains valid after split");
    expect(partition.regions().size() == 2, "split creates two regions");
    expect(partition.boundaries().size() == 5, "split creates one shared boundary");
    const auto neighbors = partition.adjacent_regions(original);
    expect(neighbors.has_value()
            && neighbors.value() == std::vector<pattern::RegionId>{split.value()},
        "adjacency query reports the neighboring region");
    const auto shared = std::find_if(
        partition.boundaries().begin(),
        partition.boundaries().end(),
        [](const auto& item) {
            return item.second.left_region && item.second.right_region;
        });
    expect(shared != partition.boundaries().end(),
        "both split regions reference one shared boundary entity");
    if (shared != partition.boundaries().end()) {
        expect(shared->second.left_region == split.value()
                && shared->second.right_region == original,
            "shared boundary records both incident regions");
    }

    const auto merged = partition.merge_regions(original, split.value());
    expect(merged.has_value(), "neighbor merge succeeds");
    if (merged) {
        expect(merged.value() == original, "merge retains the selected stable region ID");
    }
    expect(partition.validate().empty(), "partition remains valid after merge");
    expect(partition.regions().size() == 1, "merge restores one region");
    expect(partition.boundaries().size() == 4, "merge removes the internal boundary");
}

void test_boundary_insertion_and_split()
{
    pattern::Partition partition = make_partition();
    const auto inserted = partition.insert_boundary_vertex(
        pattern::BoundaryId{1}, {1.0, 0.0});
    expect(inserted.has_value(), "boundary vertex insertion succeeds");
    if (!inserted) {
        return;
    }
    expect(partition.validate().empty(), "inserted shared edge segments remain valid");
    expect(partition.vertices().size() == 5, "insertion adds one vertex");

    const auto split = partition.split_region(
        pattern::RegionId{1}, inserted.value(), pattern::VertexId{3});
    expect(split.has_value(), "split using inserted boundary vertex succeeds");
    expect(partition.validate().empty(),
        "partition remains valid after boundary insertion and split");
    expect(partition.regions().size() == 2,
        "boundary insertion and split produce two regions");
}

void test_shared_vertex_motion_is_transactional()
{
    pattern::Partition partition = make_partition();
    const auto split = partition.split_region(
        pattern::RegionId{1}, pattern::VertexId{1}, pattern::VertexId{3});
    expect(split.has_value(), "setup split creates a shared boundary");
    if (!split) {
        return;
    }

    const auto shared = std::find_if(
        partition.boundaries().begin(),
        partition.boundaries().end(),
        [](const auto& item) {
            return item.second.left_region && item.second.right_region;
        });
    if (shared == partition.boundaries().end()) {
        expect(false, "setup has an internal boundary");
        return;
    }

    const pattern::Boundary boundary = shared->second;
    const pattern::ParameterPoint from = partition.vertices().at(boundary.from);
    const pattern::ParameterPoint to = partition.vertices().at(boundary.to);
    const auto inserted = partition.insert_boundary_vertex(
        boundary.id,
        {(from.u + to.u) * 0.5, (from.v + to.v) * 0.5});
    expect(inserted.has_value(), "internal boundary vertex insertion succeeds");
    if (!inserted) {
        return;
    }

    const auto moved = partition.move_vertex(inserted.value(), {1.1, 0.5});
    expect(moved.has_value(), "shared interior vertex moves successfully");
    expect(partition.validate().empty(),
        "shared vertex motion preserves both incident region cycles");
    if (moved) {
        expect(partition.vertices().at(moved.value()).u == 1.1,
            "moved position is published");
    }

    const auto rejected = partition.move_vertex(inserted.value(), {3.0, 0.5});
    expect(!rejected, "out-of-domain vertex motion is rejected");
    expect(partition.vertices().at(inserted.value()).u == 1.1,
        "rejected vertex motion rolls back");
    expect(partition.validate().empty(),
        "partition remains valid after rejected vertex motion");

    const auto perimeter_move = partition.move_vertex(
        pattern::VertexId{1}, {0.1, 0.0});
    expect(!perimeter_move, "domain-perimeter vertices remain fixed");
}

void test_shared_boundary_insertion_and_merge()
{
    pattern::Partition partition = make_partition();
    const auto split = partition.split_region(
        pattern::RegionId{1}, pattern::VertexId{1}, pattern::VertexId{3});
    expect(split.has_value(), "setup split creates adjacent regions");
    if (!split) {
        return;
    }

    const auto shared = std::find_if(
        partition.boundaries().begin(),
        partition.boundaries().end(),
        [](const auto& item) {
            return item.second.left_region && item.second.right_region;
        });
    expect(shared != partition.boundaries().end(), "split has shared interface");
    if (shared == partition.boundaries().end()) {
        return;
    }

    const pattern::Boundary boundary = shared->second;
    const pattern::ParameterPoint from = partition.vertices().at(boundary.from);
    const pattern::ParameterPoint to = partition.vertices().at(boundary.to);
    const auto inserted = partition.insert_boundary_vertex(
        boundary.id,
        {(from.u + to.u) * 0.5, (from.v + to.v) * 0.5});
    expect(inserted.has_value(), "inserting a shared-interface vertex succeeds");
    expect(partition.validate().empty(),
        "both incident regions stay synchronized after shared-interface insertion");

    std::size_t shared_segments = 0;
    for (const auto& [id, edge] : partition.boundaries()) {
        (void)id;
        if (edge.left_region && edge.right_region) {
            ++shared_segments;
        }
    }
    expect(shared_segments == 2, "shared interface is split into two shared entities");

    const auto merged = partition.merge_regions(pattern::RegionId{1}, split.value());
    expect(merged.has_value(),
        "merge handles a shared interface made of multiple contiguous segments");
    expect(partition.validate().empty(),
        "partition remains valid after merging a segmented interface");
    expect(partition.regions().size() == 1,
        "segmented-interface merge produces one region");
}

void test_split_of_region_with_existing_neighbor()
{
    pattern::Partition partition = make_partition();
    const auto first_split = partition.split_region(
        pattern::RegionId{1}, pattern::VertexId{1}, pattern::VertexId{3});
    expect(first_split.has_value(), "setup split creates neighboring faces");
    if (!first_split) {
        return;
    }

    const auto top_boundary = std::find_if(
        partition.boundaries().begin(),
        partition.boundaries().end(),
        [](const auto& item) {
            return item.second.from == pattern::VertexId{3}
                && item.second.to == pattern::VertexId{4};
        });
    expect(top_boundary != partition.boundaries().end(),
        "setup has the expected top boundary");
    if (top_boundary == partition.boundaries().end()) {
        return;
    }

    const auto inserted = partition.insert_boundary_vertex(
        top_boundary->first, {1.0, 1.0});
    expect(inserted.has_value(), "insertion on neighboring region boundary succeeds");
    if (!inserted) {
        return;
    }

    const auto second_split = partition.split_region(
        first_split.value(), inserted.value(), pattern::VertexId{1});
    expect(second_split.has_value(), "adjacent region can be split independently");
    expect(partition.validate().empty(),
        "split updates shared-boundary incidence without breaking neighbors");
    expect(partition.regions().size() == 3,
        "splitting a neighboring face adds one region");
}

void test_small_scale_partition()
{
    const auto domain = pattern::Domain2::rectangle(0.0, 0.0, 1e-9, 1e-9);
    expect(domain.has_value(), "small but representable domain is accepted");
    if (!domain) {
        return;
    }
    const pattern::Partition partition = pattern::Partition::create(domain.value());
    expect(partition.validate().empty(),
        "scale-aware tolerances accept a small valid partition");
}

void test_deletion_and_atomic_failure()
{
    pattern::Partition partition = make_partition();
    const std::size_t original_boundary_count = partition.boundaries().size();
    const auto rejected = partition.split_region(
        pattern::RegionId{1}, pattern::VertexId{1}, pattern::VertexId{2});
    expect(!rejected, "split along an existing edge is rejected");
    expect(partition.regions().size() == 1
            && partition.boundaries().size() == original_boundary_count,
        "rejected split leaves the original partition unchanged");
    expect(partition.validate().empty(), "partition remains valid after rollback");

    const auto deleted = partition.delete_region(pattern::RegionId{1});
    expect(deleted.has_value(), "deleting a region replaces it with background");
    if (!deleted) {
        return;
    }
    expect(deleted.value() == pattern::RegionId{2},
        "rejected edits do not consume stable region IDs");

    const pattern::Region* background = partition.find_region(deleted.value());
    expect(background != nullptr && background->background,
        "replacement region is explicitly marked as background");
    expect(partition.validate().empty(),
        "background replacement preserves full-domain coverage");

    const auto second_delete = partition.delete_region(deleted.value());
    expect(!second_delete, "background region cannot be deleted under full coverage");
    expect(partition.validate().empty(),
        "failed background deletion leaves the partition valid");
}

void test_invalid_boundary_vertex_rolls_back()
{
    pattern::Partition partition = make_partition();
    const auto rejected = partition.insert_boundary_vertex(
        pattern::BoundaryId{1}, {1.0, 0.25});
    expect(!rejected, "vertex off the boundary segment is rejected");
    expect(partition.vertices().size() == 4
            && partition.boundaries().size() == 4,
        "rejected vertex insertion leaves topology unchanged");
    expect(partition.validate().empty(),
        "partition remains valid after rejected vertex insertion");

    const auto inserted = partition.insert_boundary_vertex(
        pattern::BoundaryId{1}, {1.0, 0.0});
    expect(inserted.has_value() && inserted.value() == pattern::VertexId{5},
        "rejected insertion does not consume stable vertex IDs");
}

void test_merge_preserves_other_adjacencies()
{
    pattern::Partition partition = make_partition();
    const auto first_split = partition.split_region(
        pattern::RegionId{1}, pattern::VertexId{1}, pattern::VertexId{3});
    expect(first_split.has_value(), "setup split succeeds");
    if (!first_split) {
        return;
    }
    const auto second_split = partition.split_region(
        first_split.value(), pattern::VertexId{3}, pattern::VertexId{1});
    expect(!second_split, "repeated split along the existing shared edge is rejected");
    expect(partition.validate().empty(),
        "rejected repeated split leaves neighboring regions valid");

    const auto inserted = partition.insert_boundary_vertex(
        pattern::BoundaryId{3}, {1.0, 1.0});
    expect(inserted.has_value(), "boundary vertex inserted for second split");
    if (!inserted) {
        return;
    }
    const auto second_face = partition.split_region(
        first_split.value(), inserted.value(), pattern::VertexId{1});
    expect(second_face.has_value(), "neighboring region split succeeds");
    if (!second_face) {
        return;
    }

    const auto merged = partition.merge_regions(
        first_split.value(), second_face.value());
    expect(merged.has_value(), "merge of child faces succeeds");
    expect(partition.validate().empty(),
        "merge preserves the other region's shared-boundary adjacency");
    expect(partition.regions().size() == 2,
        "merging split children leaves the untouched neighbor active");
}

}  // namespace

int main()
{
    test_initial_partition();
    test_split_and_merge();
    test_boundary_insertion_and_split();
    test_shared_vertex_motion_is_transactional();
    test_shared_boundary_insertion_and_merge();
    test_split_of_region_with_existing_neighbor();
    test_merge_preserves_other_adjacencies();
    test_small_scale_partition();
    test_deletion_and_atomic_failure();
    test_invalid_boundary_vertex_rolls_back();

    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return EXIT_FAILURE;
    }

    std::cout << "All partition tests passed\n";
    return EXIT_SUCCESS;
}
