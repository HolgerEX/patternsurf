cmake_minimum_required(VERSION 3.25)

project(PatternSurface VERSION 0.1.0 LANGUAGES CXX)

include(CTest)
include(CheckCXXSourceCompiles)

option(PATTERNSURF_ENABLE_SANITIZERS "Enable address and undefined-behavior sanitizers" OFF)

if(PATTERNSURF_ENABLE_SANITIZERS)
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        message(FATAL_ERROR
            "PATTERNSURF_ENABLE_SANITIZERS requires GCC or Clang")
    endif()

    set(CMAKE_REQUIRED_FLAGS
        "${CMAKE_REQUIRED_FLAGS} -fsanitize=address,undefined -fno-omit-frame-pointer")
    set(CMAKE_REQUIRED_LINK_OPTIONS "-fsanitize=address,undefined")
    check_cxx_source_compiles("int main() { return 0; }"
        PATTERNSURF_SANITIZERS_SUPPORTED)
    unset(CMAKE_REQUIRED_FLAGS)
    unset(CMAKE_REQUIRED_LINK_OPTIONS)

    if(NOT PATTERNSURF_SANITIZERS_SUPPORTED)
        message(FATAL_ERROR
            "The selected compiler cannot link AddressSanitizer and "
            "UndefinedBehaviorSanitizer. Install their runtime libraries or "
            "use a supported GCC/Clang toolchain.")
    endif()
endif()

add_library(pattern_core
    src/core/domain2.cpp
    src/surface/plane_patch.cpp
    src/surface/cylinder_patch.cpp
)
add_library(PatternSurface::Core ALIAS pattern_core)

add_library(pattern_partition src/partition/partition.cpp)
add_library(PatternSurface::Partition ALIAS pattern_partition)
target_link_libraries(pattern_partition PUBLIC PatternSurface::Core)
target_include_directories(pattern_partition
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_compile_features(pattern_partition PUBLIC cxx_std_20)

target_include_directories(pattern_core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_compile_features(pattern_core PUBLIC cxx_std_20)

add_executable(pattern_cli apps/cli/main.cpp)
target_link_libraries(pattern_cli PRIVATE PatternSurface::Core)
target_compile_definitions(pattern_cli PRIVATE
    PATTERNSURF_VERSION="${PROJECT_VERSION}"
)

function(pattern_enable_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive-)
    else()
        target_compile_options(${target} PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wshadow
        )
    endif()
endfunction()

function(pattern_enable_sanitizers target)
    if(PATTERNSURF_ENABLE_SANITIZERS)
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined
            -fno-omit-frame-pointer
        )
        target_link_options(${target} PRIVATE
            -fsanitize=address,undefined
            -fno-omit-frame-pointer
        )
    endif()
endfunction()

pattern_enable_warnings(pattern_core)
pattern_enable_warnings(pattern_partition)
pattern_enable_warnings(pattern_cli)
pattern_enable_sanitizers(pattern_core)
pattern_enable_sanitizers(pattern_partition)
pattern_enable_sanitizers(pattern_cli)

if(BUILD_TESTING)
    add_executable(pattern_tests tests/core_tests.cpp)
    target_link_libraries(pattern_tests PRIVATE PatternSurface::Core)
    target_compile_features(pattern_tests PRIVATE cxx_std_20)
    pattern_enable_warnings(pattern_tests)
    pattern_enable_sanitizers(pattern_tests)
    add_test(NAME pattern_core_tests COMMAND pattern_tests)

    add_executable(pattern_partition_tests tests/partition_tests.cpp)
    target_link_libraries(pattern_partition_tests PRIVATE PatternSurface::Partition)
    target_compile_features(pattern_partition_tests PRIVATE cxx_std_20)
    pattern_enable_warnings(pattern_partition_tests)
    pattern_enable_sanitizers(pattern_partition_tests)
    add_test(NAME pattern_partition_tests COMMAND pattern_partition_tests)
endif()
