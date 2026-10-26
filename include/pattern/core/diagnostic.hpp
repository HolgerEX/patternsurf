#pragma once

#include <string>

namespace pattern {

enum class Severity {
    error,
    warning
};

enum class ErrorCode {
    invalid_argument,
    outside_domain,
    not_found,
    invalid_topology,
    unsupported_operation
};

struct Diagnostic {
    ErrorCode code;
    Severity severity;
    std::string message;
};

}  // namespace pattern
