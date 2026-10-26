#pragma once

#include "pattern/core/diagnostic.hpp"

#include <utility>
#include <variant>

namespace pattern {

template <typename T>
class [[nodiscard]] Result {
public:
    static Result success(T value)
    {
        return Result(std::move(value));
    }

    static Result failure(Diagnostic diagnostic)
    {
        return Result(std::move(diagnostic));
    }

    [[nodiscard]] bool has_value() const noexcept
    {
        return std::holds_alternative<T>(storage_);
    }

    explicit operator bool() const noexcept
    {
        return has_value();
    }

    [[nodiscard]] const T& value() const&
    {
        return std::get<T>(storage_);
    }

    [[nodiscard]] T& value() &
    {
        return std::get<T>(storage_);
    }

    [[nodiscard]] T&& value() &&
    {
        return std::get<T>(std::move(storage_));
    }

    [[nodiscard]] const Diagnostic& error() const&
    {
        return std::get<Diagnostic>(storage_);
    }

private:
    explicit Result(T value)
        : storage_(std::move(value))
    {
    }

    explicit Result(Diagnostic diagnostic)
        : storage_(std::move(diagnostic))
    {
    }

    std::variant<T, Diagnostic> storage_;
};

}  // namespace pattern
