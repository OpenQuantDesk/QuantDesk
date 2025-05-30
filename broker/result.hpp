#pragma once

#include <variant>
#include <string>

namespace broker {

template<typename T>
using Result = std::variant<T, std::string>;

template<typename T>
bool isSuccess(const Result<T>& result) {
    return std::holds_alternative<T>(result);
}

template<typename T>
const T& getValue(const Result<T>& result) {
    return std::get<T>(result);
}

template<typename T>
const std::string& getError(const Result<T>& result) {
    return std::get<std::string>(result);
}

template<typename T>
Result<T> makeSuccess(T&& value) {
    return Result<T>(std::forward<T>(value));
}

template<typename T>
Result<T> makeError(const std::string& error) {
    return Result<T>(error);
}

template<typename T>
Result<T> makeError(std::string&& error) {
    return Result<T>(std::move(error));
}

}