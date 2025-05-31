/*
 * Filename: result.hpp
 * Developer: Benjamin Cance
 * Date: 5/31/2025
 * 
 * Copyright 2025 Open Quant Desk, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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