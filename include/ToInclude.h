#pragma once

#include <cstddef>
#include <compare>
#include <cstdint>
#include <random>
#include <iostream>
#include <array>

#ifdef FLOAT
#error "FLOAT is already defined"
#endif
#ifdef DOUBLE
#error "DOUBLE is already defined"
#endif
#ifdef FIXED
#error "FIXED is already defined"
#endif
#ifdef FAST_FIXED
#error "FAST_FIXED is already defined"
#endif

#define FLOAT float
#define DOUBLE double
#define FAST_FIXED(N, K) types::FastFixed<N, K>
#define FIXED(N, K) types::Fixed<N, K>

constexpr std::string_view kFloatTypeName   = "FLOAT";
constexpr std::string_view kDoubleTypeName  = "DOUBLE";

constexpr size_t T = 1'000'000;
constexpr std::array<std::pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

constexpr std::size_t max_size = std::numeric_limits<std::size_t>::max();
constexpr size_t N = 36, M = 84;
constexpr bool withprint = false;
inline std::mt19937 rnd(1337);