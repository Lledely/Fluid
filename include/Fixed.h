#pragma once

#include <cstddef>
#include <compare>
#include <cstdint>
#include <random>
#include <iostream>
#include <array>
#include <ToInclude.h>

namespace types {
    
template <size_t N, size_t K, bool is_fast=false>
class Fixed {
public:

    static constexpr size_t m_N = N;
    static constexpr size_t m_K = K;
    static constexpr bool m_is_fast = is_fast;
    int64_t v;

    constexpr Fixed( void ): v(0) {}
    constexpr Fixed( int32_t value ): v(static_cast<int64_t>(value) << 16) {}
    constexpr Fixed( float value ): v(static_cast<int64_t>(value * (1 << 16))) {}
    constexpr Fixed( double value ): v(static_cast<int64_t>(value * (1 << 16))) {}

    auto operator<=>(const Fixed& other) const = default;

    static constexpr Fixed from_row( int32_t x ) {
        Fixed ret;
        ret.v = x;
        return ret;
    }

    friend const Fixed operator+( const Fixed &a, const Fixed &b) { return from_raw(a.v + b.v); }
    friend const Fixed operator-( const Fixed &a, const Fixed &b) { return from_raw(a.v - b.v); }
    friend const Fixed operator*( const Fixed &a, const Fixed &b) { return from_raw(((int64_t) a.v * b.v) >> 16); }
    friend const Fixed operator/( const Fixed &a, const Fixed &b) { return from_raw(((int64_t) a.v << 16) / b.v); }

    constexpr Fixed& operator+=( const Fixed &b ) { return (*this) = (*this) + b; }
    constexpr Fixed& operator-=( const Fixed &b ) { return (*this) = (*this) - b; }
    constexpr Fixed& operator*=( const Fixed &b ) { return (*this) = (*this) * b; }
    constexpr Fixed& operator/=( const Fixed &b ) { return (*this) = (*this) / b; }

    Fixed operator-( void ) { return from_raw(-(*this).v); }

    friend std::ostream &operator<<(std::ostream &out, Fixed x) {
        return out << x.v / (double) (1 << 16);
    }

    friend Fixed abs(Fixed x) {
        if (x.v < 0) {
            x.v = -x.v;
        }
        return x;
    }
};

template <size_t N, size_t K>
using FastFixed = Fixed<N, K, true>;

} /* types */
