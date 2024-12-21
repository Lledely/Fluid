#ifndef FIXED_H
#define FIXED_H

#include <cstddef>
#include <compare>
#include <cstdint>

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

};

template <size_t N, size_t K>
using FastFixed = Fixed<N, K, true>;

} /* types */


#endif /* FIXED_H */