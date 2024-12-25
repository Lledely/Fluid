#pragma once

#include <iostream>
#include <limits>
#include <cstdint>

class NotDumbFixed {
public:
    int32_t v;

    constexpr NotDumbFixed(): v(0) {}
    constexpr NotDumbFixed( int val ): v(val << 16) {}
    constexpr NotDumbFixed( float val ): v(val * (1 << 16)) {}
    constexpr NotDumbFixed( double val ): v(val * (1 << 16)) {}

    auto operator<=>( const NotDumbFixed& other ) const = default;
    // bool operator==( const NotDumbFixed& other ) { return v == other.v; }
    // bool operator>( const NotDumbFixed& other ) { return v > other.v; }
    // bool operator<( const NotDumbFixed& other ) { return v < other.v; }
    // bool operator>=( const NotDumbFixed& other ) { return v >= other.v; }
    // bool operator<=( const NotDumbFixed& other ) { return v <= other.v; }
    // bool operator!=( const NotDumbFixed& other ) { return v != other.v; }

    static constexpr NotDumbFixed from_row( int32_t value ) {
        NotDumbFixed result;
        result.v = value;
        return result;
    }

    friend const NotDumbFixed operator+( const NotDumbFixed &a, const NotDumbFixed &b ) { return NotDumbFixed::from_row(a.v + b.v); }
    friend const NotDumbFixed operator-( const NotDumbFixed &a, const NotDumbFixed &b ) { return NotDumbFixed::from_row(a.v - b.v); }
    friend const NotDumbFixed operator*( const NotDumbFixed &a, const NotDumbFixed &b ) { return NotDumbFixed::from_row(((int64_t) a.v * b.v) >> 16); }
    friend const NotDumbFixed operator/( const NotDumbFixed &a, const NotDumbFixed &b ) { return NotDumbFixed::from_row(((int64_t) a.v << 16) / b.v); }

    NotDumbFixed& operator+=( const NotDumbFixed& b ) { return (*this) = (*this) + b; }
    NotDumbFixed& operator-=( const NotDumbFixed& b ) { return (*this) = (*this) - b; }
    NotDumbFixed& operator*=( const NotDumbFixed& b ) { return (*this) = (*this) * b; }
    NotDumbFixed& operator/=( const NotDumbFixed& b ) { return (*this) = (*this) / b; }

    NotDumbFixed operator-( void ) { return from_row(-(*this).v); }

    friend std::ostream& operator<<(std::ostream& out, NotDumbFixed x) {
        return out << x.v / (double) (1 << 16);
    }

    friend NotDumbFixed abs(NotDumbFixed x) {
        if (x.v < 0) {
            x.v = -x.v;
        }
        return x;
    }

};
