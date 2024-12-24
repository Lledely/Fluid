#pragma once

#include <vector>
#include <array>
#include <random>
#include <iostream>
#include <cassert>
#include <ranges>
#include <Fixed.h>

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

constexpr size_t T = 1'000'000;
constexpr std::array<std::pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

template <class MatrixElementType>
    using ArrayType = std::conditional_t<IsStatic,
            StaticMatrixStorage<MatrixElementType>,
            DynamicMatrixStorage<MatrixElementType>>;

template <typename current_type>
struct VectorField {
    std::vector<std::vector<std::array<current_type, deltas.size()>>> vec;

    VectorField( void ) = default;
    VectorField(size_t n, size_t m) {
        vec.assign(n, std::vector<std::array<current_type, deltas.size()>>(m), 0);
    }

    current_type &add(int x, int y, int dx, int dy, current_type dv) {
        return get(x, y, dx, dy) += dv;
    }

    current_type &get(int x, int y, int dx, int dy) {
        size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
        assert(i < deltas.size());
        return vec[x][y][i];
    }
};

template<typename P, typename V, typename VF>
struct ParticleParams {
    char type;
    P current_p;
    std::array<V, deltas.size()> vel;

    void swap_partincles( size_t x, size_t y, VectorField<V> velocity, ArrayType<P>, &part, ArrayType<char> &field ) {
        swap(field[x][y], type);
        swap(p[x][y], current_p);
        swap(velocity.vec[x][y], vel);
    }
};

template <typename P, typename V, typename VF, size_t rows, size_t columns>
struct FluidSimulation {
    
};

