#pragma once

#include <vector>
#include <array>
#include <random>
#include <iostream>
#include <cassert>
#include <ranges>
#include <Fixed.h>

template <typename P, typename V, typename FV, std::size_t rows_num, std::size_t cols_num>
struct FluidSimulation {
    static constexpr bool m_is_static = rows_num != max_size && cols_num != max_size;

    template <typename Type>
    using StaticStorage = std::array<std::array<Type, cols_num>, rows_num>;

    template <typename Type>
    using DynamicStorage = std::vector<std::vector<Type>>;

    template <typename Type>
    using StorageType = std::conditional_t<m_is_static, StaticStorage<Type>, DynamicStorage<Type>>;

    using Field = StorageType<char>;
    using PStorage = StorageType<P>;
    using IStorage = StorageType<int>;
    
    Field field{};
    P rho[256];
    std::size_t rows, cols;

    PStorage p{}, old_p{};
    IStorage dirs{}, last_use{};

    int UT = 0;
    V gravity;

    template <typename CurrentType>
    struct VelocityField {
        std::vector<std::vector<std::array<CurrentType, deltas.size()>>> v;

        VelocityField( void ) = default;
        VelocityField(std::size_t n, std::size_t m) {
            v.resize(n, std::vector<std::array<CurrentType, deltas.size()>>(m));
        }


        void f(std::size_t n, std::size_t m) {
            v.resize(n, std::vector<std::array<CurrentType, deltas.size()>>(m));
        }

        CurrentType& add(int x, int y, int dx, int dy, CurrentType dv) {
            return get(x, y, dx, dy) += dv;
        }

        CurrentType& get(int x, int y, int dx, int dy) {
            std::size_t index = std::ranges::find(deltas, std::pair(dx, dy)) - deltas.begin();
            assert(index < deltas.size());
            return v[x][y][index];
        }
    };

    struct ParticleParams {
        char type;
        P curr_p;
        std::array<V, deltas.size()> v;

        void swap_with(int x, int y, VelocityField<V>& vfield, PStorage& pfield, Field& fld) {
            std::swap(fld[x][y], type);
            std::swap(pfield[x][y], curr_p);
            std::swap(vfield.v[x][y], v);
        }
    };

    VelocityField<V> velocity{}, velocity_flow{};

    explicit constexpr FluidSimulation(const std::vector<std::vector<char>>& fld, float air_rho, int fluid_rho, float g) requires(m_is_static) : gravity(g) {
        assert(fld.size() == rows_num);
        assert(fld.front().size() == cols_num);
        // std::cout << "Static version with sizes: " << rows_num << " " << cols_num - 1 << std::endl;
        rows = fld.size();
        cols = fld.front().size() - 1;
        velocity.f(rows, cols);
        velocity_flow.f(rows, cols);
        for (std::size_t i = 0; i < field.size(); ++i) {
            const std::vector<char>& tmp_row = fld[i];
            assert(tmp_row.size() == field[i].size());
            for (std::size_t j = 0; j < field[i].size(); ++j) {
                field[i][j] = tmp_row[j];
            }
        }
        rho[' '] = air_rho;
        rho['.'] = fluid_rho;
    }

    explicit constexpr FluidSimulation(const Field& fld, float air_rho, int fluid_rho, float g) requires(!m_is_static): FluidSimulation(fld, air_rho, fluid_rho, g, fld.size(), fld.front().size()) {}

    explicit constexpr FluidSimulation(const Field& fld, float air_rho, int fluid_rho, float g, std::size_t rowCount, std::size_t colCount) requires(!m_is_static): gravity(g), field(fld),
                         p(rowCount, typename PStorage::value_type(colCount)),
                         old_p(rowCount, typename PStorage::value_type(colCount)),
                         dirs(rowCount, typename IStorage::value_type(colCount)),
                         last_use(rowCount, typename IStorage::value_type(colCount)) {
        // std::cout << "Dynamic version with sizes: " << fld.size() << " " << fld.front().size() - 1 << std::endl;
        rho[' '] = air_rho;
        rho['.'] = fluid_rho;
        rows = fld.size();
        cols = fld.front().size() - 1;
        velocity.f(rows, cols);
        velocity_flow.f(rows, cols);
    }

    void save_to_file( void ) {
        std::ofstream ofile;
        ofile.open("../state.txt");
        if (!ofile) {
            std::cerr << "Failed to open file! Error: " << std::strerror(errno) << std::endl;
            std::exit(-1);
        }

        ofile << rows << " " << cols << "\n";
        for (const auto& row : field) {
            for (char ch : row) {
                if (ch != '\0') {
                    ofile << ch;
                }
            }
            ofile << "\n";
        }

        ofile << rho[' '] << "\n" << rho['.'] << "\n" << gravity << "\n";
    }

    P random01( void ) {
        if constexpr (std::is_same<P, float>::value) {
            return static_cast<P>(float(rnd()) / rnd.max());
        } else if constexpr (std::is_same<P, double>::value) {
            return static_cast<P>(double(rnd()) / rnd.max());
        } else {
            return P::from_raw((rnd() & ((1 << P::getK()) - 1)));
        }
    }
    
    std::tuple<P, bool, std::pair<int, int>> propagate_flow(int x, int y, P lim) {
        last_use[x][y] = UT - 1;
        P ret = 0;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
                auto cap = velocity.get(x, y, dx, dy);
                auto flow = velocity_flow.get(x, y, dx, dy);
                if (flow == static_cast<FV>(cap)) {
                    continue;
                }
                // assert(v >= velocity_flow.get(x, y, dx, dy));
                V res = cap - static_cast<V>(flow);
                auto vp = std::min(lim, static_cast<P>(res));
                if (last_use[nx][ny] == UT - 1) {
                    velocity_flow.add(x, y, dx, dy, static_cast<FV>(vp));
                    last_use[x][y] = UT;
                    // cerr << x << " " << y << " -> " << nx << " " << ny << " " << vp << " / " << lim << "\n";
                    return {vp, 1, {nx, ny}};
                }
                auto [t, prop, end] = propagate_flow(nx, ny, vp);
                ret += t;
                if (prop) {
                    velocity_flow.add(x, y, dx, dy, static_cast<FV>(t));
                    last_use[x][y] = UT;
                    // cerr << x << " " << y << " -> " << nx << " " << ny << " " << t << " / " << lim << "\n";
                    return {t, prop && end != std::pair(x, y), end};
                }
            }
        }
        last_use[x][y] = UT;
        return {ret, 0, {0, 0}};
    }

    void propagate_stop(int x, int y, bool force = false) {
        if (!force) {
            bool stop = true;
            for (auto [dx, dy] : deltas) {
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > 0) {
                    stop = false;
                    break;
                }
            }
            if (!stop) {
                return;
            }
        }
        last_use[x][y] = UT;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT || velocity.get(x, y, dx, dy) > 0) {
                continue;
            }
            propagate_stop(nx, ny);
        }
    }

    P move_prob(int x, int y) {
        P sum = 0;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                continue;
            }
            auto v = velocity.get(x, y, dx, dy);
            if (v < static_cast<V>(0)) {
                continue;
            }
            sum += static_cast<P>(v);
        }
        return sum;
    }

    bool propagate_move(int x, int y, bool is_first) {
        last_use[x][y] = UT - is_first;
        bool ret = false;
        int nx = -1, ny = -1;
        do {
            std::array<P, deltas.size()> tres;
            P sum = 0;
            for (size_t i = 0; i < deltas.size(); ++i) {
                auto [dx, dy] = deltas[i];
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                    tres[i] = sum;
                    continue;
                }
                auto v = velocity.get(x, y, dx, dy);
                if (v < static_cast<V>(0)) {
                    tres[i] = sum;
                    continue;
                }
                sum += static_cast<P>(v);
                tres[i] = sum;
            }

            if (sum == static_cast<P>(0)) {
                break;
            }

            auto p = static_cast<P>(random01() * sum);
            size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

            auto [dx, dy] = deltas[d];
            nx = x + dx;
            ny = y + dy;
            assert(velocity.get(x, y, dx, dy) > static_cast<V>(0) && field[nx][ny] != '#' && last_use[nx][ny] < UT);

            ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
        } while (!ret);
        last_use[x][y] = UT;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < static_cast<V>(0)) {
                propagate_stop(nx, ny);
            }
        }
        if (ret) {
            if (!is_first) {
                ParticleParams pp{};
                pp.swap_with(x, y, velocity, p, field);
                pp.swap_with(nx, ny, velocity, p, field);
                pp.swap_with(x, y, velocity, p, field);
            }
        }
        return ret;
    }

    void run() {
        for (std::size_t x = 0; x < rows; ++x) {
            for (std::size_t y = 0; y < cols; ++y) {
                if (field[x][y] == '#')
                    continue;
                for (const auto& [dx, dy] : deltas) {
                    dirs[x][y] += (field[x + dx][y + dy] != '#');
                }
            }
        }

        for (std::size_t iteration = 0; iteration < T; ++iteration) {
            P totalDeltaPressure = 0;

            for (std::size_t x = 0; x < rows; ++x) {
                for (std::size_t y = 0; y < cols; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    if (field[x + 1][y] != '#')
                        velocity.add(x, y, 1, 0, gravity);
                }
            }

            for (std::size_t i = 0; i < p.size(); ++i) {
                for (std::size_t j = 0; j < p.front().size(); ++j) {
                    old_p[i][j] = p[i][j];
                }
            }

            for (std::size_t x = 0; x < rows; ++x) {
                for (std::size_t y = 0; y < cols; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    for (const auto& [dx, dy] : deltas) {
                        int newX = x + dx, newY = y + dy;
                        if (field[newX][newY] != '#' && old_p[newX][newY] < old_p[x][y]) {
                            auto deltaPressure = old_p[x][y] - old_p[newX][newY];
                            auto force = deltaPressure;
                            auto& velocityComponent = velocity.get(newX, newY, -dx, -dy);
                            if (static_cast<P>(velocityComponent) * rho[(int)field[newX][newY]] >= force) {
                                velocityComponent -= static_cast<V>(force / rho[(int)field[newX][newY]]);
                                continue;
                            }
                            force -= static_cast<P>(velocityComponent) * rho[(int)field[newX][newY]];
                            velocityComponent = 0;
                            velocity.add(x, y, dx, dy, static_cast<V>(force / rho[(int)field[x][y]]));
                            p[x][y] -= force / static_cast<P>(dirs[x][y]);
                            totalDeltaPressure -= force / static_cast<P>(dirs[x][y]);
                        }
                    }
                }
            }

            velocity_flow = {rows, cols};
            bool propagation = false;
            do {
                UT += 2;
                propagation = false;
                for (std::size_t x = 0; x < rows; ++x) {
                    for (std::size_t y = 0; y < cols; ++y) {
                        if (field[x][y] != '#' && last_use[x][y] != UT) {
                            auto [transferred, localPropagation, _] = propagate_flow(x, y, 1);
                            if (transferred > static_cast<P>(0)) {
                                propagation = true;
                            }
                        }
                    }
                }
            } while (propagation);

            for (std::size_t x = 0; x < rows; ++x) {
                for (std::size_t y = 0; y < cols; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    for (const auto& [dx, dy] : deltas) {
                        auto oldVelocity = velocity.get(x, y, dx, dy);
                        auto newVelocity = velocity_flow.get(x, y, dx, dy);
                        if (oldVelocity > static_cast<V>(0)) {
                            assert(newVelocity <= static_cast<FV>(oldVelocity));
                            velocity.get(x, y, dx, dy) = static_cast<V>(newVelocity);
                            auto force = static_cast<P>((oldVelocity - static_cast<V>(newVelocity))) * rho[(int)field[x][y]];
                            if (field[x][y] == '.')
                                force *= static_cast<P>(0.8);
                            if (field[x + dx][y + dy] == '#') {
                                p[x][y] += force / static_cast<P>(dirs[x][y]);
                                totalDeltaPressure += force / static_cast<P>(dirs[x][y]);
                            } else {
                                p[x + dx][y + dy] += force / static_cast<P>(dirs[x + dx][y + dy]);
                                totalDeltaPressure += force / static_cast<P>(dirs[x + dx][y + dy]);
                            }
                        }
                    }
                }
            }

            UT += 2;
            propagation = false;
            for (std::size_t x = 0; x < rows; ++x) {
                for (std::size_t y = 0; y < cols; ++y) {
                    if (field[x][y] != '#' && last_use[x][y] != UT) {
                        if (random01() < move_prob(x, y)) {
                            propagation = true;
                            propagate_move(x, y, true);
                        } else {
                            propagate_stop(x, y, true);
                        }
                    }
                }
            }

            if (propagation) {
                std::cout << "Tick " << iteration << ":\n";
                for (std::size_t x = 0; x < rows; ++x) {
                    for (std::size_t y = 0; y < cols; ++y) {
                        std::cout << (char)field[x][y];
                    }
                    std::cout << "\n";
                }
            }

            if (iteration % 10 == 0) {
                save_to_file();
            }
        }
    }
};