#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <thread>
#include <omp.h>
#include <ToInclude.h>
#include <NotDumbFixed.h>

using namespace std;

static constexpr NotDumbFixed INF = NotDumbFixed::from_row(std::numeric_limits<int32_t>::max());
static constexpr NotDumbFixed EPS = NotDumbFixed::from_row(4);

struct NotDumbFluidSimulation {
    char field[N][M + 1] = {
            "####################################################################################",
            "#                                                                                  #",
            "#                                                                                  #",
            "#                                                                                  #",
            "#                                                                                  #",
            "#                                                                                  #",
            "#                                       .........                                  #",
            "#..............#            #           .........                                  #",
            "#..............#            #           .........                                  #",
            "#..............#            #           .........                                  #",
            "#..............#            #                                                      #",
            "#..............#            #                                                      #",
            "#..............#            #                                                      #",
            "#..............#            #                                                      #",
            "#..............#............#                                                      #",
            "#..............#............#                                                      #",
            "#..............#............#                                                      #",
            "#..............#............#                                                      #",
            "#..............#............#                                                      #",
            "#..............#............#                                                      #",
            "#..............#............#                                                      #",
            "#..............#............#                                                      #",
            "#..............#............################                     #                 #",
            "#...........................#....................................#                 #",
            "#...........................#....................................#                 #",
            "#...........................#....................................#                 #",
            "##################################################################                 #",
            "#                                                                                  #",
            "#                                                                                  #",
            "#                                                                                  #",
            "#                                                                                  #",
            "#                                                                                  #",
            "#                                                                                  #",
            "#                                                                                  #",
            "#                                                                                  #",
            "####################################################################################",
    };

    NotDumbFixed rho[256];

    NotDumbFixed p[N][M]{}, old_p[N][M];

    struct VectorField {
        array<NotDumbFixed, deltas.size()> v[N][M];
        NotDumbFixed &add(int x, int y, int dx, int dy, NotDumbFixed dv) {
            return get(x, y, dx, dy) += dv;
        }

        NotDumbFixed &get(int x, int y, int dx, int dy) {
            if (deltas[0] == pair(dx, dy)) {
                return v[x][y][0];
            } else if (deltas[1] == pair(dx, dy)) {
                return v[x][y][1];
            } else if (deltas[2] == pair(dx, dy)) {
                return v[x][y][2];
            } else {
                return v[x][y][3];
            }
        }
    };

    VectorField velocity{}, velocity_flow{};
    int last_use[N][M]{};
    int UT = 0;
    inline tuple<NotDumbFixed, bool, pair<int, int>> propagate_flow(int x, int y, NotDumbFixed lim) {
        last_use[x][y] = UT - 1;
        NotDumbFixed ret = 0;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
                auto cap = velocity.get(x, y, dx, dy);
                auto flow = velocity_flow.get(x, y, dx, dy);
                if (flow == cap) {
                    continue;
                }
                auto vp = min(lim, cap - flow);
                if (last_use[nx][ny] == UT - 1) {
                    velocity_flow.add(x, y, dx, dy, vp);
                    last_use[x][y] = UT;
                    return {vp, 1, {nx, ny}};
                }
                auto [t, prop, end] = propagate_flow(nx, ny, vp);
                ret += t;
                if (prop) {
                    velocity_flow.add(x, y, dx, dy, t);
                    last_use[x][y] = UT;
                    return {t, prop && end != pair(x, y), end};
                }
            }
        }
        last_use[x][y] = UT;
        return {ret, 0, {0, 0}};
    }

    NotDumbFixed random01() {
        return NotDumbFixed::from_row((rnd() & ((1 << 16) - 1)));
    }

    inline void propagate_stop(int x, int y, bool force = false) {
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

    inline NotDumbFixed move_prob(int x, int y) {
        NotDumbFixed sum = 0;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                continue;
            }
            auto v = velocity.get(x, y, dx, dy);
            if (v < 0) {
                continue;
            }
            sum += v;
        }
        return sum;
    }

    struct ParticleParams {
        char type;
        NotDumbFixed cur_p;
        array<NotDumbFixed, deltas.size()> v;

        void swap_with(int x, int y, char (&field1)[N][M + 1], NotDumbFixed (&p1)[N][M], VectorField &velocity1) {
            swap(field1[x][y], type);
            swap(p1[x][y], cur_p);
            swap(velocity1.v[x][y], v);
        }
    };

    inline bool propagate_move(int x, int y, bool is_first) {
        last_use[x][y] = UT - is_first;
        bool ret = false;
        int nx = -1, ny = -1;
        do {
            std::array<NotDumbFixed, deltas.size()> tres;
            NotDumbFixed sum = 0;
            for (size_t i = 0; i < deltas.size(); ++i) {
                auto [dx, dy] = deltas[i];
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                    tres[i] = sum;
                    continue;
                }
                auto v = velocity.get(x, y, dx, dy);
                if (v < 0) {
                    tres[i] = sum;
                    continue;
                }
                sum += v;
                tres[i] = sum;
            }

            if (sum == 0) {
                break;
            }

            NotDumbFixed p = random01() * sum;
            size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

            auto [dx, dy] = deltas[d];
            nx = x + dx;
            ny = y + dy;
            assert(velocity.get(x, y, dx, dy) > 0 && field[nx][ny] != '#' && last_use[nx][ny] < UT);

            ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
        } while (!ret);
        last_use[x][y] = UT;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < 0) {
                propagate_stop(nx, ny);
            }
        }
        if (ret) {
            if (!is_first) {
                ParticleParams pp{};
                pp.swap_with(x, y, field, p, velocity);
                pp.swap_with(nx, ny, field, p, velocity);
                pp.swap_with(x, y, field, p, velocity);
            }
        }
        return ret;
    }

    int dirs[N][M]{};

    inline void func(int left, int right) {
        for (size_t x = left; x < right; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    dirs[x][y] += (field[x + dx][y + dy] != '#');
                }
            }
        }
    }

    inline int run(int ticks, int numThreads) {
        omp_set_num_threads(numThreads);
        rho[' '] = 0.01;
        rho['.'] = 1000;
        NotDumbFixed g = 0.1;

        vector<thread> threads;
        size_t chunkSize = N / numThreads;

        for (int i = 0; i < numThreads; ++i) {
            int left = i * chunkSize;
            int right = (i == numThreads - 1) ? N : left + chunkSize;

            threads.emplace_back([this, left, right]() {
                func(left, right);
            });
        }

        for (auto& t : threads) {
            t.join();
        }



        for (size_t i = 0; i < ticks; ++i) {
            NotDumbFixed total_delta_p = 0;
            // Apply external forces
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    if (field[x + 1][y] != '#')
                        velocity.add(x, y, 1, 0, g);
                }
            }

            // Apply forces from p
            memcpy(old_p, p, sizeof(p));
            #pragma omp parallel for
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    for (auto [dx, dy] : deltas) {
                        int nx = x + dx, ny = y + dy;
                        if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
                            auto delta_p = old_p[x][y] - old_p[nx][ny];
                            auto force = delta_p;
                            auto &contr = velocity.get(nx, ny, -dx, -dy);
                            if (contr * rho[(int) field[nx][ny]] >= force) {
                                contr -= force / rho[(int) field[nx][ny]];
                                continue;
                            }
                            force -= contr * rho[(int) field[nx][ny]];
                            contr = 0;
                            velocity.add(x, y, dx, dy, force / rho[(int) field[x][y]]);
                            p[x][y] -= force / dirs[x][y];
                            total_delta_p -= force / dirs[x][y];
                        }
                    }
                }
            }

            // Make flow from velocities
            velocity_flow = {};
            bool prop = false;
            do {
                UT += 2;
                prop = 0;
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] != '#' && last_use[x][y] != UT) {
                            auto [t, local_prop, _] = propagate_flow(x, y, 1);
                            if (t > 0) {
                                prop = 1;
                            }
                        }
                    }
                }
            } while (prop);

            // Recalculate p with kinetic energy
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    for (auto [dx, dy] : deltas) {
                        auto old_v = velocity.get(x, y, dx, dy);
                        auto new_v = velocity_flow.get(x, y, dx, dy);
                        if (old_v > 0) {
                            assert(new_v <= old_v);
                            velocity.get(x, y, dx, dy) = new_v;
                            auto force = (old_v - new_v) * rho[(int) field[x][y]];
                            if (field[x][y] == '.')
                                force *= 0.8;
                            if (field[x + dx][y + dy] == '#') {
                                p[x][y] += force / dirs[x][y];
                                total_delta_p += force / dirs[x][y];
                            } else {
                                p[x + dx][y + dy] += force / dirs[x + dx][y + dy];
                                total_delta_p += force / dirs[x + dx][y + dy];
                            }
                        }
                    }
                }
            }

            UT += 2;
            prop = false;
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] != '#' && last_use[x][y] != UT) {
                        if (random01() < move_prob(x, y)) {
                            prop = true;
                            propagate_move(x, y, true);
                        } else {
                            propagate_stop(x, y, true);
                        }
                    }
                }
            }

            if (prop) {
                if (withprint) {
                    cout << "Tick " << i << ":\n";

                    for (size_t x = 0; x < N; ++x) {
                        cout << field[x] << "\n";
                    }
                }
            }
        }
    }
};
