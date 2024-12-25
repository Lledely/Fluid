#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <string_view>
#include <cstddef>
#include <type_traits>
#include <cctype>
#include <cerrno>
#include <limits>
#include <stdexcept>
#include <FluidSimulation.h>

using namespace std;

template <class... Types>
struct TypeList {};

template <std::size_t I, class Type, class... Types>
struct TypeIterator {
    static_assert(I < 1 + sizeof...(Types));
    using type = typename TypeIterator<I - 1, Types...>::type;
};

template <class Type, class... Types>
struct TypeIterator<0, Type, Types...> {
    using type = Type;
};

template <std::size_t I, class... Types>
using TypeAt = typename TypeIterator<I, Types...>::type;

struct Dimensions {
    size_t rows = 0;
    size_t columns = 0;
};

#define DIM(N, M) Dimensions { .rows = N, .columns = M }

struct SimulationParams {
    float rho_air;
    int rho_fluid;
    float g;
};

template <class P, class V, class FV, size_t rows = max_size, size_t cols = max_size>
void run_fluid_simulation(const vector<vector<char>>& field, const SimulationParams params) {
    FluidSimulation<P, V, FV, rows, cols>(field, params.rho_air, params.rho_fluid, params.g).run();
}

template <class TypeP, class V, class FV, Dimensions curSize, Dimensions... Sizes>
void select_static_sizes_impl(const vector<vector<char>>& field, const SimulationParams params) {
    if (curSize.rows == field.size() && curSize.columns == field.front().size() - 1) {
        run_fluid_simulation<TypeP, V, FV, curSize.rows, curSize.columns + 1>(field, params);
    } else if constexpr (sizeof...(Sizes) == 0) {
        run_fluid_simulation<TypeP, V, FV>(field, params);
    } else {
        select_static_sizes_impl<TypeP, V, FV, Sizes...>(field, params);
    }
}

template <class TypeP, class V, class FV, Dimensions... Sizes>
void select_static_sizes(const vector<vector<char>>& field, const SimulationParams params) {
    if constexpr (sizeof...(Sizes) > 0) {
        select_static_sizes_impl<TypeP, V, FV, Sizes...>(field, params);
    } else {
        run_fluid_simulation<TypeP, V, FV>(field, params);
    }
}

#define DIMENSIONS DIM(10, 10), DIM(36, 84)
#define TYPE_LIST FIXED(32, 16), FLOAT, DOUBLE, FAST_FIXED(32, 16)

template <class DefinedTypes, class SelectedTypes, Dimensions... Sizes>
class TypeSelector;

template <class... DefinedTypes, class... SelectedTypes, Dimensions... Sizes>
class TypeSelector<TypeList<DefinedTypes...>, TypeList<SelectedTypes...>, Sizes...> {
    static constexpr bool HasFloat = disjunction_v<is_same<DefinedTypes, FLOAT>...>;
    static constexpr bool HasDouble = disjunction_v<is_same<DefinedTypes, DOUBLE>...>;

public:
    static bool parse_fixed_string(string_view input, size_t& n, size_t& k, bool isFast) {
        string_view prefix_fixed = "FIXED(";
        string_view prefix_fast_fixed = "FAST_FIXED(";

        if (input.substr(0, prefix_fixed.size()) == prefix_fixed) {
            if (isFast) return false;
            input.remove_prefix(prefix_fixed.size());
        } else if (input.substr(0, prefix_fast_fixed.size()) == prefix_fast_fixed) {
            if (!isFast) return false;
            input.remove_prefix(prefix_fast_fixed.size());
        } else {
            return false;
        }

        size_t pos = 0;
        n = atoi(input.data());
        while (pos < input.size() && isdigit(input[pos])) ++pos;
        if (pos >= input.size() || input[pos] != ',') return false;
        ++pos;
        while (pos < input.size() && isspace(input[pos])) ++pos;
        k = atoi(input.data() + pos);
        return true;
    }

    template <class... Args>
    static void process_types(const vector<vector<char>>& field, SimulationParams& params, std::string_view type_name, Args... type_names) {
        if (HasFloat && type_name == kFloatTypeName) {
            next_type<float, Args...>(field, params, type_names...);
            return;
        }
        if (HasDouble && type_name == kDoubleTypeName) {
            next_type<double, Args...>(field, params, type_names...);
            return;
        }

        size_t n, k;
        if (parse_fixed_string(type_name, n, k, true)) {
            if (handle_fixed<0, true, Args...>(n, k, field, params, type_names...)) return;
        }
        if (parse_fixed_string(type_name, n, k, false)) {
            if (handle_fixed<0, false, Args...>(n, k, field, params, type_names...)) return;
        }

        throw invalid_argument("Could not parse type " + std::string{type_name});
    }

    template <class FloatType, class... Args>
    static void next_type(const vector<vector<char>>& field, SimulationParams& params, const Args&... type_names) {
        if constexpr (sizeof...(type_names) > 0) {
            using DefinedTypesList = TypeList<DefinedTypes...>;
            using SelectedTypesList = TypeList<SelectedTypes..., FloatType>;
            TypeSelector<DefinedTypesList, SelectedTypesList, Sizes...>::process_types(field, params, type_names...);
        } else {
            select_static_sizes<SelectedTypes..., FloatType, Sizes...>(field, params);
        }
    }

    template <size_t I, bool IsFast, class... Args>
    static bool handle_fixed(size_t n, size_t k, const vector<vector<char>>& field, SimulationParams& params, Args... type_names) {
        using FloatType = TypeAt<I, DefinedTypes...>;

        if constexpr (requires {
            { FloatType::Nval == size_t{} } -> std::same_as<bool>;
            { FloatType::Kval == size_t{} } -> std::same_as<bool>;
            { FloatType::IsFast == bool{} } -> std::same_as<bool>;
        }) {
            if constexpr (FloatType::IsFast == IsFast) {
                if (FloatType::Nval == n && FloatType::Kval == k) {
                    next_type<FloatType, Args...>(field, params, type_names...);
                    return true;
                }
            }
        }

        if constexpr (I + 1 < sizeof...(DefinedTypes)) {
            return handle_fixed<I + 1, IsFast, Args...>(n, k, field, params, type_names...);
        }

        return false;
    }
};

int main() {
    ifstream fin{};
    fin.open("../input.txt", ios_base::in);
    if (!fin) {
        cerr << "Failed to open file! Error: " << strerror(errno) << endl;
        return 1;
    }

    int N, M;
    fin >> N >> M;
    float rho_air, g;
    int rho_fluid;
    string line;
    vector<vector<char>> field(N, vector<char>(M + 1));
    fin.ignore(numeric_limits<streamsize>::max(), '\n');
    for (int i = 0; i < N; ++i) {
        getline(fin, line);
        for (int j = 0; j < M + 1; ++j) {
            field[i][j] = line[j];
            if (j == M) field[i][j] = '\000';
        }
    }

    fin >> rho_air >> rho_fluid >> g;
    fin.close();
    SimulationParams params = {rho_air, rho_fluid, g};

    TypeSelector<TypeList<TYPE_LIST>, TypeList<>, DIMENSIONS>::process_types(field, params, "FLOAT", "FIXED(32, 16)", "FAST_FIXED(32, 16)");
    return 0;
}