#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>

#include "linear_algebra.hpp"
#include "quantum_utilities.hpp"
#include "clifford.hpp"


StabState ground_state(int n) {
    StabState res;

    res.n = n;
    res.phase = 0;
    res.lin_part = std::vector<int>(n, 0);
    res.quad_part = BMatrix(n, n);

    res.A = identity(n);
    res.b = BVector(n);

    return res;
}


// generate n entries of the Bell state, with indices (0, n), (1, n+1), ..., (n-1, 2n-1)
StabState bell_states(int n) {
    StabState res;

    res.n = 2 * n;
    res.phase = 0;
    res.lin_part = std::vector<int>(2 * n);
    res.quad_part = BMatrix(2 * n, 2 * n);

    res.A = BMatrix(n, 2 * n);
    res.b = BVector(n);

    for (int i = 0; i < n; ++i) {
        res.A.set(i, i, 1);
        res.A.set(i, i + n, 1);
    }

    return res;
};


// based on assumption at page 11 in my qec notebook
StabState normal_form(StabState state) {
    std::vector<int> pivots;

    to_row_echelon(state.A, state.b);
    pivots = get_pivots(state.A, range(0, state.n));

    for (int piv = 0; piv < pivots.size(); ++piv) {
        int i, j, bi;
        std::vector<int> freevars;

        j = pivots[piv];
        for (int k = 0; k < state.A.n; ++k) {
            if (state.A.get(k, j)) {
                i = k;
                break;
            }
        }

        bi = state.b.get(i);
        for (int k = j + 1; k < state.n; ++k)
            if (state.A.get(i, k))
                freevars.push_back(k);


        if (state.lin_part[j] & 1) {
            state.phase = (state.phase + 2 * bi) % 8;
            for (auto k: freevars)
                state.lin_part[k] = (state.lin_part[k] + 1 + 2 * bi) % 4;

            for (auto a: freevars)
                for (auto b: freevars)
                    state.quad_part.flip(a, b);

            state.lin_part[j] ^= 1;
        }

        if (state.lin_part[j] & 2) {
            state.phase = (state.phase + 4 * bi) % 8;
            for (auto k: freevars)
                state.lin_part[k] = (state.lin_part[k] + 2) % 4;

            state.lin_part[j] ^= 2;
        }
    }


    for (int piv = 0; piv < pivots.size(); ++piv) {
        int i, j, bi;
        std::vector<int> freevars;

        j = pivots[piv];
        for (int k = 0; k < state.A.n; ++k) {
            if (state.A.get(k, j)) {
                i = k;
                break;
            }
        }

        bi = state.b.get(i);
        for (int k = j + 1; k < state.n; ++k)
            if (state.A.get(i, k))
                freevars.push_back(k);
        
        std::vector<int> cz_targs;
        for (int k = 0; k < state.n; ++k) if (k != j) {
            if (state.quad_part.get(j, k)) {
                cz_targs.push_back(k);
                state.quad_part.set(j, k, 0);
                state.quad_part.set(k, j, 0);
            }
        }

        for (auto targ: cz_targs) {
            if (bi)
                state.lin_part[targ] = (state.lin_part[targ] + 2) % 4;

            for (auto replacement: freevars) {
                if (replacement == targ) {
                    state.lin_part[targ] = (state.lin_part[targ] + 2) % 4;
                }
                else {
                    state.quad_part.flip(targ, replacement);
                    state.quad_part.flip(replacement, targ);
                }
            }

        }
    }

    for (int i = 0; i < state.n; ++i)
        state.quad_part.set(i, i, 0);

    return state;
}

void apply_cz(StabState &state, int i, int j) {
    state.quad_part.flip(i, j);
    state.quad_part.flip(j, i);
}

void apply_cx(StabState &state, int i, int j) {
    apply_h(state, j);
    apply_cz(state, i, j);
    apply_h(state, j);
}

void apply_s(StabState &state, int i) {
    state.lin_part[i] = (state.lin_part[i] + 1) % 4;
}

void apply_swap(StabState &state, int i, int j) {
    state.A.swap_cols(i, j);
    state.quad_part.swap_cols(i, j);
    state.quad_part.swap_rows(i, j);
    std::swap(state.lin_part[i], state.lin_part[j]);
}

void apply_x(StabState &state, int i) {
    // TO REPLACE THIS BULLSHIT:
    apply_h(state, i);
    apply_z(state, i);
    apply_h(state, i);
}

void apply_z(StabState &state, int i) {
    state.lin_part[i] = (state.lin_part[i] + 2) % 4;
}

void apply_h(StabState &state, int h_target) {
    apply_swap(state, state.n - 1, h_target); // applies hadamard TO THE LAST QUBIT after conjugating with swap(h_target, state.n)

    int n = state.n;
    BVector w0(n), w1(n);

    int a = state.lin_part[n - 1];
    int a0 = (a & 1);
    int a1 = (a & 2) >> 1;

    for (int k = 0; k < n - 1; ++k)
        if (state.quad_part.get(k, n - 1))
            w0.set(k, true);
    if (a0)
        w0.set(n - 1, true);
    w1 = w0;
    w1.flip(n - 1);

    if ((state.A * versor(state.n, state.n - 1)).is_zero()) {
        int p = (a1 + int(w0.get(n - 1))) % 2;
        
        state.lin_part[n - 1] = (state.lin_part[n - 1] + 2) % 4;
        if (!a0) {
            state.A.append_row(w1);
            state.b.push_back(p);
        }
        else {
            std::vector<int> w1_support;

            for (int i = 0; i < state.n; ++i)
                if (w1.get(i))
                    w1_support.push_back(i);

            state.phase = (state.phase + 7) % 8;
            state.phase = (state.phase + 2 * p) % 8;
            for (int i = 0; i < w1_support.size(); ++i) { // full KCZ
                for (int j = i + 1; j < w1_support.size(); ++j) {
                    apply_cz(state, w1_support[i], w1_support[j]);
                    state.quad_part.flip(w1_support[i], w1_support[j]);
                    state.quad_part.flip(w1_support[j], w1_support[i]);
                }
            }
            for (auto i : w1_support)
                state.lin_part[i] = (state.lin_part[i] + 1) % 4;
            if (!p) {
                for (auto i : w1_support)
                    state.lin_part[i] = (state.lin_part[i] + 2) % 4;
            }
        }
    }
    else {
        BMatrix kerspace = find_kernel(state.A);
        std::vector<int> pivots;
        std::vector<int> pivot_rows;
        std::vector<int> I_basis;


        int  m = 0;
        for (int i = 0; i < kerspace.n; ++i) {
            if (kerspace.get(i, n - 1)) {
                kerspace.swap_rows(i, m);
                m+= 1;
            }
        }
        pivots = get_pivots(kerspace);
        for (int k = 0; k < pivots.size(); ++k) {
            int j = pivots[k];
            for (int i = 0; i < kerspace.n; ++i) {
                if (kerspace.get(i, j))  {
                    pivot_rows.push_back(i);
                    break;
                }
            }
        }

        for (int i = 0; i < kerspace.n; ++i)
            if (kerspace.row(i) * w1)
                I_basis.push_back(i);

        BVector sol = solve(state.A, state.b);

        auto apply_cz = [&](int i, int j) -> void {
            state.quad_part.flip(i, j);
            state.quad_part.flip(j, i);
            if (sol.get(i))
                state.lin_part[j] = (state.lin_part[j] + 2) % 4;
            if (sol.get(j))
                state.lin_part[i] = (state.lin_part[i] + 2) % 4;
            if (sol.get(i) && sol.get(j))
                state.phase = (state.phase + 4) % 8;
        };

        auto apply_s = [&](int i) -> void {
            state.lin_part[i] = (state.lin_part[i] + 1) % 4;
            if (sol.get(i)) {
                state.lin_part[i] = (state.lin_part[i] + 2) % 4; /// WHY DID IT WORK WITHOUT THIS?!
                state.phase = (state.phase + 2) % 8;
            }
        };

        auto apply_z = [&](int i) -> void {
            state.lin_part[i] = (state.lin_part[i] + 2) % 4;
            if (sol.get(i)) {
                state.phase = (state.phase + 4) % 8;
            }
        };

        auto amazing_gadget = [&](int t) -> void {
            apply_cz(t, n - 1);
            for (int i = 0; i < m; ++i) {
                int pi = pivots[i];
                if (pi == t)
                    apply_z(t);
                else
                    apply_cz(t, pi);
            }
        };


//        apply_z(n - 1);
        state.lin_part[n - 1] = (state.lin_part[n - 1] + 2) % 4; // this is the only time when we pull a phase gadget before shifting the affine space to the origin
        for (auto t : I_basis)
            amazing_gadget(pivots[t]);
        
        if ((sol * w1) ^ a1 ^ w1.get(n - 1)) {
            apply_z(n - 1);
            for (int i = 0; i < m; ++i)
                apply_z(pivots[i]);
        }
        
        if (a0) {
            std::vector<int> cz_domain;
            for (int i = 0; i < m; ++i)
                cz_domain.push_back(pivots[i]);
            cz_domain.push_back(n - 1);

            for (int i = 0; i < cz_domain.size(); ++i)
                for (int j = i + 1; j < cz_domain.size(); ++j)
                    apply_cz(cz_domain[i], cz_domain[j]);

            apply_s(n - 1);
            for (int i = 0; i < m; ++i)
                apply_s(pivots[i]);
        }
        
        std::tie(state.A, state.b) = affine_extension(state.A, state.b, versor(n, n - 1));
    }

    apply_swap(state, state.n - 1, h_target);
}

void print_compact(StabState state) {
    std::cout << "(";
    std::cout << state.n << ",";
    for (int i = 0; i < state.n; ++i)
        std::cout << state.lin_part[i];
    std::cout << ',';
    for (int i = 0; i < state.n; ++i)
        for (int j = 0; j < state.n; ++j)
            std::cout << state.quad_part.get(i, j);
    std::cout << ',';
    for (int i = 0; i < state.A.n; ++i)
        for (int j = 0; j < state.A.m; ++j)
            std::cout << state.A.get(i, j);
    std::cout << ',';
    for (int i = 0; i < state.b.n; ++i)
        std::cout << state.b.get(i);
    std::cout<< ")\n";
}

void print(StabState state) {
    std::cout << "Phase polynomial matrix:\n";
    for (int i = 0; i < state.n; ++i) {
        for (int j = 0; j < state.n; ++j) {
            if (i == j)
                std::cout << state.lin_part[i] << ' ';
            else
                std::cout << 2 * int(state.quad_part.get(i, j)) << ' ';
        }
        std::cout << '\n';
    }
    std::cout << "Affine pair:\n";
    std::cout << "A:\n";
    for (int i = 0; i < state.A.n; ++i) {
        for (int j = 0; j < state.A.m; ++j)
            std::cout << state.A.get(i, j) << ' ';
        std::cout << '\n';
    }
    std::cout<< "B: ";
    for (int i = 0; i < state.b.n; ++i)
        std::cout << state.b.get(i) << " ";
    std::cout << std::endl;
    std::cout << "Global phase: " << state.phase << std::endl;
    std::cout << std::endl;
}

bool is_ground(StabState state) {
    bool flag = true;
 
    flag&= state.phase == 0;
    flag&= state.A == identity(state.n);
    flag&= state.b.is_zero();
    flag&= state.lin_part == std::vector<int>(state.n, 0);
    flag&= is_zero(state.quad_part);
    return flag;
}

bool operator == (StabState u, StabState v) {
    u = normal_form(u);
    v = normal_form(v);

    bool equal = true;
    equal&= u.phase == v.phase;
    equal&= u.A == v.A;
    equal&= u.b == v.b;
    equal&= u.lin_part == v.lin_part;
    equal&= u.quad_part == v.quad_part;

    return equal;
}


void print_superposition(StabState state) {
    std::function<void(std::vector<bool>)> bkt = [&](std::vector<bool> arr) -> void {
        if (arr.size() == state.n) {
            if (state.A * BVector(arr) != state.b)
                return;

            int phase = 0;
            for (int i = 0; i < state.n; ++i)
                for (int j = i + 1; j < state.n; ++j)
                    if (arr[i] && arr[j])
                        phase+= 2 * int(state.quad_part.get(i, j));
            for (int i = 0; i < state.n; ++i)
                if (arr[i])
                    phase+= state.lin_part[i];

            phase%= 4;

            switch (phase) {
                case 0:
                    std::cout << "  +";
                    break;
                case 1:
                    std::cout << " +i";
                    break;
                case 2:
                    std::cout << "  -";
                    break;
                case 3:
                    std::cout << " -i";
                    break;
            }
            std::cout << "|";
            for (auto i: arr)
                std::cout << int(i);
            std::cout << ">";

            return;
        }

        arr.push_back(0);
        bkt(arr);
        arr.pop_back();

        arr.push_back(1);
        bkt(arr);
        arr.pop_back();
    };

    bkt(std::vector<bool>());
    std::cout << '\n';
}
