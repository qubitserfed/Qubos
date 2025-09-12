#include "clifford_map.hpp"
#include "clifford_state.hpp"
#include "linear_algebra.hpp"

#include <iostream>


CliffordMap tensor(CliffordMap f, CliffordMap g) {
    CliffordMap h;

    h.state = tensor(f.state, g.state);
    h.in_wires = f.in_wires + g.in_wires;
    h.out_wires = f.out_wires + g.out_wires;

    std::vector<int> perm(h.state.n);
    // product in_wires
    for (int i = 0; i < f.in_wires; ++i)
        perm[i] = i;
    for (int i = 0; i < g.in_wires; ++i)
        perm[i + f.in_wires] = i + f.in_wires + f.out_wires;
    
    for (int i = 0; i < f.out_wires; ++i)
        perm[i + f.in_wires + g.in_wires] = i + f.in_wires;
    for (int i = 0; i < g.out_wires; ++i)
        perm[i + f.in_wires + g.in_wires + f.out_wires] = i + f.in_wires + f.out_wires + g.in_wires;

    permute(h.state, perm);
    return h;
}

CliffordMap compose(CliffordMap f, CliffordMap g) {
    if (f.out_wires != g.in_wires) {
        my_assert(0);
    }


    StabState h_state = tensor(f.state, g.state);

    std::vector<int> perm(h_state.n);
    for (int i = 0; i < f.in_wires; ++i)
        perm[i] = i;

    for (int i = 0; i < g.out_wires; ++i)
        perm[i + f.in_wires] = i + f.in_wires + f.out_wires + g.in_wires;

    for (int i = 0; i < g.in_wires; ++i) {
        perm[2 * i + f.in_wires + g.out_wires] = f.in_wires + i;
        perm[2 * i + 1 + f.in_wires + g.out_wires] = f.in_wires + f.out_wires + i;
    }

    permute(h_state, perm); h_state = normal_form(h_state);

    for (int i = 0; i < f.out_wires; ++i) {
        apply_cx(h_state, h_state.n - 2 * i - 2, h_state.n - 2 * i - 1);
        h_state = normal_form(h_state);
        apply_h(h_state, h_state.n - 2 * i - 2);
        h_state = normal_form(h_state);
    }

    for (int i = 0; i < f.out_wires; ++i) {
        pop_qubit(h_state);
        pop_qubit(h_state);
    }

    CliffordMap h;
    h.state = normal_form(h_state);
    h.in_wires = f.in_wires;
    h.out_wires = g.out_wires;
    return h;
}

CliffordMap id_map(int n) {
    CliffordMap h;
    h.state = ground_state(2 * n);
    h.in_wires = n;
    h.out_wires = n;
    for (int i = 0; i < n; ++i) {
        apply_h(h.state, i);
        apply_cx(h.state, i, i + n);
    }
    h.state = normal_form(h.state);
    
    return h;
}

CliffordMap h_map() {
    CliffordMap h = id_map(1);
    apply_h(h.state, 1);
    h.state = normal_form(h.state);
    return h;
}

CliffordMap x_map() {
    CliffordMap h = id_map(1);
    apply_x(h.state, 1);
    h.state = normal_form(h.state);
    return h;
}

CliffordMap z_map() {
    CliffordMap h = id_map(1);
    apply_z(h.state, 1);
    h.state = normal_form(h.state);
    return h;
}

CliffordMap s_map() {
    CliffordMap h = id_map(1);
    apply_s(h.state, 1);
    h.state = normal_form(h.state);
    return h;
}

CliffordMap cnot_map() {
    CliffordMap h = id_map(2);
    apply_cx(h.state, 0, 1);
    h.state = normal_form(h.state);
    return h;
}

CliffordMap cx_map() {
    CliffordMap h = id_map(2);
    apply_cx(h.state, 2, 3);
    h.state = normal_form(h.state);
    return h;
}

CliffordMap cz_map() {
    CliffordMap h = id_map(2);
    apply_cz(h.state, 2, 3);
    h.state = normal_form(h.state);
    return h;
}
