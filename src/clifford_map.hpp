#pragma once

#include "clifford_state.hpp"

#include <algorithm>
#include <vector>

struct CliffordMap {
    StabState state;
    int in_wires, out_wires;
};

CliffordMap compose(CliffordMap f, CliffordMap g);
CliffordMap tensor(CliffordMap f, CliffordMap g);
CliffordMap id_map(int n);
CliffordMap h_map();
CliffordMap x_map();
CliffordMap z_map();
CliffordMap s_map();
CliffordMap cnot_map();
CliffordMap cx_map();
CliffordMap cz_map();

