#pragma once

#include "linear_algebra.hpp"

struct StabState {
    int n;

    int phase; // global phase, not used for now

    std::vector<int> lin_part;
    BMatrix quad_part; // diagonal gets ignored

    BMatrix A;
    BVector b; // define affine space Ax=b
};

bool operator == (StabState, StabState);

StabState   ground_state    (int);
StabState   bell_states     (int);
StabState   normal_form     (StabState);
void        apply_cz        (StabState &, int, int);
void        apply_cx        (StabState &, int, int);
void        apply_swap      (StabState &, int, int);
void        apply_x         (StabState &, int);
void        apply_z         (StabState &, int);
void        apply_h         (StabState &, int);
void        apply_s         (StabState &, int);
void        print           (StabState);
void        print_compact   (StabState);
bool        is_ground       (StabState);
void        print_superposition (StabState);