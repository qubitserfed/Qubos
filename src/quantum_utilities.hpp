#pragma once

#include <algorithm>
#include <functional>
#include <vector>

#include "linear_algebra.hpp"
#include "combinatorics.hpp"
#include "quantum_utilities.hpp"


int stabilizer_weight(BVector);


void                        iterate_words_of_weight         (int, int, std::function<void(std::vector<bool>) >);
int                         bruteforce_distance0            (BMatrix);
int                         bruteforce_distance1            (BMatrix);
std::pair<int, int>         bruteforce_zx_distance0         (BMatrix);
std::pair<BMatrix, BMatrix> zx_parts                        (BMatrix);
BMatrix                     logical_operators               (BMatrix);
int                         symplectic_weight               (BVector);
bool                        is_css                          (BMatrix);

