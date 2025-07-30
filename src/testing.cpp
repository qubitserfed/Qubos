#include <algorithm>
#include <fstream>
#include <iostream>
#include <functional>
#include <string>
#include <tuple>
#include <vector>

#include "linear_algebra.hpp"
#include "quantum_utilities.hpp"
#include "clifford.hpp"

struct Gate {
    std::string name;
    int i, j;
};


std::vector<std::string> gate_types = {
    "H",
    "S",
    "CZ",
    "SWAP"
};

std::vector<Gate> rand_mirror(int qubits, int depth) {
    std::vector<Gate> res;
    for (int i = 0; i < depth; ++i) {
        int gate, a, b;
        if (qubits == 1) {
            gate = rand() % 2;
            a = 0;
            b = 0;
        }
        else {
            gate = rand() % gate_types.size();
            b = rand() % (qubits - 1) + 1;
            a = rand() % b;
        }
        res.push_back({gate_types[gate], a, b});
    }

    for (int i = 0; i < depth; ++i) {
        int revi = depth - i - 1;
        if (res[revi].name == "CZ") {
            res.push_back({ "CZ", res[revi].i, res[revi].j });
        }
        else if (res[revi].name == "SWAP") {
            res.push_back({ "SWAP", res[revi].i, res[revi].j });
        }
        else if (res[revi].name == "X") {
            res.push_back({ "X", res[revi].i, res[revi].j });
        }
        else if (res[revi].name == "Z") {
            res.push_back({ "Z", res[revi].i, res[revi].j });
        }
        else if (res[revi].name == "S") {
            res.push_back({ "SD", res[revi].i, res[revi].j });
        }
        else if (res[revi].name == "H") {
            res.push_back({ "H", res[revi].i, res[revi].j });
        }
    }
    return res;
}

int main() {
    const int N = 10;
    const int depth = 1000;

    if (false) {
        srand(32313);

        std::vector<Gate> gates = rand_mirror(N, depth);
        StabState state = ground_state(N);
        std::vector<StabState> states;

        states.push_back(state);
        for (int ptr = 0; ptr < gates.size(); ++ptr) {
            auto gate = gates[ptr];
            std::cout << ptr << ": " << gate.name << " " << gate.i;
            if (gate.name == "S" || gate.name == "H" || gate.name == "X" || gate.name == "Z")
                std::cout << std::endl;
            else
                std::cout << " " << gate.j << std::endl;

            if (gate.name == "H") {
                apply_h(state, gate.i);
            }
            else if (gate.name == "CZ") {
                apply_cz(state, gate.i, gate.j);
            }
            else if (gate.name == "SWAP") {
                apply_swap(state, gate.i, gate.j);
            }
            else if (gate.name == "X") {
                apply_x(state, gate.i);
            }
            else if (gate.name == "Z") {
                apply_z(state, gate.i);
            }
            else if (gate.name == "S") {
                apply_s(state, gate.i);
            }
            else if (gate.name == "SD") {
                apply_s(state, gate.i);
                apply_s(state, gate.i);
                apply_s(state, gate.i);
            }
            state = normal_form(state);
//            print(state);
//            print_superposition(state);
            states.push_back(state);
            std::cout << "--------------------------------" << std::endl;
        }

        for (int i = 0; i < depth; ++i) {
            if (!(states[depth + i + 1] == states[depth - i])) {
                std :: cout << depth + i << std::endl;
                break;
            }
        }
        print(state);
    }
    else {
        for (int seed = 0; seed < int(1e6); ++seed) {
            srand(seed);

            if (seed % 100 == 0)
                std::cout << seed << '\n';

            std::vector<Gate> gates = rand_mirror(N, depth);
            StabState state = ground_state(N);

            for (int ptr = 0; ptr < gates.size(); ++ptr) {
                auto gate = gates[ptr];


                if (gate.name == "H") {
                    apply_h(state, gate.i);
                }
                else if (gate.name == "CZ") {
                    apply_cz(state, gate.i, gate.j);
                }
                else if (gate.name == "SWAP") {
                    apply_swap(state, gate.i, gate.j);
                }
                else if (gate.name == "X") {
                    apply_x(state, gate.i);
                }
                else if (gate.name == "Z") {
                    apply_z(state, gate.i);
                }
                else if (gate.name == "S") {
                    apply_s(state, gate.i);
                }
                else if (gate.name == "SD") {
                    apply_s(state, gate.i);
                    apply_s(state, gate.i);
                    apply_s(state, gate.i);
                }
                state = normal_form(state);
            }

            if (!is_ground(state)) {
                state = ground_state(N);
                for (int ptr = 0; ptr < gates.size(); ++ptr) {
                    auto gate = gates[ptr];
                    std::cout << ptr << ": " << gate.name << " " << gate.i;
                    if (gate.name == "S" || gate.name == "H" || gate.name == "X" || gate.name == "Z")
                        std::cout << std::endl;
                    else
                        std::cout << " " << gate.j << std::endl;

                    if (gate.name == "H") {
                        apply_h(state, gate.i);
                    }
                    else if (gate.name == "CZ") {
                        apply_cz(state, gate.i, gate.j);
                    }
                    else if (gate.name == "SWAP") {
                        apply_swap(state, gate.i, gate.j);
                    }
                    else if (gate.name == "X") {
                        apply_x(state, gate.i);
                    }
                    else if (gate.name == "Z") {
                        apply_z(state, gate.i);
                    }
                    else if (gate.name == "S") {
                        apply_s(state, gate.i);
                    }
                    state = normal_form(state);
                    print(state);
                    std::cout << "--------------------------------" << std::endl;
                }
                
                std::cout << "Not ground!" << std::endl;
                std::cout << "Seed: " <<  seed << std::endl;
                break;
            }
        }
    }

    return 0;
}
