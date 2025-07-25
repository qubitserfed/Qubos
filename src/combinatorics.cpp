#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>

#include "combinatorics.hpp"


// input: two integers a and b
// output: the list [a, a+1, ..., b-1]
std::vector<int> range(int a, int b) {
    std::vector<int> res;
    for (int i = a; i < b; ++i)
        res.push_back(i);
    return res;
}

// input:  a number n and a function whose sole argument is a vector of bools
// effect: the function gets called on all length n bitstrings in lexicographical order
void partitions(int n, std::function<void(std::vector<bool>)> f) {
    std::vector<bool> stack;
    stack.reserve(n);

    std::function<void()> bkt = [&]() {
        if (stack.size() == n) {
            f(stack);
            return;
        }

        stack.push_back(0);
        bkt();
        stack.pop_back();

        stack.push_back(1);
        bkt();
        stack.pop_back();
    };

    bkt();
}

// input:  a number n and a function whose sole argument is a vector of bools
// effect: the function gets called on all length n bitstrings of Hamming weight k in lexicographical order
void combinations(int n, int k, std::function<void(std::vector<bool>) > f) {
    std::vector<bool> stack;
    stack.reserve(n);

    std::function<void(int)> bkt = [&](int weight) {
        if (stack.size() == n) {
            f(stack);
            return;
        }

        if (weight <= n - (int(stack.size()) + 1)) {
            stack.push_back(0);
            bkt(weight);
            stack.pop_back();
        }

        if (weight > 0) {
            stack.push_back(1);
            bkt(weight - 1);
            stack.pop_back();
        }
    };

    bkt(k);
}

bool parallel_combinations(int n, int k, std::function<bool(BVector &)> f, int no_threads) {
    std::vector<u64> delimiters;
    std::vector<std::future<bool>> threads;

    u64 no_combinations = 1;
    for (int i = 1; i <= k; ++i)
        no_combinations = no_combinations * (n - i + 1) / i;

    for (int i = 0; i <= no_threads; ++i)
        delimiters.push_back(i * no_combinations / no_threads);

    ith_lexicographic_permutation(1, 1, 0);
    for (int i = 1; i < delimiters.size(); ++i) {
        int iv_start = delimiters[i - 1];
        int iv_end = delimiters[i];

        threads.push_back(std::async(std::launch::async, [n, k, iv_start, iv_end, f] () -> bool {
            BVector comb = ith_lexicographic_permutation(n, k, iv_start);
            for (u64 it = 0; it < iv_end - iv_start; ++it) {
                if (f(comb))
                    return true;
                advance_iterator(comb);
            }
            return false;
        }));
    }

    for (std::future<bool> &th: threads)
        if (th.get())
            return true;
    return false;
}

void symplectic_combinations(int n, int k, std::function<void(std::vector<bool> &) > f) {
    std::vector<bool> stack;
    stack.reserve(2 * n);

    std::function<void(int)> bkt = [&](int weight) {
        if (stack.size() / 2 == n) {
            f(stack);
            return;
        }

        if (weight <= n - (int(stack.size()) / 2 + 1)) {
            stack.push_back(0);
            stack.push_back(0);
            bkt(weight);
            stack.pop_back();
            stack.pop_back();
        }

        if (weight > 0) {
            stack.push_back(1);
            stack.push_back(0);
            bkt(weight - 1);
            stack.pop_back();
            stack.pop_back();

            stack.push_back(0);
            stack.push_back(1);
            bkt(weight - 1);
            stack.pop_back();
            stack.pop_back();

            stack.push_back(1);
            stack.push_back(1);
            bkt(weight - 1);
            stack.pop_back();
            stack.pop_back();
        }
    };

    bkt(k);
}

BVector ith_lexicographic_permutation(int n, int k, u64 num) {
    const u64 infty = 0x3f3f3f3f3f3f3f3fULL;

    static bool initialized = false;
    static u64 c[200][200];

    if (!initialized) {
        initialized = true;
        for (int i = 0; i < 200; ++i)
            c[i][0] = 1;
        for (int i = 1; i < 200; ++i)
            for (int j = 1; j <= i; ++j)
                c[i][j] = std::min(infty, c[i - 1][j] + c[i - 1][j - 1]);
    }

    std::vector<bool> res;
    for (int i = 0; i < n; ++i) { // VERIFY
        if (c[n - i - 1][k] > num) {
            res.push_back(0);
        }
        else {
            num-= c[n - i - 1][k];
            k-= 1;
            res.push_back(1);
        }
    }
    my_assert(k == 0);

    return BVector(res);
}

bool advance_iterator(BVector &it) {
    static bool initialized = false;
    static u64 first_bits[65];
    static u64 last_bits[65];

    if (!initialized) {
        first_bits[0] = last_bits[0] = 0;
        initialized = true;

        for (int bit = 0; bit < 64; ++bit) {
            first_bits[bit + 1] = (1ULL << bit) | first_bits[bit];
            last_bits[bit + 1] = (1ULL << 63 - bit) | last_bits[bit];
        }
    }

    /*
    Algorithm:
    1) find the last one bit b which doesn't have a one immediately before it
    2) swap bits b and b - 1
    3) move all the 1 bits before b to the end of the string 
    */

    // !!! All the implemented steps can be further optimized !!!
    // 1)
    const int n = it.n;

    int bit = n, one_count = 0;
    while (bit > 0) {
        if (it.get(bit)) {
            one_count+= 1;
        }
        else {
            if (one_count > 0) {
                bit+= 1;
                break;
            }
        }
        bit-= 1;
    }
    if (bit == 0)
        return false;

    my_assert(it.get(bit));
    // 2)
    it.set(bit - 1, true);
    it.set(bit, false);

    // 3)
    for (int i = bit; i < n; ++i)
        it.set(i, false);
    for (int i = n - 1; i >= n - one_count + 1; --i)
        it.set(i, true);

    return true;
}
