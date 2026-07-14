#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using std::vector;

static int mod_pow(long long a, long long e, int p) {
    long long r = 1 % p;
    a %= p;
    while (e) {
        if (e & 1) r = r * a % p;
        a = a * a % p;
        e >>= 1;
    }
    return static_cast<int>(r);
}

static int det_mod(vector<vector<int>> a, int p) {
    const int n = static_cast<int>(a.size());
    long long det = 1;
    for (int i = 0; i < n; ++i) {
        int pivot = i;
        while (pivot < n && a[pivot][i] == 0) ++pivot;
        if (pivot == n) return 0;
        if (pivot != i) {
            std::swap(a[pivot], a[i]);
            det = (p - det) % p;
        }
        const int v = a[i][i];
        det = det * v % p;
        const int inv = mod_pow(v, p - 2, p);
        for (int row = i + 1; row < n; ++row) {
            const int factor = static_cast<int>((long long)a[row][i] * inv % p);
            if (!factor) continue;
            for (int col = i; col < n; ++col) {
                int x = a[row][col] - static_cast<int>((long long)factor * a[i][col] % p);
                if (x < 0) x += p;
                a[row][col] = x;
            }
        }
    }
    return static_cast<int>(det);
}

static int determinant_value(int n, int p_index, int a0, int q, int modulus) {
    const int d = n - 1;
    vector<vector<int>> matrix(d, vector<int>(d));
    int column = 0;
    for (int a = 1; a <= n; ++a) {
        if (a == a0) continue;
        const int omitted = (a < a0) ? p_index - 1 : p_index;
        vector<int> polynomial(1, 1);
        for (int j = 0; j <= n - 1; ++j) {
            if (j == omitted) continue;
            int constant = static_cast<int>((a - (long long)j * q) % modulus);
            if (constant < 0) constant += modulus;
            vector<int> next(polynomial.size() + 1, 0);
            for (std::size_t k = 0; k < polynomial.size(); ++k) {
                next[k] = static_cast<int>((next[k] + (long long)polynomial[k] * constant) % modulus);
                next[k + 1] += polynomial[k];
                if (next[k + 1] >= modulus) next[k + 1] -= modulus;
            }
            polynomial.swap(next);
        }
        for (int k = 1; k <= d; ++k) matrix[k - 1][column] = polynomial[k];
        ++column;
    }
    return det_mod(matrix, modulus);
}

struct Family { int n, p, a0; };

int main(int argc, char** argv) {
    const std::string output = argc > 1 ? argv[1] : "../certificates/global_residual_certificates_n5_n30.csv";
    const vector<int> primes = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107};
    vector<Family> families;
    for (int n = 5; n <= 30; ++n)
        for (int p = 1; p <= n - 1; ++p)
            for (int a0 = 2; a0 <= n - 1; ++a0)
                families.push_back({n,p,a0});

    vector<int> certificate(families.size(), -1);
    const auto start = std::chrono::steady_clock::now();
#pragma omp parallel for schedule(dynamic)
    for (int index = 0; index < static_cast<int>(families.size()); ++index) {
        const Family f = families[index];
        for (int prime : primes) {
            bool no_root = true;
            for (int q = 0; q < prime; ++q) {
                if (determinant_value(f.n, f.p, f.a0, q, prime) == 0) {
                    no_root = false;
                    break;
                }
            }
            if (no_root) {
                certificate[index] = prime;
                break;
            }
        }
    }

    std::ofstream out(output);
    if (!out) {
        std::cerr << "Cannot open output file: " << output << "\n";
        return 2;
    }
    out << "n,p,a0,prime\n";
    int certified = 0, exceptional = 0, max_prime = 0;
    std::map<int,int> distribution;
    for (std::size_t i = 0; i < families.size(); ++i) {
        out << families[i].n << ',' << families[i].p << ',' << families[i].a0 << ',' << certificate[i] << '\n';
        if (certificate[i] > 0) {
            ++certified;
            max_prime = std::max(max_prime, certificate[i]);
            ++distribution[certificate[i]];
        } else {
            ++exceptional;
        }
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << "families=" << families.size()
              << " certified=" << certified
              << " exceptional=" << exceptional
              << " max_prime=" << max_prime
              << " seconds=" << seconds << "\n";
    for (const auto& [prime,count] : distribution) std::cout << prime << ':' << count << ' ';
    std::cout << '\n';
    return (certified == 8110 && exceptional == 2 && max_prime == 107) ? 0 : 1;
}
