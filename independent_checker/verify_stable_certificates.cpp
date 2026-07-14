#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

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

int main(int argc, char** argv) {
    const std::string input = argc > 1 ? argv[1] : "../certificates/global_residual_certificates_n5_n30.csv";
    std::ifstream in(input);
    if (!in) {
        std::cerr << "Cannot open certificate file: " << input << "\n";
        return 2;
    }
    std::string line;
    std::getline(in, line);
    int rows = 0, certified = 0, exceptions = 0, failures = 0, max_prime = 0;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream parser(line);
        int n, p, a0, prime;
        if (!(parser >> n >> p >> a0 >> prime)) {
            ++failures;
            continue;
        }
        ++rows;
        if (prime < 0) {
            const bool expected = (n == 5 && ((p == 1 && a0 == 4) || (p == 4 && a0 == 2)));
            if (!expected) ++failures;
            ++exceptions;
            continue;
        }
        ++certified;
        max_prime = std::max(max_prime, prime);
        for (int q = 0; q < prime; ++q) {
            if (determinant_value(n, p, a0, q, prime) == 0) {
                std::cerr << "Certificate failure at n=" << n << " p=" << p << " a0=" << a0
                          << " prime=" << prime << " q=" << q << "\n";
                ++failures;
                break;
            }
        }
    }
    std::cout << "rows=" << rows
              << " certified=" << certified
              << " exceptions=" << exceptions
              << " max_prime=" << max_prime
              << " failures=" << failures << "\n";
    return (rows == 8112 && certified == 8110 && exceptions == 2 && max_prime == 107 && failures == 0) ? 0 : 1;
}
