#include <algorithm>
#include <atomic>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using std::vector;

namespace {
constexpr int kModulus = 1'000'000'007;

int mod_pow(long long a, long long exponent) {
    long long result = 1;
    while (exponent) {
        if (exponent & 1) result = result * a % kModulus;
        a = a * a % kModulus;
        exponent >>= 1;
    }
    return static_cast<int>(result);
}

int det_mod(vector<vector<int>> matrix) {
    const int n = static_cast<int>(matrix.size());
    long long determinant = 1;
    for (int column = 0; column < n; ++column) {
        int pivot = column;
        while (pivot < n && matrix[pivot][column] == 0) ++pivot;
        if (pivot == n) return 0;
        if (pivot != column) {
            std::swap(matrix[pivot], matrix[column]);
            determinant = kModulus - determinant;
        }
        const int pivot_value = matrix[column][column];
        determinant = determinant * pivot_value % kModulus;
        const int inverse = mod_pow(pivot_value, kModulus - 2);
        for (int row = column + 1; row < n; ++row) {
            const int factor = static_cast<int>((long long)matrix[row][column] * inverse % kModulus);
            if (!factor) continue;
            for (int j = column; j < n; ++j) {
                int value = matrix[row][j] - static_cast<int>((long long)factor * matrix[column][j] % kModulus);
                if (value < 0) value += kModulus;
                matrix[row][j] = value;
            }
        }
    }
    return static_cast<int>(determinant);
}

vector<int> digit_sum_coefficients(int q, int n) {
    vector<int> coefficients(1, 1);
    int current_degree = 0;
    const int N = q - 1;
    for (int power = 0; power < n; ++power) {
        vector<int> next(current_degree + N + 1, 0);
        long long window = 0;
        for (int u = 0; u <= current_degree + N; ++u) {
            if (u <= current_degree) window += coefficients[u];
            if (u - q >= 0 && u - q <= current_degree) window -= coefficients[u - q];
            window %= kModulus;
            if (window < 0) window += kModulus;
            next[u] = static_cast<int>(window);
        }
        coefficients.swap(next);
        current_degree += N;
    }
    return coefficients;
}

struct AuditResult {
    long long total = 0;
    long long predicted = 0;
    long long predicted_zero = 0;
    long long nonpredicted = 0;
    long long nonpredicted_nonzero = 0;
    long long classification_mismatch = 0;
    long long route_mismatch = 0;
};

}  // namespace

int main(int argc, char** argv) {
    const std::string output = argc > 1 ? argv[1] : "../logs/unstable_audit_summary.txt";

    vector<std::pair<int, int>> parameter_pairs;
    for (int n = 2; n <= 30; ++n)
        for (int q = 2; q <= n; ++q)
            parameter_pairs.emplace_back(n, q);

    vector<AuditResult> results(parameter_pairs.size());
    const auto start = std::chrono::steady_clock::now();

#pragma omp parallel for schedule(dynamic)
    for (int index = 0; index < static_cast<int>(parameter_pairs.size()); ++index) {
        const auto [n, q] = parameter_pairs[index];
        const int N = q - 1;
        const int r = n + 2;
        const int maximum_degree = n * N;
        const int q_to_n = mod_pow(q, n);
        const vector<int> coefficients = digit_sum_coefficients(q, n);

        vector<int> tail(maximum_degree + 2, 0);
        for (int u = maximum_degree; u >= 0; --u) {
            tail[u] = tail[u + 1] + coefficients[u];
            if (tail[u] >= kModulus) tail[u] -= kModulus;
        }
        auto tail_count = [&](int threshold) {
            if (threshold <= 0) return q_to_n;
            if (threshold > maximum_degree) return 0;
            return tail[threshold];
        };

        AuditResult audit;
        for (int t = 1; t <= r * N; ++t) {
            ++audit.total;
            const int t_star = r * N - t + 1;
            const bool self_dual = q % 2 == 0 && r % 2 == 1 && 2 * t == r * N + 1;
            const bool predicted = t < N || t_star < N || self_dual;
            if (predicted) ++audit.predicted;
            else ++audit.nonpredicted;

            vector<vector<int>> matrix(N, vector<int>(N));
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    const int u1 = t - 2 - i - j;
                    const int u2 = t - 1 - i - j;
                    const int c1 = (0 <= u1 && u1 <= maximum_degree) ? coefficients[u1] : 0;
                    const int c2 = (0 <= u2 && u2 <= maximum_degree) ? coefficients[u2] : 0;
                    int value = c1 - c2;
                    if (value < 0) value += kModulus;
                    matrix[i][j] = value;

                    long long alternate = static_cast<long long>(tail_count(t - i - j - 2))
                                        - 2LL * tail_count(t - i - j - 1)
                                        + tail_count(t - i - j);
                    alternate %= kModulus;
                    if (alternate < 0) alternate += kModulus;
                    if (alternate != value) ++audit.route_mismatch;
                }
            }

            const bool determinant_zero = det_mod(matrix) == 0;
            if (predicted && determinant_zero) ++audit.predicted_zero;
            if (!predicted && !determinant_zero) ++audit.nonpredicted_nonzero;
            if (predicted != determinant_zero) ++audit.classification_mismatch;
        }
        results[index] = audit;
    }

    AuditResult total;
    for (const AuditResult& result : results) {
        total.total += result.total;
        total.predicted += result.predicted;
        total.predicted_zero += result.predicted_zero;
        total.nonpredicted += result.nonpredicted;
        total.nonpredicted_nonzero += result.nonpredicted_nonzero;
        total.classification_mismatch += result.classification_mismatch;
        total.route_mismatch += result.route_mismatch;
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    std::cout << "total=" << total.total
              << " predicted=" << total.predicted
              << " predicted_zero=" << total.predicted_zero
              << " nonpredicted=" << total.nonpredicted
              << " nonpredicted_nonzero=" << total.nonpredicted_nonzero
              << " mismatch=" << total.classification_mismatch
              << " route_mismatch=" << total.route_mismatch
              << " seconds=" << seconds << '\n';

    std::ofstream out(output);
    if (!out) {
        std::cerr << "Cannot open output file: " << output << '\n';
        return 2;
    }
    out << "modulus=" << kModulus << '\n'
        << "n_range=2..30\n"
        << "q_range=2..n\n"
        << "total=" << total.total << '\n'
        << "predicted_singular=" << total.predicted << '\n'
        << "predicted_zero=" << total.predicted_zero << '\n'
        << "nonpredicted=" << total.nonpredicted << '\n'
        << "nonpredicted_nonzero=" << total.nonpredicted_nonzero << '\n'
        << "mismatch=" << total.classification_mismatch << '\n'
        << "route_mismatch=" << total.route_mismatch << '\n'
        << "seconds=" << seconds << '\n';

    return (total.total == 112375 && total.predicted == 8225 && total.predicted_zero == 8225 &&
            total.nonpredicted == 104150 && total.nonpredicted_nonzero == 104150 &&
            total.classification_mismatch == 0 && total.route_mismatch == 0) ? 0 : 1;
}
