#pragma once

#ifdef __APPLE__
#include <mach/mach.h>
#endif

#include <unistd.h>
#include <cstring>
#include <fstream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <string_view>

#include "abort_if.hpp"
#include "basics.hpp"

namespace rcomp::utils {

template <class T>
[[maybe_unused]] inline T get_memory_block(const uint8_t* p, size_type n) {
    T v = {};
    std::memcpy(reinterpret_cast<uint8_t*>(&v), p, n);
    return v;
}
template <class T>
[[maybe_unused]] inline T get_memory_block(const uint8_t* p) {
    return get_memory_block<T>(p, sizeof(T));
}

template <class T>
[[maybe_unused]] inline void set_memory_block(uint8_t* p, T v, size_type n) {
    std::memcpy(p, reinterpret_cast<uint8_t*>(&v), n);
}
template <class T>
[[maybe_unused]] inline void set_memory_block(uint8_t* p, T v) {
    set_memory_block(p, v, sizeof(T));
}

// TODO: Make it branch-Free
// https://stackoverflow.com/questions/2274428/how-to-determine-how-many-bytes-an-integer-needs
[[maybe_unused]] inline uint8_t get_nbytes(size_type v) {
    uint8_t n = 1;
    for (; 256 <= v; v >>= 8) n++;
    return n;
}

template <class t_Rlbwt>
[[maybe_unused]] static void extend_in_reverse(const uchar_type* text, size_type size, t_Rlbwt& rlbwt) {
    for (size_type i = 1; i <= size; i++) {
        rlbwt.extend(text[size - i]);
    }
}

template <class t_Rlbwt>
[[maybe_unused]] static text_type decode_text(const t_Rlbwt& rlbwt) {
    text_type decoded;
    decoded.reserve(rlbwt.get_num_chars());
    rlbwt.decode_text([&](uchar_type c) { decoded.push_back(c); });
    return decoded;
}

template <class t_Rlbwt>
[[maybe_unused]] static std::vector<run_type> output_runs(t_Rlbwt& rlbwt) {
    std::vector<run_type> runs;
    rlbwt.output_runs([&](const run_type& rn) { runs.push_back(rn); });
    return runs;
}

template <class t_Rindex>
[[maybe_unused]] static std::vector<size_type> locate(range_type<const uchar_type*> pat, t_Rindex& rindex) {
    std::vector<size_type> sa_entries;
    rindex.locate(pat, [&](const size_type& sae) { sa_entries.push_back(sae); });
    return sa_entries;
}

template <class t_Rindex>
[[maybe_unused]] static std::vector<size_type> extract_sa_entries(t_Rindex& rindex) {
    std::vector<size_type> sa_entries;
    rindex.extract_sa_entries([&](const size_type& sae) {
        sa_entries.push_back(sae);
        DEBUG_ABORT_IF_LT(rindex.get_num_chars(), sa_entries.size());
    });
    return sa_entries;
}

template <class t_Rlbwt>
[[maybe_unused]] static size_type compute_num_runs(t_Rlbwt& rlbwt) {
    size_type num_runs = 0;
    rlbwt.output_runs([&](const run_type&) { num_runs += 1; });
    return num_runs;
}

template <class t_Rlbwt>
[[maybe_unused]] static size_type compute_max_exponent(t_Rlbwt& rlbwt) {
    size_type max_exp = 0;
    rlbwt.output_runs([&](const run_type& rn) { max_exp = std::max(max_exp, rn.exp); });
    return max_exp;
}

[[maybe_unused]] static text_type to_text(std::string_view str) {
    text_type text(str.size());
    for (size_type i = 0; i < str.size(); i++) {
        text[i] = static_cast<uchar_type>(str[i]);
    }
    return text;
}

[[maybe_unused]] static std::vector<text_type>  //
sample_subtexts(const text_type& text, size_type samples, size_type length, size_type seed = 13) {
    ABORT_IF_LT(text.size(), length);

    std::mt19937_64 gen(seed);
    std::uniform_int_distribution<size_type> dist(0, text.size() - length);

    std::vector<text_type> substrs(samples);
    for (size_type i = 0; i < samples; i++) {
        const size_type beg = dist(gen);
        substrs[i].resize(length);
        for (size_type j = 0; j < length; j++) {
            substrs[i][j] = text[beg + j];
        }
    }
    return substrs;
}

template <class T>
[[maybe_unused]] inline double get_average(const std::vector<T>& vec) {
    return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
}

[[maybe_unused]] static bool is_equivalent_runs(const std::vector<run_type>& a, const std::vector<run_type>& b) {
    if (a.size() != b.size()) {
        return false;
    }
    for (size_type i = 0; i < a.size(); i++) {
        if ((a[i].chr != b[i].chr) or (a[i].exp != b[i].exp)) {
            return false;
        }
    }
    return true;
}

template <class T>
[[maybe_unused]] static bool is_equivalent_vec(const std::vector<T>& a, const std::vector<T>& b) {
    if (a.size() != b.size()) {
        return false;
    }
    for (size_type i = 0; i < a.size(); i++) {
        if (a[i] != b[i]) {
            return false;
        }
    }
    return true;
}

[[maybe_unused]] inline double to_KiB(double bytes) {
    return bytes / 1024.0;
}

[[maybe_unused]] inline double to_MiB(double bytes) {
    return bytes / (1024.0 * 1024.0);
}

[[maybe_unused]] inline double ns_to_sec(double ns) {
    return ns / 1'000'000'000;
}

}  // namespace rcomp::utils
