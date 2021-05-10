/**
 * @file utils.hpp
 * @author Shunsuke Kanda (kampersanda)
 */
#pragma once

#include <cstring>
#include <fstream>
#include <map>
#include <random>
#include <sstream>

#include "abort_if.hpp"
#include "basics.hpp"

namespace rcomp::utils {

[[maybe_unused]] static std::ifstream make_ifstream(const std::string& filepath) {
    std::ifstream ifs(filepath);
    ABORT_IF(!ifs);
    return ifs;
}

[[maybe_unused]] static std::ofstream make_ofstream(const std::string& filepath) {
    std::ofstream ofs(filepath);
    ABORT_IF(!ofs);
    return ofs;
}

[[maybe_unused]] static size_type get_filesize(const std::string& filepath) {
    std::ifstream is(filepath, std::ios::binary | std::ios::ate);
    ABORT_IF(!is.good());
    const auto bytes = static_cast<size_type>(is.tellg());
    is.close();
    return bytes;
}

[[maybe_unused]] static std::vector<uchar_type> load_text(const std::string& filepath) {
    const size_type         filesize = get_filesize(filepath);
    std::vector<uchar_type> text(filesize);
    {
        auto ifs = make_ifstream(filepath);
        ifs.read(reinterpret_cast<char*>(text.data()), filesize);
    }
    return text;
}

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

template <class t_Compressor>
[[maybe_unused]] static void extend_in_reverse(const uchar_type* text, size_type size, t_Compressor& cmpr) {
    for (size_type i = 1; i <= size; i++) {
        cmpr.extend(text[size - i]);
    }
}

template <class t_Compressor>
[[maybe_unused]] static std::vector<run_type> output_runs(t_Compressor& cmpr) {
    std::vector<run_type> runs;
    cmpr.output_runs([&](const run_type& rn) { runs.push_back(rn); });
    return runs;
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

template <class t_Compressor>
[[maybe_unused]] static std::vector<uchar_type> decode(const t_Compressor& cmpr) {
    std::vector<uchar_type> decoded;
    decoded.reserve(cmpr.get_num_chars());
    cmpr.decode([&](uchar_type c) { decoded.push_back(c); });
    return decoded;
}

[[maybe_unused]] static bool is_equivalent_text(const std::vector<uchar_type>& a, const std::vector<uchar_type>& b) {
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

[[maybe_unused]] static void print_runs(const std::vector<run_type>& runs, const char* msg = "") {
    tfm::printf("%s:", msg);
    for (size_type i = 0; i < runs.size(); i++) {
        tfm::printf(" %s", runs[i]);
    }
    tfm::printfln("");
}

[[maybe_unused]] static void print_text(const std::vector<uchar_type>& text, const char* msg = "") {
    tfm::printf("%s:", msg);
    for (size_type i = 0; i < text.size(); i++) {
        tfm::printf(" %c", to_print(text[i]));
    }
    tfm::printfln("");
}

template <class t_Compressor>
[[maybe_unused]] static size_type compute_num_runs(t_Compressor& cmpr) {
    size_type num_runs = 0;
    cmpr.output_runs([&](const run_type&) { num_runs++; });
    return num_runs;
}
template <class t_Compressor>
[[maybe_unused]] static size_type compute_max_exp(t_Compressor& cmpr) {
    size_type max_exp = 0;
    cmpr.output_runs([&](const run_type& rn) { max_exp = std::max(max_exp, rn.exp); });
    return max_exp;
}

[[maybe_unused]] static std::vector<uchar_type> to_uchar_vec(const char* text) {
    std::vector<uchar_type> vec;
    while (*text) {
        vec.push_back(*text++);
    }
    return vec;
}

[[maybe_unused]] static std::vector<uchar_type> gen_tiny_text() {
    return utils::to_uchar_vec("abaababaab");
}

[[maybe_unused]] static std::vector<uchar_type> gen_random_text(size_type size, uchar_type min_ch, uchar_type max_ch,
                                                                size_type seed = 13) {
    std::mt19937                              gen(seed);
    std::uniform_int_distribution<uchar_type> dist_ch(min_ch, max_ch);

    std::vector<uchar_type> text(size);
    for (uint64_t i = 0; i < text.size(); i++) {
        text[i] = dist_ch(gen);
    }
    return text;
}

[[maybe_unused]] static std::vector<uchar_type> gen_fibonacci_text(size_type n) {
    std::string sn_1 = "a";
    std::string sn   = "ab";

    std::string tmp;
    for (size_type i = 2; i <= n; i++) {
        tmp = sn;
        sn += sn_1;
        sn_1 = tmp;
    }
    return to_uchar_vec(sn.c_str());
}

[[maybe_unused]] inline double ns_to_sec(double ns) {
    return ns / 1'000'000'000;
}

}  // namespace rcomp::utils
