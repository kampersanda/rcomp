/**
 * @file test_common.hpp
 */
#pragma once

#include "utils.hpp"

namespace rcomp::test_common {

[[maybe_unused]] static void print_runs(const std::vector<run_type>& runs, const char* msg = "") {
    tfm::printf("%s:", msg);
    for (size_type i = 0; i < runs.size(); i++) {
        tfm::printf(" %s", runs[i]);
    }
    tfm::printfln("");
}

[[maybe_unused]] static void print_text(const text_type& text, const char* msg = "") {
    tfm::printf("%s:", msg);
    for (size_type i = 0; i < text.size(); i++) {
        tfm::printf(" %c", to_print(text[i]));
    }
    tfm::printfln("");
}

template <class T>
[[maybe_unused]] static void print_ints(const std::vector<T>& text, const char* msg = "") {
    tfm::printf("%s:", msg);
    for (size_type i = 0; i < text.size(); i++) {
        tfm::printf(" %d", text[i]);
    }
    tfm::printfln("");
}

[[maybe_unused]] static text_type gen_tiny_text() {
    return utils::to_text("abaababaab");
}

[[maybe_unused]] static text_type gen_random_text(size_type length, uchar_type min_ch, uchar_type max_ch,
                                                  size_type seed = 13) {
    std::mt19937 gen(seed);
    std::uniform_int_distribution<uchar_type> dist(min_ch, max_ch);

    text_type text(length);
    for (uint64_t i = 0; i < text.size(); i++) {
        text[i] = dist(gen);
    }
    return text;
}

[[maybe_unused]] static auto gen_random_patterns(size_type num, size_type min_length, size_type max_length,
                                                 uchar_type min_ch, uchar_type max_ch, size_type seed = 13) {
    std::mt19937 gen(seed);
    std::uniform_int_distribution<size_type> dist(min_length, max_length);

    std::vector<text_type> patterns(num);
    for (uint64_t i = 0; i < patterns.size(); i++) {
        patterns[i] = gen_random_text(dist(gen), min_ch, max_ch, gen());
    }
    return patterns;
}

[[maybe_unused]] static text_type gen_fibonacci_text(size_type n) {
    std::string sn_1 = "a";
    std::string sn = "ab";

    std::string tmp;
    for (size_type i = 2; i <= n; i++) {
        tmp = sn;
        sn += sn_1;
        sn_1 = tmp;
    }
    return utils::to_text(sn);
}

}  // namespace rcomp::test_common