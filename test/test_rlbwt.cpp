/**
 * @file test_rlbwt.cpp
 */
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest/doctest.h"
#include "nameof/nameof.hpp"
#include "tinyformat/tinyformat.h"

#include "Bwt_Naive.hpp"
#include "rlbwt_types.hpp"

#include "test_common.hpp"

using namespace rcomp;

template <class t_Rlbwt>
void test_rlbwt(const text_type& text) {
    auto rlbwt1 = Bwt_Naive();
    auto rlbwt2 = t_Rlbwt();

    utils::extend_in_reverse(text.data(), text.size(), rlbwt1);
    utils::extend_in_reverse(text.data(), text.size(), rlbwt2);

    rlbwt2.test();

    const auto runs1 = utils::output_runs(rlbwt1);
    const auto runs2 = utils::output_runs(rlbwt2);
    REQUIRE(utils::is_equivalent_runs(runs1, runs2));

    text_type decoded = utils::decode_text(rlbwt2);
    std::reverse(decoded.begin(), decoded.end());

    REQUIRE(utils::is_equivalent_vec(text, decoded));
}

static constexpr size_type SHORT_LENGTH = 20;
static constexpr size_type SHORT_NUM_SEED = 1000;
static constexpr size_type SHORT_BEG_SEED = 46;
static constexpr size_type SHORT_END_SEED = SHORT_BEG_SEED + SHORT_NUM_SEED;

static constexpr size_type LONG_LENGTH = 1000;
static constexpr size_type LONG_NUM_SEED = 20;
static constexpr size_type LONG_BEG_SEED = 46;
static constexpr size_type LONG_END_SEED = LONG_BEG_SEED + LONG_NUM_SEED;

static constexpr size_type DIV_BOUND = 7;

#ifdef LFIG_NAIVE
using rlbwt_type = rlbwt_types::lfig_naive<DIV_BOUND>::type;
#define RLBWT_NAME "lfig_naive"
#elif GLFIG_NAIVE_4
using rlbwt_type = rlbwt_types::glfig_naive<4, DIV_BOUND>::type;
#define RLBWT_NAME "glfig_naive_4"
#elif GLFIG_NAIVE_16
using rlbwt_type = rlbwt_types::glfig_naive<16, DIV_BOUND>::type;
#define RLBWT_NAME "glfig_naive_16"
#elif GLFIG_SERIALIZED_4
using rlbwt_type = rlbwt_types::glfig_serialized<4, DIV_BOUND>::type;
#define RLBWT_NAME "glfig_serialized_4"
#elif GLFIG_SERIALIZED_16
using rlbwt_type = rlbwt_types::glfig_serialized<16, DIV_BOUND>::type;
#define RLBWT_NAME "glfig_serialized_16"
#endif

TEST_CASE("Test RLBWT " RLBWT_NAME " (Tiny)") {
    const auto text = test_common::gen_tiny_text();
    test_rlbwt<rlbwt_type>(text);
}

TEST_CASE("Test RLBWT " RLBWT_NAME " (Random, Short, [a-b])") {
    for (size_type seed = SHORT_BEG_SEED; seed < SHORT_END_SEED; seed++) {
        const auto text = test_common::gen_random_text(SHORT_LENGTH, 'a', 'b', seed);
        test_rlbwt<rlbwt_type>(text);
    }
}

TEST_CASE("Test RLBWT " RLBWT_NAME " (Random, Short, [a-g])") {
    for (size_type seed = SHORT_BEG_SEED; seed < SHORT_END_SEED; seed++) {
        const auto text = test_common::gen_random_text(SHORT_LENGTH, 'a', 'g', seed);
        test_rlbwt<rlbwt_type>(text);
    }
}

TEST_CASE("Test RLBWT " RLBWT_NAME " (Random, Short, [a-z])") {
    for (size_type seed = SHORT_BEG_SEED; seed < SHORT_END_SEED; seed++) {
        const auto text = test_common::gen_random_text(SHORT_LENGTH, 'a', 'z', seed);
        test_rlbwt<rlbwt_type>(text);
    }
}

TEST_CASE("Test RLBWT " RLBWT_NAME " (Random, Long, [a-b])") {
    for (size_type seed = LONG_BEG_SEED; seed < LONG_END_SEED; seed++) {
        const auto text = test_common::gen_random_text(LONG_LENGTH, 'a', 'b', seed);
        test_rlbwt<rlbwt_type>(text);
    }
}

TEST_CASE("Test RLBWT " RLBWT_NAME " (Random, Long, [a-g])") {
    for (size_type seed = LONG_BEG_SEED; seed < LONG_END_SEED; seed++) {
        const auto text = test_common::gen_random_text(LONG_LENGTH, 'a', 'g', seed);
        test_rlbwt<rlbwt_type>(text);
    }
}

TEST_CASE("Test RLBWT " RLBWT_NAME " (Random, Long, [a-z])") {
    for (size_type seed = LONG_BEG_SEED; seed < LONG_END_SEED; seed++) {
        const auto text = test_common::gen_random_text(LONG_LENGTH, 'a', 'z', seed);
        test_rlbwt<rlbwt_type>(text);
    }
}

TEST_CASE("Test RLBWT " RLBWT_NAME " (Fibonacci, n=10)") {
    const auto text = test_common::gen_fibonacci_text(10);
    test_rlbwt<rlbwt_type>(text);
}

TEST_CASE("Test RLBWT " RLBWT_NAME " (Fibonacci, n=15)") {
    const auto text = test_common::gen_fibonacci_text(15);
    test_rlbwt<rlbwt_type>(text);
}
