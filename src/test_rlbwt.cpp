/**
 * @file test_rlbwt.cpp
 * @author Shunsuke Kanda (kampersanda)
 */
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest/doctest.h"
#include "nameof/nameof.hpp"
#include "tinyformat/tinyformat.h"

#include "Compressor_Naive.hpp"
#include "compressors.hpp"

#include "utils.hpp"

using namespace rcomp;

template <class t_Compressor>
void test_rlbwt(const std::vector<uchar_type>& text) {
    auto cmpr1 = Compressor_Naive();
    auto cmpr2 = t_Compressor();

    utils::extend_in_reverse(text.data(), text.size(), cmpr1);
    utils::extend_in_reverse(text.data(), text.size(), cmpr2);

    cmpr2.test_all();

    const auto runs1 = utils::output_runs(cmpr1);
    const auto runs2 = utils::output_runs(cmpr2);
    REQUIRE(utils::is_equivalent_runs(runs1, runs2));

    std::vector<uchar_type> decoded = utils::decode(cmpr2);
    REQUIRE(utils::is_equivalent_text(text, decoded));
}

static constexpr size_type SHORT_LENGTH   = 20;
static constexpr size_type SHORT_NUM_SEED = 1000;
static constexpr size_type SHORT_BEG_SEED = 46;
static constexpr size_type SHORT_END_SEED = SHORT_BEG_SEED + SHORT_NUM_SEED;

static constexpr size_type LONG_LENGTH   = 1000;
static constexpr size_type LONG_NUM_SEED = 20;
static constexpr size_type LONG_BEG_SEED = 46;
static constexpr size_type LONG_END_SEED = LONG_BEG_SEED + LONG_NUM_SEED;

static constexpr size_type DIV_BOUND = 7;

#ifdef LFIG_NAIVE
using compressor_type = compressors::lfig_naive<DIV_BOUND>::compressor_type;
#define COMPRESSOR_NAME "lfig_naive"
#elif GLFIG_NAIVE_4
using compressor_type = compressors::glfig_naive<4, DIV_BOUND, true, true>::compressor_type;
#define COMPRESSOR_NAME "glfig_naive_4"
#elif GLFIG_NAIVE_16
using compressor_type = compressors::glfig_naive<16, DIV_BOUND, true, true>::compressor_type;
#define COMPRESSOR_NAME "glfig_naive_16"
#elif GLFIG_SERIALIZED_4
using compressor_type = compressors::glfig_serialized<4, DIV_BOUND, true, true>::compressor_type;
#define COMPRESSOR_NAME "glfig_serialized_4"
#elif GLFIG_SERIALIZED_16
using compressor_type = compressors::glfig_serialized<16, DIV_BOUND, true, true>::compressor_type;
#define COMPRESSOR_NAME "glfig_serialized_16"
#endif

TEST_CASE("Testing " COMPRESSOR_NAME " (Tiny)") {
    const auto text = utils::gen_tiny_text();
    test_rlbwt<compressor_type>(text);
}

TEST_CASE("Testing " COMPRESSOR_NAME " (Random, Short, [a-b])") {
    for (size_type seed = SHORT_BEG_SEED; seed < SHORT_END_SEED; seed++) {
        // tfm::printfln("seed=%d", seed);
        const auto text = utils::gen_random_text(SHORT_LENGTH, 'a', 'b', seed);
        test_rlbwt<compressor_type>(text);
    }
}

TEST_CASE("Testing " COMPRESSOR_NAME " (Random, Short, [a-g])") {
    for (size_type seed = SHORT_BEG_SEED; seed < SHORT_END_SEED; seed++) {
        // tfm::printfln("seed=%d", seed);
        const auto text = utils::gen_random_text(SHORT_LENGTH, 'a', 'g', seed);
        test_rlbwt<compressor_type>(text);
    }
}

TEST_CASE("Testing " COMPRESSOR_NAME " (Random, Short, [a-z])") {
    for (size_type seed = SHORT_BEG_SEED; seed < SHORT_END_SEED; seed++) {
        // tfm::printfln("seed=%d", seed);
        const auto text = utils::gen_random_text(SHORT_LENGTH, 'a', 'z', seed);
        test_rlbwt<compressor_type>(text);
    }
}

TEST_CASE("Testing " COMPRESSOR_NAME " (Random, Long, [a-b])") {
    for (size_type seed = LONG_BEG_SEED; seed < LONG_END_SEED; seed++) {
        const auto text = utils::gen_random_text(LONG_LENGTH, 'a', 'b', seed);
        test_rlbwt<compressor_type>(text);
    }
}

TEST_CASE("Testing " COMPRESSOR_NAME " (Random, Long, [a-g])") {
    for (size_type seed = LONG_BEG_SEED; seed < LONG_END_SEED; seed++) {
        const auto text = utils::gen_random_text(LONG_LENGTH, 'a', 'g', seed);
        test_rlbwt<compressor_type>(text);
    }
}

TEST_CASE("Testing " COMPRESSOR_NAME " (Random, Long, [a-z])") {
    for (size_type seed = LONG_BEG_SEED; seed < LONG_END_SEED; seed++) {
        const auto text = utils::gen_random_text(LONG_LENGTH, 'a', 'z', seed);
        test_rlbwt<compressor_type>(text);
    }
}

TEST_CASE("Testing " COMPRESSOR_NAME " (Fibonacci, n=10)") {
    const auto text = utils::gen_fibonacci_text(10);
    test_rlbwt<compressor_type>(text);
}

TEST_CASE("Testing " COMPRESSOR_NAME " (Fibonacci, n=15)") {
    const auto text = utils::gen_fibonacci_text(15);
    test_rlbwt<compressor_type>(text);
}
