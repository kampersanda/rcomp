/**
 * @file test_rindex.cpp
 */
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest/doctest.h"
#include "nameof/nameof.hpp"
#include "tinyformat/tinyformat.h"

#include "Bwt_Naive.hpp"
#include "rindex_types.hpp"

#include "test_common.hpp"

using namespace rcomp;

template <class t_Rindex>
void test_rindex(const text_type& text, const std::vector<text_type>& queries) {
    auto rindex1 = Bwt_Naive();
    auto rindex2 = t_Rindex();

    utils::extend_in_reverse(text.data(), text.size(), rindex1);
    utils::extend_in_reverse(text.data(), text.size(), rindex2);

    rindex2.test();

    const auto runs1 = utils::output_runs(rindex1);
    const auto runs2 = utils::output_runs(rindex2);
    REQUIRE(utils::is_equivalent_runs(runs1, runs2));

    auto decoded = utils::decode_text(rindex2);
    std::reverse(decoded.begin(), decoded.end());
    REQUIRE(utils::is_equivalent_vec(text, decoded));

    auto sa_entries1 = utils::extract_sa_entries(rindex1);
    auto sa_entries2 = utils::extract_sa_entries(rindex2);

    REQUIRE(utils::is_equivalent_vec(sa_entries1, sa_entries2));

    for (const auto& q : queries) {
        sa_entries1 = utils::locate(make_range(q), rindex1);
        sa_entries2 = utils::locate(make_range(q), rindex2);
        REQUIRE(utils::is_equivalent_vec(sa_entries1, sa_entries2));

        const size_type occ = rindex2.count(make_range(q));
        REQUIRE_EQ(sa_entries1.size(), occ);
    }
}

static constexpr size_type SHORT_LENGTH = 20;
static constexpr size_type SHORT_NUM_SEED = 1000;
static constexpr size_type SHORT_BEG_SEED = 46;
static constexpr size_type SHORT_END_SEED = SHORT_BEG_SEED + SHORT_NUM_SEED;
static constexpr size_type SHORT_PAT_MIN_LENGTH = 2;
static constexpr size_type SHORT_PAT_MAX_LENGTH = 6;

static constexpr size_type LONG_LENGTH = 1000;
static constexpr size_type LONG_NUM_SEED = 20;
static constexpr size_type LONG_BEG_SEED = 46;
static constexpr size_type LONG_END_SEED = LONG_BEG_SEED + LONG_NUM_SEED;
static constexpr size_type LONG_PAT_MIN_LENGTH = 2;
static constexpr size_type LONG_PAT_MAX_LENGTH = 20;

static constexpr size_type DIV_BOUND = 7;

#ifdef LFIG_NAIVE
using rindex_type = rindex_types::lfig_naive<DIV_BOUND>::type;
#define RINDEX_NAME "lfig_naive"
#elif GLFIG_NAIVE_4
using rindex_type = rindex_types::glfig_naive<4, DIV_BOUND>::type;
#define RINDEX_NAME "glfig_naive_4"
#elif GLFIG_NAIVE_16
using rindex_type = rindex_types::glfig_naive<16, DIV_BOUND>::type;
#define RINDEX_NAME "glfig_naive_16"
#elif GLFIG_SERIALIZED_4
using rindex_type = rindex_types::glfig_serialized<4, DIV_BOUND>::type;
#define RINDEX_NAME "glfig_serialized_4"
#elif GLFIG_SERIALIZED_16
using rindex_type = rindex_types::glfig_serialized<16, DIV_BOUND>::type;
#define RINDEX_NAME "glfig_serialized_16"
#endif

auto gen_short_random_queries(uchar_type min_ch, uchar_type max_ch) {
    return test_common::gen_random_patterns(10, SHORT_PAT_MIN_LENGTH, SHORT_PAT_MAX_LENGTH, min_ch, max_ch,
                                            SHORT_END_SEED);
}

auto gen_long_random_queries(uchar_type min_ch, uchar_type max_ch) {
    return test_common::gen_random_patterns(10, LONG_PAT_MIN_LENGTH, LONG_PAT_MAX_LENGTH, min_ch, max_ch,
                                            LONG_END_SEED);
}

TEST_CASE("Test RINDEX " RINDEX_NAME " (Tiny)") {
    std::vector<text_type> queries = {
        utils::to_text("aba"), utils::to_text("baa"),  utils::to_text("abab"),
        utils::to_text("a"),   utils::to_text("aaaa"), utils::to_text("abaababaab"),
    };
    const auto text = test_common::gen_tiny_text();
    test_rindex<rindex_type>(text, queries);
}

TEST_CASE("Test RINDEX " RINDEX_NAME " (Random, Short, [a-b])") {
    const auto queries = gen_short_random_queries('a', 'b');
    for (size_type seed = SHORT_BEG_SEED; seed < SHORT_END_SEED; seed++) {
        // tfm::printfln("seed = %d", seed);
        const auto text = test_common::gen_random_text(SHORT_LENGTH, 'a', 'b', seed);
        test_rindex<rindex_type>(text, queries);
    }
}

TEST_CASE("Test RINDEX " RINDEX_NAME " (Random, Short, [a-g])") {
    const auto queries = gen_short_random_queries('a', 'g');
    for (size_type seed = SHORT_BEG_SEED; seed < SHORT_END_SEED; seed++) {
        // tfm::printfln("seed = %d", seed);
        const auto text = test_common::gen_random_text(SHORT_LENGTH, 'a', 'g', seed);
        test_rindex<rindex_type>(text, queries);
    }
}

TEST_CASE("Test RINDEX " RINDEX_NAME " (Random, Short, [a-z])") {
    const auto queries = gen_short_random_queries('a', 'z');
    for (size_type seed = SHORT_BEG_SEED; seed < SHORT_END_SEED; seed++) {
        // tfm::printfln("seed = %d", seed);
        const auto text = test_common::gen_random_text(SHORT_LENGTH, 'a', 'z', seed);
        test_rindex<rindex_type>(text, queries);
    }
}

TEST_CASE("Test RINDEX " RINDEX_NAME " (Random, Long, [a-b])") {
    const auto queries = gen_long_random_queries('a', 'b');
    for (size_type seed = LONG_BEG_SEED; seed < LONG_END_SEED; seed++) {
        const auto text = test_common::gen_random_text(LONG_LENGTH, 'a', 'b', seed);
        test_rindex<rindex_type>(text, queries);
    }
}

TEST_CASE("Test RINDEX " RINDEX_NAME " (Random, Long, [a-g])") {
    const auto queries = gen_long_random_queries('a', 'g');
    for (size_type seed = LONG_BEG_SEED; seed < LONG_END_SEED; seed++) {
        const auto text = test_common::gen_random_text(LONG_LENGTH, 'a', 'g', seed);
        test_rindex<rindex_type>(text, queries);
    }
}

TEST_CASE("Test RINDEX " RINDEX_NAME " (Random, Long, [a-z])") {
    const auto queries = gen_long_random_queries('a', 'z');
    for (size_type seed = LONG_BEG_SEED; seed < LONG_END_SEED; seed++) {
        const auto text = test_common::gen_random_text(LONG_LENGTH, 'a', 'z', seed);
        test_rindex<rindex_type>(text, queries);
    }
}

TEST_CASE("Test RINDEX " RINDEX_NAME " (Fibonacci, n=10)") {
    const auto queries = gen_long_random_queries('a', 'b');
    const auto text = test_common::gen_fibonacci_text(10);
    test_rindex<rindex_type>(text, queries);
}

TEST_CASE("Test RINDEX " RINDEX_NAME " (Fibonacci, n=15)") {
    const auto queries = gen_long_random_queries('a', 'b');
    const auto text = test_common::gen_fibonacci_text(15);
    test_rindex<rindex_type>(text, queries);
}
