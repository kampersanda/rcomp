/**
 * @file perf_rindex.cpp
 */
#include "Timer.hpp"
#include "io.hpp"
#include "utils.hpp"

#include "rindex_types.hpp"

#include "cmd_line_parser/parser.hpp"

using namespace rcomp;

static constexpr size_type DIV_BOUND = 7;
static constexpr size_type MILESTONE_BEG = 10000;

static constexpr size_type NUM_QUERIES = 1000;
static constexpr size_type QUERY_LENGTH = 8;
static constexpr size_type QUERY_SEED = 13;
static constexpr size_type NUM_TRIALS = 10;

cmd_line_parser::parser make_parser(int argc, char** argv) {
    cmd_line_parser::parser p(argc, argv);
    p.add("input_path", "Input file path of text");
    p.add("rindex_type", "Rindex data structure type", "-t", false);
    p.add("reverse_mode", "Loading the text in reverse? (default=1)", "-r", false);
    p.add("enable_test", "Testing the data structure? (default=1)", "-T", false);
    return p;
}

template <class t_Rindex>
int perf_rindex(const cmd_line_parser::parser& p) {
    const auto input_path = p.get<std::string>("input_path");
    const auto reverse_mode = p.get<bool>("reverse_mode", true);
    const auto enable_test = p.get<bool>("enable_test", true);

    tfm::printfln("[Input_Params]");
    tfm::reportfln("input_path:\t%s", input_path);
    tfm::reportfln("rindex_type:\t%s", nameof::nameof_type<t_Rindex>());
    tfm::reportfln("reverse_mode:\t%d", int(reverse_mode));

    Timer timer;
    t_Rindex rindex;

    const text_type text = io::load_text(input_path, reverse_mode);

    /**
     *  Construction Phase
     */
    double construction_sec = 0.0;
    size_type rss_bytes = io::get_rss_in_bytes();
    size_type milestone = MILESTONE_BEG;

    timer.start();
    for (size_type i = 0; i < text.size(); i++) {
        rindex.extend(text[i]);
        if (rindex.get_num_chars() == milestone) {
            timer.stop();
            construction_sec += timer.get_elapsed_sec();
            tfm::printfln("[Progress_Report]");
            tfm::reportfln("num_chars:\t%d", rindex.get_num_chars());
            tfm::reportfln("construction_sec:\t%g", construction_sec);
            milestone *= 10;
            timer.start();
        }
    }
    timer.stop();

    construction_sec += timer.get_elapsed_sec();
    rss_bytes = io::get_rss_in_bytes() - rss_bytes;

    tfm::printfln("[Final_Report]");
    tfm::reportfln("construction_sec:\t%g", construction_sec);

    const size_type num_runs = utils::compute_num_runs(rindex);
    const size_type num_chars = rindex.get_num_chars();
    tfm::reportfln("num_runs:\t%d", num_runs);
    tfm::reportfln("num_chars:\t%d", num_chars);
    tfm::reportfln("compression_ratio:\t%g", num_runs / double(num_chars));

    const size_type alloc_bytes = rindex.get_memory_in_bytes();
    tfm::reportfln("alloc_memory_in_bytes:\t%d", alloc_bytes);
    tfm::reportfln("alloc_memory_in_MiB:\t%g", utils::to_MiB(alloc_bytes));
    tfm::reportfln("peak_memory_in_bytes:\t%d", rss_bytes);
    tfm::reportfln("peak_memory_in_MiB:\t%g", utils::to_MiB(rss_bytes));

    // rindex.show_memory_statistics();
    // rindex.show_detailed_statistics();

    /**
     *  Search Phase
     */
    tfm::printfln("[Search_Settings]");
    tfm::reportfln("num_trials:\t%s", NUM_TRIALS);
    tfm::reportfln("num_queries:\t%s", NUM_QUERIES);
    tfm::reportfln("query_length:\t%s", QUERY_LENGTH);
    tfm::reportfln("query_seed:\t%s", QUERY_SEED);

    const std::vector<text_type> queries = utils::sample_subtexts(text, NUM_QUERIES, QUERY_LENGTH, QUERY_SEED);

    {
        tfm::printfln("Warming up now...");
        size_type dummy = 0;
        for (size_type i = 0; i < queries.size(); i++) {
            rindex.locate(make_range(queries[i]), [&](size_type) { dummy += 1; });
        }
        tfm::printfln("dummy:\t%d", dummy);
    }

    size_type occ = 0;
    std::vector<double> times(NUM_TRIALS);  // in micro sec

    // Count query
    {
        for (size_type r = 0; r < NUM_TRIALS; r++) {
            timer.start();
            for (size_type i = 0; i < queries.size(); i++) {
                occ += rindex.count(make_range(queries[i]));
            }
            timer.stop();
            times[r] = timer.get_elapsed_us();
        }

        const double occ_per_query = double(occ / NUM_TRIALS) / queries.size();
        const double count_us_per_query = utils::get_average(times) / queries.size();

        tfm::printfln("[Count_Query]");
        tfm::reportfln("occ_per_query:\t%g", occ_per_query);
        tfm::reportfln("microsec_per_query:\t%g", count_us_per_query);
    }

    // Locate query
    {
        for (size_type r = 0; r < NUM_TRIALS; r++) {
            timer.start();
            for (size_type i = 0; i < queries.size(); i++) {
                rindex.locate(make_range(queries[i]), [&](size_type) { occ += 1; });
            }
            timer.stop();
            times[r] = timer.get_elapsed_us();
        }

        const double occ_per_query = double(occ / NUM_TRIALS) / queries.size();
        const double locate_us_per_query = utils::get_average(times) / queries.size();
        const double locate_us_per_occ = locate_us_per_query / occ_per_query;

        tfm::printfln("[Locate_Query]");
        tfm::reportfln("occ_per_query:\t%g", occ_per_query);
        tfm::reportfln("microsec_per_query:\t%g", locate_us_per_query);
        tfm::reportfln("microsec_per_occ:\t%g", locate_us_per_occ);
    }

    if (enable_test) {
        tfm::printfln("Testing the data structure now...");
        rindex.test();
        tfm::printfln("No Problem!");

        tfm::printfln("Testing the decoded text now...");
        auto dec_text = utils::decode_text(rindex);
        if (utils::is_equivalent_vec(text, dec_text)) {
            tfm::printfln("No Problem!");
        } else {
            tfm::errorfln("What's up!?");
        }
    }

    return 0;
}

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
#ifndef NDEBUG
    tfm::warnfln("The code is running in debug mode.");
#endif

    auto p = make_parser(argc, argv);
    if (!p.parse()) {
        return 1;
    }

    const auto rindex_type_str = p.get<std::string>("rindex_type", "lfig_naive");

    if (rindex_type_str == "lfig_naive") {
        return perf_rindex<rindex_types::lfig_naive<DIV_BOUND>::type>(p);
    }

    p.help();
    return 1;
}