/**
 * @file perf_rlbwt.cpp
 */
#include "Timer.hpp"
#include "io.hpp"
#include "rlbwt_types.hpp"
#include "utils.hpp"

#include "cmd_line_parser/parser.hpp"

using namespace rcomp;

static constexpr size_type DIV_BOUND = 7;
static constexpr size_type MILESTONE_BEG = 10000;

cmd_line_parser::parser make_parser(int argc, char** argv) {
    cmd_line_parser::parser p(argc, argv);
    p.add("input_path", "Input file path of text");
    p.add("rlbwt_type", "Rlbwt data structure type: lfig | glfig_[8|16|32|64] (default=glfig_16)", "-t", false);
    p.add("reverse_mode", "Load the text in reverse? (default=1)", "-r", false);
    p.add("enable_test", "Test the data structure? (default=1)", "-T", false);
    return p;
}

template <class t_Rlbwt>
int perf_rlbwt(const cmd_line_parser::parser& p) {
    const auto input_path = p.get<std::string>("input_path");
    const auto reverse_mode = p.get<bool>("reverse_mode", true);
    const auto enable_test = p.get<bool>("enable_test", true);

    tfm::printfln("[Input_Params]");
    tfm::reportfln("input_path:\t%s", input_path);
    tfm::reportfln("rlbwt_type:\t%s", nameof::nameof_type<t_Rlbwt>());
    tfm::reportfln("reverse_mode:\t%d", int(reverse_mode));

    const text_type text = io::load_text(input_path, reverse_mode);

    double construction_sec = 0.0;
    size_type rss_bytes = io::get_rss_in_bytes();
    size_type milestone = MILESTONE_BEG;

    t_Rlbwt rlbwt;

    Timer timer;
    timer.start();

    for (size_type i = 0; i < text.size(); i++) {
        rlbwt.extend(text[i]);
        if (rlbwt.get_num_chars() == milestone) {
            timer.stop();
            construction_sec += timer.get_elapsed_sec();
            tfm::printfln("[Progress_Report]");
            tfm::reportfln("num_chars:\t%d", rlbwt.get_num_chars());
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

    const size_type num_runs = utils::compute_num_runs(rlbwt);
    const size_type num_chars = rlbwt.get_num_chars();
    tfm::reportfln("num_runs:\t%d", num_runs);
    tfm::reportfln("num_chars:\t%d", num_chars);
    tfm::reportfln("compression_ratio:\t%g", num_runs / double(num_chars));

    const size_type alloc_bytes = rlbwt.get_memory_in_bytes();
    tfm::reportfln("alloc_memory_in_bytes:\t%d", alloc_bytes);
    tfm::reportfln("alloc_memory_in_MiB:\t%g", utils::to_MiB(alloc_bytes));
    tfm::reportfln("peak_memory_in_bytes:\t%d", rss_bytes);
    tfm::reportfln("peak_memory_in_MiB:\t%g", utils::to_MiB(rss_bytes));

    // rlbwt.show_memory_statistics();
    // rlbwt.show_detailed_statistics();
    // rlbwt.show_monitored_statistics();

    if (enable_test) {
        tfm::printfln("Testing the data structure now...");
        rlbwt.test();
        tfm::printfln("No Problem!");

        tfm::printfln("Testing the decoded text now...");
        auto decoded = utils::decode_text(rlbwt);
        if (utils::is_equivalent_vec(text, decoded)) {
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

    const auto rlbwt_type_str = p.get<std::string>("rlbwt_type", "glfig_16");

    if (rlbwt_type_str == "lfig") {
        return perf_rlbwt<rlbwt_types::lfig_naive<DIV_BOUND>::type>(p);
    } else if (rlbwt_type_str == "glfig_8") {
        return perf_rlbwt<rlbwt_types::glfig_serialized<8, DIV_BOUND>::type>(p);
    } else if (rlbwt_type_str == "glfig_16") {
        return perf_rlbwt<rlbwt_types::glfig_serialized<16, DIV_BOUND>::type>(p);
    } else if (rlbwt_type_str == "glfig_32") {
        return perf_rlbwt<rlbwt_types::glfig_serialized<32, DIV_BOUND>::type>(p);
    } else if (rlbwt_type_str == "glfig_64") {
        return perf_rlbwt<rlbwt_types::glfig_serialized<64, DIV_BOUND>::type>(p);
    }

    p.help();
    return 1;
}