/**
 * @file perf_rlbwt.cpp
 * @author Shunsuke Kanda (kampersanda)
 */
#include "TextLoader.hpp"
#include "Timer.hpp"
#include "compressors.hpp"
#include "utils.hpp"

#include "cmd_line_parser/parser.hpp"

using namespace rcomp;

static constexpr size_type DIV_BOUND = 7;

cmd_line_parser::parser make_parser(int argc, char** argv) {
    cmd_line_parser::parser p(argc, argv);
    p.add("input_path", "input file path of text");
    p.add("compressor", "compressor type", "-c", false);
    p.add("is_in_online", "load the text in online? (default=0)", "-l", false);
    p.add("enable_progrep", "report the progress? (default=1)", "-p", false);
    p.add("enable_test", "test the data structure? (default=1)", "-t", false);
    p.add("enable_check", "check the decode? (default=1)", "-t", false);
    return p;
}

template <class t_Compressor, class t_TextLoader>
int perf_rlbwt(const cmd_line_parser::parser& p) {
    const auto input_path     = p.get<std::string>("input_path");
    const auto enable_progrep = p.get<bool>("enable_progrep", true);
    const auto enable_test    = p.get<bool>("enable_test", true);
    const auto enable_check   = p.get<bool>("enable_check", true);

    tfm::printfln("[Basic]");
    tfm::reportfln("input_path:\t%s", input_path);
    tfm::reportfln("compressor_type:\t%s", nameof::nameof_type<t_Compressor>());
    tfm::reportfln("is_in_online:\t%d", int(t_TextLoader::IN_ONLINE));
    tfm::reportfln("enable_progrep:\t%d", int(enable_progrep));

    Timer        t;
    t_Compressor cmpr;

    double    compression_sec = 0.0;
    size_type milestone       = enable_progrep ? 10000 : MAX_SIZE_INT;

    {
        t_TextLoader text_loader(input_path);

        t.start();
        while (!text_loader.eof()) {
            cmpr.extend(text_loader.get());
            if (cmpr.get_num_chars() == milestone) {
                t.stop();
                compression_sec += t.get_elapsed_sec();
                tfm::printfln("[Progress_Report]");
                tfm::reportfln("num_process_chars:\t%d", cmpr.get_num_chars());
                tfm::reportfln("compression_sec:\t%g", compression_sec);
                milestone *= 10;
                t.start();
            }
        }
        t.stop();
    }

    compression_sec += t.get_elapsed_sec();
    tfm::printfln("[Final_Report]");
    tfm::reportfln("compression_sec:\t%g", compression_sec);

    const size_type num_runs  = utils::compute_num_runs(cmpr);
    const size_type num_chars = cmpr.get_num_chars();
    tfm::reportfln("num_runs:\t%d", num_runs);
    tfm::reportfln("num_chars:\t%d", num_chars);
    tfm::reportfln("compression_ratio:\t%g", num_runs / double(num_chars));

    const size_type mem_bytes = cmpr.get_memory_in_bytes();
    tfm::reportfln("memory_in_bytes:\t%d", mem_bytes);
    tfm::reportfln("memory_in_MiB:\t%g", mem_bytes / (1024.0 * 1024.0));

    cmpr.show_memory_statistics();
    cmpr.show_detailed_statistics();
    cmpr.show_monitored_statistics();

    if (enable_test) {
        tfm::printfln("Now teesting the data structure...");
        cmpr.test_all();
        tfm::printfln("No Problem!");
    }

    if (enable_check) {
        tfm::printfln("Now checking the decoded text...");

        std::vector<uchar_type> orig_text = utils::load_text(input_path);
        std::vector<uchar_type> dec_text  = utils::decode(cmpr);

        std::reverse(orig_text.begin(), orig_text.end());
        if (utils::is_equivalent_text(orig_text, dec_text)) {
            tfm::printfln("No Problem!");
        } else {
            tfm::errorfln("What's up!?");
        }
    }

    return 0;
}

template <class t_Compressor>
int perf_rlbwt(const cmd_line_parser::parser& p, const bool is_in_online) {
    if (is_in_online) {
        return perf_rlbwt<t_Compressor, TextLoader<true>>(p);
    } else {
        return perf_rlbwt<t_Compressor, TextLoader<false>>(p);
    }
}

int main(int argc, char** argv) {
#ifndef NDEBUG
    tfm::warnfln("The code is running in debug mode.");
#endif

    auto p = make_parser(argc, argv);
    if (!p.parse()) {
        return 1;
    }

    const auto compressor_str = p.get<std::string>("compressor", "glfig_serialized_16");
    const auto is_in_online   = p.get<bool>("is_in_online", false);

    // LFIG_Naive
    if (compressor_str == "lfig_naive") {
        return perf_rlbwt<compressors::lfig_naive<DIV_BOUND>::compressor_type>(p, is_in_online);
    }
    // Grouped LFIG_Serialized with hints
    else if (compressor_str == "glfig_serialized_8") {
        return perf_rlbwt<compressors::glfig_serialized<8, DIV_BOUND>::compressor_type>(p, is_in_online);
    } else if (compressor_str == "glfig_serialized_16") {
        return perf_rlbwt<compressors::glfig_serialized<16, DIV_BOUND>::compressor_type>(p, is_in_online);
    } else if (compressor_str == "glfig_serialized_32") {
        return perf_rlbwt<compressors::glfig_serialized<32, DIV_BOUND>::compressor_type>(p, is_in_online);
    } else if (compressor_str == "glfig_serialized_64") {
        return perf_rlbwt<compressors::glfig_serialized<64, DIV_BOUND>::compressor_type>(p, is_in_online);
    }

    p.help();
    return 1;
}