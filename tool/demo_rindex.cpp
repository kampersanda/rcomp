/**
 * @file demo_rindex.cpp
 */
#include "Timer.hpp"
#include "io.hpp"
#include "rindex_types.hpp"
#include "utils.hpp"

#include "cmd_line_parser/parser.hpp"

using namespace rcomp;

static constexpr size_type DIV_BOUND = 7;

cmd_line_parser::parser make_parser(int argc, char** argv) {
    cmd_line_parser::parser p(argc, argv);
    p.add("input_path", "Input file path of text");
    p.add("rindex_type", "Rindex data structure type: lfig | glfig_[8|16|32|64] (default=glfig_16)", "-t", false);
    p.add("step", "Number of characters to index in a single step (default=1000000)", "-s", false);
    return p;
}

template <class t_Rindex>
int demo_rindex(const cmd_line_parser::parser& p) {
    const auto input_path = p.get<std::string>("input_path");
    const auto step = p.get<size_type>("step", 1000000);

    auto ifs = io::make_ifstream(input_path);
    size_type milestone = step;

    t_Rindex rindex;
    Timer timer;

    auto interactive = [&]() {
        const size_type mem_byte = rindex.get_memory_in_bytes();
        tfm::reportfln("%d characters indexed in %d bytes = %g KiB = %g MiB.",  //
                       rindex.get_num_chars(), mem_byte, utils::to_KiB(mem_byte), utils::to_MiB(mem_byte));

        while (true) {
            tfm::printfln("1. Enter query string to search.");
            tfm::printfln("2. Enter \"exit\" to continue indexing.");

            tfm::printf("> ");

            std::string str;
            std::cin >> str;

            if (str == "exit") {
                break;
            }

            const auto query = make_range(str);

            timer.start();
            const size_type occ = rindex.count(query);
            timer.stop();

            tfm::reportfln("Count(\"%s\") = %d, done in %g micro sec.", str, occ, timer.get_elapsed_us());

            if (occ == 0) {
                continue;
            }

            tfm::printfln("1. Enter '1' to run locate with print.");
            tfm::printfln("2. Enter '2' to run locate without print.");
            tfm::printfln("3. Enter another not to run locate.");

            tfm::printf("> ");

            std::string mode;
            std::cin >> mode;

            timer.start();
            if (mode == "1") {
                tfm::reportf("Locate(\"%s\") = {", str);
                rindex.locate(query, [&](size_type pos) { tfm::reportf("%d, ", pos - str.size()); });
                tfm::reportfln("}");
            } else if (mode == "2") {
                rindex.locate(query, [&](size_type pos) {});
            } else {
                continue;
            }
            timer.stop();

            tfm::reportfln("Locate query, done in %g micro sec, %g micro sec per occ.",  //
                           timer.get_elapsed_us(), timer.get_elapsed_us() / occ);
        }
    };

    tfm::printfln("Constructing r-index...");

    while (ifs.peek() != std::ios::traits_type::eof()) {
        char c;
        ifs.get(c);
        rindex.extend(c);
        if (rindex.get_num_chars() == milestone) {
            interactive();
            milestone += step;
        }
    }

    interactive();
    tfm::printfln("Thanks!");

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

    const auto rindex_type_str = p.get<std::string>("rindex_type", "glfig_16");

    if (rindex_type_str == "lfig") {
        return demo_rindex<rindex_types::lfig_naive<DIV_BOUND>::type>(p);
    } else if (rindex_type_str == "glfig_8") {
        return demo_rindex<rindex_types::glfig_serialized<8, DIV_BOUND>::type>(p);
    } else if (rindex_type_str == "glfig_16") {
        return demo_rindex<rindex_types::glfig_serialized<16, DIV_BOUND>::type>(p);
    } else if (rindex_type_str == "glfig_32") {
        return demo_rindex<rindex_types::glfig_serialized<32, DIV_BOUND>::type>(p);
    } else if (rindex_type_str == "glfig_64") {
        return demo_rindex<rindex_types::glfig_serialized<64, DIV_BOUND>::type>(p);
    }

    p.help();
    return 1;
}
