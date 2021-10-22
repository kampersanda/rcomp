/**
 * @file transform.cpp
 */
#include "TextLoader.hpp"
#include "Timer.hpp"
#include "io.hpp"
#include "rlbwt_types.hpp"
#include "utils.hpp"

#include "cmd_line_parser/parser.hpp"

using namespace rcomp;

static constexpr size_type DIV_BOUND = 7;

cmd_line_parser::parser make_parser(int argc, char** argv) {
    cmd_line_parser::parser p(argc, argv);
    p.add("input_path", "Input file path of text");
    p.add("output_path", "Output file path of BWT-text");
    p.add("rlbwt_type", "Rlbwt data structure type: lfig | glfig_[8|16|32|64] (default=glfig_16)", "-t", false);
    p.add("reverse_mode", "Input the file in reverse order? (default=1)", "-r", false);
    return p;
}

template <class t_Rlbwt>
int transform(const cmd_line_parser::parser& p) {
    const auto input_path = p.get<std::string>("input_path");
    const auto output_path = p.get<std::string>("output_path");
    const auto reverse_mode = p.get<bool>("reverse_mode", true);

    if (!reverse_mode) {
        tfm::warnfln(
            "Since the text will be input in the original order, "
            "the BWT for the reversed text will be constructed.");
    }

    tfm::printfln("Constructing now...");
    t_Rlbwt rlbwt;
    {
        TextLoader loader(input_path, reverse_mode);
        while (!loader.eof()) {
            rlbwt.extend(loader.next());
            if (rlbwt.get_num_chars() % 10'000'000 == 0) {
                tfm::printfln("%d...", rlbwt.get_num_chars());
            }
        }
    }
    tfm::printfln("RLBWT was constructed for %d chars.", rlbwt.get_num_chars());

    tfm::printfln("Outputting now...");
    size_type num_runs = 0;
    {
        auto ofs = io::make_ofstream(output_path);
        rlbwt.output_runs([&](const run_type& rn) {
            num_runs += 1;
            for (size_type i = 0; i < rn.exp; i++) {
                ofs << rn.chr;
            }
        });
    }
    tfm::printfln("BWT-text was output.");
    tfm::printfln("The number of resulting runs was %d.", num_runs);

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
        return transform<rlbwt_types::lfig_naive<DIV_BOUND>::type>(p);
    } else if (rlbwt_type_str == "glfig_8") {
        return transform<rlbwt_types::glfig_serialized<8, DIV_BOUND>::type>(p);
    } else if (rlbwt_type_str == "glfig_16") {
        return transform<rlbwt_types::glfig_serialized<16, DIV_BOUND>::type>(p);
    } else if (rlbwt_type_str == "glfig_32") {
        return transform<rlbwt_types::glfig_serialized<32, DIV_BOUND>::type>(p);
    } else if (rlbwt_type_str == "glfig_64") {
        return transform<rlbwt_types::glfig_serialized<64, DIV_BOUND>::type>(p);
    }

    p.help();
    return 1;
}
