#include <string>

#include "rlbwt_types.hpp"
#include "utils.hpp"

using namespace rcomp;

int main(int argc, char** argv) {
    // Input text
    const std::string text = "abaababaab";

    // Construct the r-index by appending characters (with end-marker $) in reverse.
    // Note that '\0' is used for the end marker (i.e., the text should not contain '\0').
    rlbwt_types::glfig_serialized<8>::type rlbwt;
    for (size_t i = 1; i <= text.size(); i++) {
        rlbwt.extend(text[text.size() - i]);
    }

    // Extract the resulted BWT-runs
    std::cout << "BWT-runs: ";
    rlbwt.output_runs([](const run_type& r) { std::cout << r << ","; });
    std::cout << std::endl;

    // Decode the original text (except $) from the index
    std::string decoded;
    rlbwt.decode_text([&](uchar_type c) { decoded.push_back(c); });
    std::reverse(decoded.begin(), decoded.end());  // need to be reversed
    std::cout << "text:    " << text << std::endl;
    std::cout << "decoded: " << decoded << std::endl;

    return 0;
}