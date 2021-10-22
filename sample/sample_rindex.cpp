#include <string>

#include "rindex_types.hpp"
#include "utils.hpp"

using namespace rcomp;

int main(int argc, char** argv) {
    // Input text
    const std::string text = "abaababaab";

    // Construct the r-index by appending characters (with end-marker $) in reverse.
    // Note that '\0' is used for the end marker (i.e., the text should not contain '\0').
    rindex_types::glfig_naive<8>::type rindex;
    for (size_t i = 1; i <= text.size(); i++) {
        rindex.extend(text[text.size() - i]);
    }

    // Extract the resulted BWT-runs
    std::cout << "BWT-runs: ";
    rindex.output_runs([](const run_type& r) { std::cout << r << ","; });
    std::cout << std::endl;

    // Decode the original text (except $) from the index
    std::string decoded;
    rindex.decode_text([&](uchar_type c) { decoded.push_back(c); });
    std::reverse(decoded.begin(), decoded.end());  // need to be reversed
    std::cout << "text:    " << text << std::endl;
    std::cout << "decoded: " << decoded << std::endl;

    // Count the occurrences of the (reversed) query
    const std::string query = "aaba";  // i.e., "abaa" in the original order
    const size_type occ = rindex.count(make_range(query));
    std::cout << "count(" << query << ") = " << occ << std::endl;

    // Locate the (reversed) query
    std::cout << "locate(" << query << ") = {";
    rindex.locate(make_range(query), [&](size_type pos) {
        // We will get pos such that text_r[..pos] = "..aaba", where text_r = reverse(text+'$').
        // In other words, we can extract the original position starting at "abaa" as follows.
        std::cout << text.size() - pos << ",";
    });
    std::cout << "}" << std::endl;

    return 0;
}