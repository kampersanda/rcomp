# R-comp: Online RLBWT compression in optimal-time and BWT-runs bounded space

This is an experimental library of R-comp, an online RLBWT compression algorithm in optimal-time and BWT-runs bounded space.

## Build instructions

You can download and compile the library with the following commands:

```shell
$ git clone https://github.com/kampersanda/rcomp.git
$ cd rcomp
$ mkdir build
$ cd build
$ cmake ..
$ make -j
```

The code is written in C++17, so please install g++ >= 7.0 or clang >= 4.0. For the build system, CMake >= 3.0 have to be installed to compile the library.

The library employs the third-party libraries [cmd\_line\_parser](https://github.com/jermp/cmd_line_parser), [doctest](https://github.com/onqtam/doctest), [nameof](https://github.com/Neargye/nameof) and [tinyformat](https://github.com/c42f/tinyformat), whose header files are contained in this repository.

The code has been tested only on Mac OS X and Linux. That is, this library considers only UNIX-compatible OS.

## Implementations

The library implements several data structures and provides the following variants of R-comp defined in `rlbwt_types`:

- `rlbwt_types::lfig_naive` is a straightforward implementation of the LF-interval graph with `O(r)` nodes, and
- `rlbwt_types::glfig_serialized<g>` is a spece-efficient implementation of the LF-interval graph with `O(r/g)` nodes,

where `r` is the number of BWT-runs.

Also, the library implements r-index on these data structures, providing `count` and `locate` queries in the compressed space. In the same manner as `rlbwt_types`, the variants are defined in `rindex_types`.

## Limitations

- An input text must NOT contain the `0x00` character because it is used as a special end marker.
- In the current version, class `GroupedFIndex` resorts to static global variables. Please do NOT create multiple instances of `glfig_serialized` in a single process.

## Sample usage

### RLBWT

`sample/sample_rlbwt.cpp` provides a sample usage.

```c++
#include <string>

#include "rlbwt_types.hpp"
#include "utils.hpp"

using namespace rcomp;

int main(int argc, char** argv) {
    // Input text
    const std::string text = "abaababaab";

    // Construct the RLBWT by appending characters (with end-marker $) in reverse.
    // Note that '\0' is used for the end marker (i.e., the text should not contain '\0').
    rlbwt_types::glfig_serialized<8>::type rlbwt;
    for (size_t i = 1; i <= text.size(); i++) {
        rlbwt.extend(text[text.size() - i]);
    }

    // Extract the resulted BWT-runs
    std::cout << "BWT-runs: ";
    rlbwt.output_runs([](const run_type& r) { std::cout << r << ","; });
    std::cout << std::endl;

    // Decode the original text (except $) from the RLBWT
    std::string decoded;
    rlbwt.decode_text([&](uchar_type c) { decoded.push_back(c); });
    std::reverse(decoded.begin(), decoded.end());  // need to be reversed
    std::cout << "text:    " << text << std::endl;
    std::cout << "decoded: " << decoded << std::endl;

    return 0;
}
```

The output will be

```
BWT-runs: (b,3),(a,1),(b,1),($,1),(a,5),
text:    abaababaab
decoded: abaababaab
```

### r-index

`sample/sample_rindex.cpp` provides a sample usage.

```c++
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
```

The output will be

```
BWT-runs: (b,3),(a,1),(b,1),($,1),(a,5),
text:    abaababaab
decoded: abaababaab
count(aaba) = 2
locate(aaba) = {5,0,}
```

## Performance test

### RLBWT

The executable `perf/perf_rlbwt` measures the performance of R-comp. The command line options are printed by specifying the parameter `-h`.

```
$ ./perf/perf_rlbwt -h
Usage: ./perf/perf_rlbwt [-h,--help] input_path [-t rlbwt_type] [-r reverse_mode] [-T enable_test]

 input_path
    Input file path of text
 [-t rlbwt_type]
    Rlbwt data structure type: lfig | glfig_[8|16|32|64] (default=glfig_16)
 [-r reverse_mode]
    Load the text in reverse? (default=1)
 [-T enable_test]
    Test the data structure? (default=1)
 [-h,--help]
    Print this help text and silently exits.
```

- For data structure type `t`, `lfig` is a straight forward implementation, and `glfig_g` is a memory-efficient implementation by the grouping technique with group size `g`. 
- When `r` is set to `1`, the text will be input in reverse order to build the RLBWT for the text in the original order.
- When `T` is set to `1`, the correctness of the result will be tested.

For example, for dataset  `alice29.txt`, the following command measures the performace of R-comp with data structure `glfig_16` and shows the detailed statistics.

```
$ ./perf/perf_rlbwt alice29.txt
[Input_Params]
input_path:     alice29.txt
rlbwt_type:     rcomp::Rlbwt_GLFIG<rcomp::GroupedLFIntervalGraph<rcomp::GroupedLData_Serialized<16, false, 2, true, true>, rcomp::GroupedFData_Serialized<rcomp::GroupedLData_Serialized<16, false, 2, true, true> >, 7> >
reverse_mode:   1
[Progress_Report]
num_chars:      10000
construction_sec:       0.002
[Progress_Report]
num_chars:      100000
construction_sec:       0.03
[Final_Report]
construction_sec:       0.048
num_runs:       66903
num_chars:      152090
compression_ratio:      0.439891
alloc_memory_in_bytes:  2892218
alloc_memory_in_MiB:    2.75823
peak_memory_in_bytes:   3604480
peak_memory_in_MiB:     3.4375
Testing the data structure now...
No Problem!
Testing the decoded text now...
No Problem!
```

### r-index

The executable `perf/perf_rindex` measures the performance of r-index on R-comp. The command line options are printed by specifying the parameter `-h`.

```
$ ./perf/perf_rindex -h
Usage: ./perf/perf_rindex [-h,--help] input_path [-t rindex_type] [-r reverse_mode] [-T enable_test]

 input_path
        Input file path of text
 [-t rindex_type]
        Rindex data structure type: lfig | glfig_[8|16|32|64] (default=glfig_16)
 [-r reverse_mode]
        Loading the text in reverse? (default=1)
 [-T enable_test]
        Testing the data structure? (default=1)
 [-h,--help]
        Print this help text and silently exits.
```

The parameter settings are the same as `perf_rlbwt`, and the following command measures the performace of r-index on R-comp with data structure `glfig_16`.

```
$ ./perf/perf_rindex alice29.txt 
[Input_Params]
input_path:     alice29.txt
rindex_type:    rcomp::Rindex_GLFIG<rcomp::GroupedLFIntervalGraph<rcomp::GroupedLData_Serialized<16, true, 2, true, true>, rcomp::GroupedFData_Serialized<rcomp::GroupedLData_Serialized<16, true, 2, true, true> >, 7> >
reverse_mode:   1
[Progress_Report]
num_chars:      10000
construction_sec:       0.005
[Progress_Report]
num_chars:      100000
construction_sec:       0.07
[Final_Report]
construction_sec:       0.116
num_runs:       66903
num_chars:      152090
compression_ratio:      0.439891
alloc_memory_in_bytes:  6408190
alloc_memory_in_MiB:    6.11133
peak_memory_in_bytes:   8417280
peak_memory_in_MiB:     8.02734
[Search_Settings]
num_trials:     10
num_queries:    1000
query_length:   8
query_seed:     13
Warming up now...
dummy:  18858
[Count_Query]
occ_per_query:  18.858
microsec_per_query:     4.19685
[Locate_Query]
occ_per_query:  37.716
microsec_per_query:     7.14312
microsec_per_occ:       0.189392
Testing the data structure now...
No Problem!
Testing the decoded text now...
No Problem!
```

## BWT tool

The executable `tool/transform` constructs the BWT text from a given text. You need to set `r` to `1` to output the BWT for the text in the original order.

```
$ ./tool/transform alice29.txt alice29.bwt -r 1
Constructing now...
RLBWT was constructed for 152090 chars.
Outputting now...
BWT-text was output.
The number of resulting runs was 66903.
```

Note that, when `r` is set to `0`, the BWT for the reversed text will be constructed.

```
$ ./tool/transform alice29.txt alice29.bwt -r 0
WARNING: Since the text will be input in the original order, the BWT for the reversed text will be constructed.
Constructing now...
RLBWT was constructed for 152090 chars.
Outputting now...
BWT-text was output.
The number of resulting runs was 66186.
```

## r-index demo

The executable `tool/demo_rindex` offers a demo of `count`/`locate` queries using r-index.

```
$ ./tool/demo_rindex alice29.txt
Constructing r-index...
152090 characters indexed in 15465584 bytes = 15103.1 KiB = 14.7491 MiB.
1. Enter query string to search.
2. Enter "exit" to continue indexing.
> Dinah
Count("Dinah") = 14, done in 11.3 micro sec.
1. Enter '1' to run locate with print.
2. Enter '2' to run locate without print.
3. Enter another not to run locate.
> 1
Locate("Dinah") = {43681, 4612, 5237, 33563, 35845, 5189, 33713, 32627, 32751, 21324, 36155, 4532, 4475, 32894, }
Locate query, done in 67.9 micro sec, 4.85 micro sec per occ.
1. Enter query string to search.
2. Enter "exit" to continue indexing.
> exit
Thanks!
```

## Unit test

The unit tests are written using [doctest](https://github.com/onqtam/doctest). After compiling, you can run tests with the following command.

```
$ make test
```

## Authors

- [Takaaki Nishimoto](https://github.com/TNishimoto)
- [Shunsuke Kanda](https://github.com/kampersanda) (Creator)
- [Yasuo Tabei](https://github.com/tb-yasu)

## Licensing

This program is available for only academic use, basically. For the academic use, please keep [MIT License](https://github.com/kampersanda/rcomp/blob/main/LICENSE). For the commercial use, please keep GPL 2.0 and make a contact to one of the authors.

If you use the library, please cite the following paper:

```
@inproceedings{
    ...
}
```

## Related software

- [renum](https://github.com/TNishimoto/renum) is a C++ implementation of enumeration of characteristic substrings in BWT-runs bounded space.
- [rlbwt\_iterator](https://github.com/TNishimoto/rlbwt_iterator) is a C++ implementation of some iterators in BWT-runs bounded space.