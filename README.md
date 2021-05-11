# R-comp: Online RLBWT compression in optimal-time and BWT-runs bounded space

[![experimental](http://badges.github.io/stability-badges/dist/experimental.svg)](http://github.com/badges/stability-badges)

This library provides a C++17 implementation of **R-comp**, an online RLBWT compression algorithm in optimal-time and BWT-runs bounded space.

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

After the commands, the executables will be produced in `build/bin` directory.

The code is written in C++17, so please install g++ >= 7.0 or clang >= 4.0. For the build system, CMake >= 3.0 have to be installed to compile the library.

The library employs the third-party libraries [cmd\_line\_parser](https://github.com/jermp/cmd_line_parser), [doctest](https://github.com/onqtam/doctest), [nameof](https://github.com/Neargye/nameof) and [tinyformat](https://github.com/c42f/tinyformat), whose header files are contained in this repository.

The code has been tested only on Mac OS X and Linux. That is, this library considers only UNIX-compatible OS.

## Sample usage

The executable `perf_rlbwt` provides a benchmark of RLBWT construction, as follows.

```shell
$ echo -n GATCAATGAGGTGGACACCAGAGGCGGGGACTTGT > sample.txt
$ ./bin/perf_rlbwt sample.txt
[Basic]
input_path:	sample.txt
compressor_type:	rcomp::Compressor_GroupLFIG<rcomp::GroupLIndex<rcomp::GroupLData_Serialized<16, true, true, 2> >, rcomp::GroupFIndex<rcomp::GroupFData_Serialized<rcomp::GroupLData_Serialized<16, true, true, 2> > >, 7>
is_in_online:	0
enable_progrep:	1
[Final_Report]
compression_sec:	0
num_runs:	27
num_chars:	36
compression_ratio:	0.75
memory_in_bytes:	1763
memory_in_MiB:	0.00168133
[Memory_Compressor]
m_lindex:	648
m_findex:	803
[Memory_LIndex]
m_list:	648
[Memory_LList]
size_info:	9
link_info:	75
unit_info:	336
sumexp_hint:	24
lookup_hint:	54
[Memory_FIndex]
m_bsts:	180
m_list:	623
[Memory_FList]
size_info:	15
link_info:	50
unit_info:	348
sumexp_hint:	40
lookup_hint:	90
[Detail_LIndex]
num_hlinks:	3
num_lnodes:	3
ratio_hlinks:	1
size_units:	25
capa_units:	26
ave_size_units:	8.33333
ave_capa_units:	8.66667
[Detail_FIndex]
num_tlinks:	3
num_fnodes:	5
ratio_tlinks:	0.6
size_units:	26
capa_units:	28
ave_size_units:	5.2
ave_capa_units:	5.6
Now teesting the data structure...
No Problem!
Now checking the decoded text...
No Problem!
```

## Test

The unit tests are written using [doctest](https://github.com/onqtam/doctest). After compiling, you can run tests with the following command:

```shell
$ make test
```

## Authors

- [Takaaki Nishimoto](https://github.com/TNishimoto)
- [Shunsuke Kanda](https://github.com/kampersanda) (Creator)
- [Yasuo Tabei](https://github.com/tb-yasu)

## Licensing

This program is available for only academic use, basically. For the academic use, please keep [MIT License](https://github.com/kampersanda/rcomp/blob/main/LICENSE). For the commercial use, please keep GPL 2.0 and make a contact to one of the authors.

## Related software

- [renum](https://github.com/TNishimoto/renum) is a C++ implementation of enumeration of characteristic substrings in BWT-runs bounded space.
- [rlbwt\_iterator](https://github.com/TNishimoto/rlbwt_iterator) is a C++ implementation of some iterators in BWT-runs bounded space.