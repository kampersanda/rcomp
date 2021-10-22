#pragma once

#ifdef __APPLE__
#include <mach/mach.h>
#endif

#include <unistd.h>
#include <fstream>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <string_view>

#include "abort_if.hpp"
#include "basics.hpp"

namespace rcomp::io {

//! Make the input stream.
[[maybe_unused]] static std::ifstream make_ifstream(const std::string& filepath) {
    std::ifstream ifs(filepath);
    ABORT_IF(!ifs);
    return ifs;
}

//! Make the output stream.
[[maybe_unused]] static std::ofstream make_ofstream(const std::string& filepath) {
    std::ofstream ofs(filepath);
    ABORT_IF(!ofs);
    return ofs;
}

//! Get the file size in bytes.
[[maybe_unused]] static size_type get_filesize(const std::string& filepath) {
    std::ifstream ifs(filepath, std::ios::binary | std::ios::ate);
    ABORT_IF(!ifs.good());
    const auto bytes = static_cast<size_type>(ifs.tellg());
    ifs.close();
    return bytes;
}

//! Load the text in forward or backward order.
[[maybe_unused]] static text_type load_text(const std::string& filepath, bool is_reverse = false) {
    const size_type filesize = get_filesize(filepath);
    text_type text(filesize);
    {
        auto ifs = make_ifstream(filepath);
        ifs.read(reinterpret_cast<char*>(text.data()), filesize);
    }
    if (is_reverse) {
        std::reverse(text.begin(), text.end());
    }
    return text;
}

template <typename T>
[[maybe_unused]] static void load_pod(std::istream& is, T& val) {
    static_assert(std::is_pod<T>::value);
    is.read(reinterpret_cast<char*>(&val), sizeof(T));
}

template <typename T, typename Allocator>
[[maybe_unused]] static void load_vec(std::istream& is, std::vector<T, Allocator>& vec) {
    size_t n;
    load_pod(is, n);
    vec.resize(n);
    is.read(reinterpret_cast<char*>(vec.data()), static_cast<std::streamsize>(sizeof(T) * n));
}

template <typename T>
[[maybe_unused]] static void save_pod(std::ostream& os, T const& val) {
    static_assert(std::is_pod<T>::value);
    os.write(reinterpret_cast<char const*>(&val), sizeof(T));
}

template <typename T, typename Allocator>
[[maybe_unused]] static void save_vec(std::ostream& os, std::vector<T, Allocator> const& vec) {
    static_assert(std::is_pod<T>::value);
    size_t n = vec.size();
    save_pod(os, n);
    os.write(reinterpret_cast<char const*>(vec.data()), static_cast<std::streamsize>(sizeof(T) * n));
}

//! Get the (maximum) resident set size to compute the real memory usage (in this process).
//! c.f. http://www.tkl.iis.u-tokyo.ac.jp/~ynaga/cedar/bench.cc
[[maybe_unused]] static size_t get_rss_in_bytes() {
#ifdef __APPLE__
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
    task_info(current_task(), TASK_BASIC_INFO, reinterpret_cast<task_info_t>(&t_info), &t_info_count);
    return t_info.resident_size;
#else
    FILE* fp = std::fopen("/proc/self/statm", "r");
    size_t dummy(0), vm(0);
    ABORT_IF_EQ(std::fscanf(fp, "%ld %ld ", &dummy, &vm), EOF);  // get resident (see procfs)
    std::fclose(fp);
    return vm * ::getpagesize();
#endif
}

}  // namespace rcomp::io