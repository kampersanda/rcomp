/**
 * @file basics.hpp
 */
#pragma once

#include <stdint.h>
#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <string_view>
#include <vector>

#include "nameof/nameof.hpp"
#include "tinyformat/tinyformat.h"

#include "abort_if.hpp"

namespace rcomp {

using size_type = uint64_t;
using offset_type = int64_t;
using uchar_type = uint8_t;
using loint_type = uint64_t;  // for list-order labeling
using text_type = std::vector<uchar_type>;

static constexpr size_type MAX_SIZE_INT = std::numeric_limits<size_type>::max();
static constexpr offset_type MAX_OFFSET_INT = std::numeric_limits<offset_type>::max();
static constexpr uchar_type END_MARKER = '\0';

constexpr char to_print(uchar_type c) {
    return c == END_MARKER ? '$' : c;
}

inline std::ostream& operator<<(std::ostream& out, const text_type& text) {
    for (size_type i = 0; i < text.size(); i++) {
        out << to_print(text[i]);
    }
    return out;
}

struct run_type {
    uchar_type chr;
    size_type exp;
};
inline run_type make_run(uchar_type chr) {
    return {chr, size_type(1)};
}
inline run_type make_run(uchar_type chr, size_type exp) {
    return {chr, exp};
}
inline std::ostream& operator<<(std::ostream& out, const run_type& rn) {
    tfm::format(out, "(%c,%d)", to_print(rn.chr), rn.exp);
    return out;
}

template <class T>
struct range_type {
    T beg;
    T end;
};
inline range_type<const uchar_type*> make_range(std::string_view vec) {
    auto data = reinterpret_cast<const uchar_type*>(vec.data());
    return {data, data + vec.size()};
}
template <class T>
inline range_type<const T*> make_range(const std::vector<T>& vec) {
    return {vec.data(), vec.data() + vec.size()};
}

// https://ny23.hatenadiary.org/entry/20111129/p1
template <class T, class Compare>
[[maybe_unused]] static size_type get_bst_memory_in_bytes(const std::set<T, Compare>& bst) {
    return bst.size() * (sizeof(int) +  // Node color (enum; int)
                         3 * sizeof(void*) +  // Pointers to parent, left, and right
                         sizeof(T)  // Key
                        );
}
template <class T, class V, class Compare>
[[maybe_unused]] static size_type get_bst_memory_in_bytes(const std::map<T, V, Compare>& bst) {
    return bst.size() * (sizeof(int) +  // Node color (enum; int)
                         3 * sizeof(void*) +  // Pointers to parent, left, and right
                         sizeof(std::pair<T, V>)  // Key, Value
                        );
}

namespace lo_common {

static constexpr loint_type LO_MIN_VALUE = 0;
static constexpr loint_type LO_MAX_VALUE = std::numeric_limits<loint_type>::max();

template <class T>
inline loint_type get_basic_order(const T* x) {
    return x->get_data().get_order();
}

template <class T>
inline loint_type get_star_order(const T* x) {
    const loint_type xv = x->get_data().get_order();
    return xv == LO_MIN_VALUE ? LO_MAX_VALUE : xv;
}

template <class T>
inline loint_type get_distance(const T* x, const T* y) {
    return get_star_order(y) - get_basic_order(x);
}

template <class T>
inline bool compare_grouped_order(const T* x, const T* y) {
    // First, compare with the L-node order
    const loint_type x_order = lo_common::get_basic_order(x);
    const loint_type y_order = lo_common::get_basic_order(y);
    if (x_order != y_order) {
        return x_order < y_order;
    }

    // Second, compare with the offset position in the L-node
    const size_type x_second_order = x->get_data().get_second_order();
    const size_type y_second_order = y->get_data().get_second_order();
    return x_second_order < y_second_order;
}

}  // namespace lo_common

}  // namespace rcomp