/**
 * @file basics.hpp
 * @author Shunsuke Kanda (kampersanda)
 */
#pragma once

#include <stdint.h>
#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

#include "nameof/nameof.hpp"
#include "tinyformat/tinyformat.h"

#include "abort_if.hpp"

namespace rcomp {

using size_type   = uint64_t;
using offset_type = int64_t;
using uchar_type  = uint8_t;
using loint_type  = uint64_t;  // for list labeling

static constexpr uchar_type END_MARKER = '\0';

inline char to_print(uchar_type c) {
    return c == END_MARKER ? '$' : c;
}

static constexpr size_type   MAX_SIZE_INT   = std::numeric_limits<size_type>::max();
static constexpr offset_type MAX_OFFSET_INT = std::numeric_limits<offset_type>::max();

struct run_type {
    uchar_type chr;
    size_type  exp;
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

}  // namespace lo_common

}  // namespace rcomp