/**
 * @file intrinsics.hpp
 * @author Shunsuke Kanda (kampersanda)
 */
#pragma once

#include <xmmintrin.h>

#if defined(__GNUC__) || defined(__clang__)
#define __INTRIN_INLINE inline __attribute__((__always_inline__))
#elif defined(_MSC_VER)
#define __INTRIN_INLINE inline __forceinline
#else
#define __INTRIN_INLINE inline
#endif

namespace rcomp::intrinsics {

template <typename T>
__INTRIN_INLINE void prefetch(const T* p) {
    _mm_prefetch(reinterpret_cast<const char*>(p), _MM_HINT_T0);
}

}  // namespace rcomp::intrinsics