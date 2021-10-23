/**
 * @file StaticVector.hpp
 */
#pragma once

#include <array>

#include "abort_if.hpp"
#include "basics.hpp"

namespace rcomp {

/**
 * Vector with static capacity.
 *
 * @tparam T data type
 * @tparam N capacity
 */
template <class T, size_type N>
class StaticVector {
  private:
    std::array<T, N> m_vec;
    size_type m_size = 0;

  public:
    //! Default constructor
    StaticVector() = default;

    //! Get if the data structure is empty.
    inline bool is_empty() const {
        return m_size == 0;
    }
    inline bool is_full() const {
        return m_size == N;
    }

    //! Get the number of used elements.
    inline size_type get_size() const {
        return m_size;
    }

    inline void clear() {
        m_size = 0;
    }

    inline void push_back(T v) {
        DEBUG_ABORT_IF(is_full());
        m_vec[m_size++] = v;
    }

    inline T operator[](size_type i) const {
        DEBUG_ABORT_IF_LE(m_size, i);
        return m_vec[i];
    }

    inline bool is_member(T v) const {
        for (size_type i = 0; i < m_size; i++) {
            if (m_vec[i] == v) {
                return true;
            }
        }
        return false;
    }

    inline T front() const {
        DEBUG_ABORT_IF(is_empty());
        return m_vec[0];
    }
    inline T back() const {
        DEBUG_ABORT_IF(is_empty());
        return m_vec[m_size - 1];
    }

    inline T* begin() {
        return m_vec.data();
    }
    inline const T* begin() const {
        return m_vec.data();
    }
    inline T* end() {
        return m_vec.data() + m_size;
    }
    inline const T* end() const {
        return m_vec.data() + m_size;
    }
};

}  // namespace rcomp