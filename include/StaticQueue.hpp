/**
 * @file StaticQueue.hpp
 */
#pragma once

#include <array>

#include "abort_if.hpp"
#include "basics.hpp"

namespace rcomp {

/**
 * Ring-queue with static capacity.
 *
 * @tparam T data type
 * @tparam N capacity
 */
template <class T, size_type N>
class StaticQueue {
  private:
    std::array<T, N> m_queue;
    size_type m_bpos = 0;
    size_type m_size = 0;

  public:
    //! Default constructor
    StaticQueue() = default;

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

    inline void push(T v) {
        DEBUG_ABORT_IF(is_full());
        m_queue[get_pos(m_bpos + m_size)] = v;
        m_size += 1;
    }

    inline T pop() {
        DEBUG_ABORT_IF(is_empty());
        const T v = m_queue[m_bpos];
        m_bpos = get_pos(m_bpos + 1);
        m_size = m_size - 1;
        return v;
    }

    inline bool is_member(T v) const {
        for (size_type i = 0; i < m_size; i++) {
            if (m_queue[get_pos(m_bpos + i)] == v) {
                return true;
            }
        }
        return false;
    }

  private:
    inline static size_type get_pos(size_type i) {
        return i % N;
    }
};

}  // namespace rcomp