/**
 * @file Successor.hpp
 */
#pragma once

#include <map>

#include "basics.hpp"

namespace rcomp {

/**
 * A successor class using std::map.
 *
 * @tparam T The key type.
 * @tparam V The value type.
 */
template <class T, class V>
class Successor {
  public:
    using this_type = Successor<T, V>;
    using key_type = T;
    using value_type = V;

  private:
    std::map<key_type, value_type> m_bst;

  public:
    //! Default constructor
    Successor() = default;

    //! Default destructor
    virtual ~Successor() = default;

    //! Copy constructor (deleted)
    Successor(const Successor&) = delete;

    //! Copy constructor (deleted)
    Successor& operator=(const Successor&) = delete;

    //! Move constructor
    Successor(Successor&&) noexcept = default;

    //! Move constructor
    Successor& operator=(Successor&&) noexcept = default;

    //! Get the allocated memory in bytes.
    size_type get_memory_in_bytes(bool include_this = true) const {
        const size_type this_bytes = sizeof(*this) * include_this;
        return this_bytes + get_bst_memory_in_bytes(m_bst);
    }

    //! Insert key k associated with value v.
    inline void insert(key_type k, value_type v) {
        auto itr = m_bst.find(k);
        if (itr == m_bst.end()) {
            m_bst.insert(std::make_pair(k, v));
        } else {
            itr->second = v;
        }
    }

    //! Remove the item with key k.
    inline void remove(key_type k) {
        m_bst.erase(k);
    }

    //! Return the value associated with key k.
    //! If such an item is not stored, return ng_value.
    inline value_type find(key_type k, value_type ng_value) const {
        auto itr = m_bst.find(k);
        if (itr == m_bst.end()) {
            return ng_value;
        } else {
            return itr->second;
        }
    }

    //! Search the smallest key k' such that k <= k'.
    //! Return the item if found or ng_key otherwise.
    inline std::tuple<key_type, value_type> search(key_type k, key_type ng_key) const {
        auto itr = m_bst.lower_bound(k);
        if (itr == m_bst.end()) {
            return {ng_key, value_type{}};
        }
        return {itr->first, itr->second};
    }

    //! Get the begin iterator.
    inline auto begin() const {
        return m_bst.begin();
    }

    //! Get the end iterator.
    inline auto end() const {
        return m_bst.end();
    }
};

}  // namespace rcomp