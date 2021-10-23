/**
 * @file FData_Naive.hpp
 */
#pragma once

#include "LinkedList.hpp"

namespace rcomp {

/**
 * Naive F-node data type.
 *
 * @tparam t_LData The data type of L-node.
 */
template <class t_LData>
class FData_Naive {
  public:
    using lnode_type = typename LinkedList<t_LData>::node_type;

  private:
    lnode_type* m_lnode = nullptr;
    lnode_type* m_tnode = nullptr;
    uint8_t m_weight = 0;

  public:
    //! Default constructor
    FData_Naive() = default;

    //! Default destructor
    virtual ~FData_Naive() = default;

    //! Copy constructor (deleted)
    FData_Naive(const FData_Naive&) = delete;

    //! Copy constructor (deleted)
    FData_Naive& operator=(const FData_Naive&) = delete;

    //! Move constructor
    FData_Naive(FData_Naive&&) noexcept = default;

    //! Move constructor
    FData_Naive& operator=(FData_Naive&&) noexcept = default;

    //! Get the allocated memory in bytes.
    inline size_type get_memory_in_bytes(bool include_this = true) const {
        return sizeof(*this) * include_this;
    }

    inline void reset_tlink() {
        m_tnode = nullptr;
    }

    inline void set_lnode(lnode_type* v) {
        m_lnode = v;
    }
    inline void set_tnode(lnode_type* v) {
        m_tnode = v;
    }

    inline uchar_type get_chr() const {
        return m_lnode->get_data().get_chr();
    }
    inline size_type get_exp() const {
        return m_lnode->get_data().get_exp();
    }
    inline loint_type get_order() const {
        return m_lnode->get_data().get_order();
    }
    inline lnode_type* get_lnode() const {
        return m_lnode;
    }
    inline lnode_type* get_tnode() const {
        return m_tnode;
    }
    inline offset_type get_tofst() const {
        return m_tnode->get_data().get_hofst();
    }

    inline size_type get_weight() const {
        return static_cast<size_type>(m_weight);
    }
    inline void set_weight(size_type w) {
        DEBUG_ABORT_IF_LT(UINT8_MAX, w);
        m_weight = static_cast<uint8_t>(w);
    }
    inline void add_weight(size_type w) {
        set_weight(m_weight + w);
    }
    inline void sub_weight(size_type w) {
        set_weight(m_weight - w);
    }

    inline size_type get_sae() const {
        return m_lnode->get_data().get_sae();
    }
};

}  // namespace rcomp
