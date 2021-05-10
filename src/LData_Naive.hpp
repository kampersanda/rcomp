/**
 * @file LData_Naive.hpp
 * @author Shunsuke Kanda (kampersanda)
 */
#pragma once

#include "FData_Naive.hpp"
#include "LinkedList.hpp"

namespace rcomp {

/**
 * A class for naive data type of L-node.
 */
class LData_Naive {
  public:
    using fnode_type = typename LinkedList<FData_Naive<LData_Naive>>::node_type;

  private:
    uchar_type  m_chr    = '\0';
    size_type   m_exp    = 0;  // compressible
    loint_type  m_order  = 0;
    fnode_type* m_fnode  = nullptr;
    fnode_type* m_hnode  = nullptr;  // eliminatable
    offset_type m_hofst  = 0;  // compressible
    uint8_t     m_weight = 0;

  public:
    LData_Naive() = default;

    //! Default destructor
    virtual ~LData_Naive() = default;

    //! Copy constructor (deleted)
    LData_Naive(const LData_Naive&) = delete;

    //! Copy constructor (deleted)
    LData_Naive& operator=(const LData_Naive&) = delete;

    //! Move constructor
    LData_Naive(LData_Naive&&) noexcept = default;

    //! Move constructor
    LData_Naive& operator=(LData_Naive&&) noexcept = default;

    size_type get_memory_in_bytes() const {
        return sizeof(m_chr) + sizeof(m_exp) + sizeof(m_order) + sizeof(m_fnode) + sizeof(m_hnode) + sizeof(m_hofst);
    }

    inline void add_exp(size_type v) {
        m_exp += v;
    }
    inline void reset_hlink() {
        m_hnode = nullptr;
        m_hofst = 0;
    }
    inline void set_hlink(fnode_type* nd, offset_type of) {
        m_hnode = nd;
        m_hofst = of;
    }

    inline void set_chr(uchar_type v) {
        m_chr = v;
    }
    inline void set_exp(size_type v) {
        m_exp = v;
    }
    inline void set_order(loint_type v) {
        m_order = v;
    }
    inline void set_fnode(fnode_type* v) {
        m_fnode = v;
    }
    inline void set_hnode(fnode_type* v) {
        m_hnode = v;
    }
    inline void set_hofst(offset_type v) {
        m_hofst = v;
    }

    inline uchar_type get_chr() const {
        return m_chr;
    }
    inline size_type get_exp() const {
        return m_exp;
    }
    inline loint_type get_order() const {
        return m_order;
    }
    inline fnode_type* get_fnode() const {
        return m_fnode;
    }
    inline fnode_type* get_hnode() const {
        return m_hnode;
    }
    inline offset_type get_hofst() const {
        return m_hofst;
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
};

}  // namespace rcomp
