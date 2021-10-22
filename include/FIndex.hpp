/**
 * @file FIndex.hpp
 */
#pragma once

#include <set>

#include "LinkedList.hpp"
#include "abort_if.hpp"

namespace rcomp {

/**
 * A class for self-balancing binary search tree on F-nodes.
 *
 * @tparam t_FData The data type of F-node.
 */
template <class t_FData>
class FIndex {
    static_assert(sizeof(uchar_type) == 1);
    static_assert(END_MARKER == uchar_type(0));

  public:
    using this_type = FIndex<t_FData>;
    using list_type = LinkedList<t_FData>;

    using fdata_type = t_FData;
    using fnode_type = typename LinkedList<t_FData>::node_type;
    using lnode_type = typename fdata_type::lnode_type;

    struct comparer_type {
        bool operator()(const fnode_type* x, const fnode_type* y) const {
            return lo_common::get_basic_order(x) < lo_common::get_basic_order(y);
        }
    };

    using bst_type = std::set<const fnode_type*, comparer_type>;

  private:
    list_type m_list;
    bst_type m_bsts[256];

  public:
    //! Default constructor
    FIndex() = default;

    //! Default destructor
    virtual ~FIndex() = default;

    //! Copy constructor (deleted)
    FIndex(const FIndex&) = delete;

    //! Copy constructor (deleted)
    FIndex& operator=(const FIndex&) = delete;

    //! Move constructor
    FIndex(FIndex&&) noexcept = default;

    //! Move constructor
    FIndex& operator=(FIndex&&) noexcept = default;

    //! Get the head pointer of the list.
    fnode_type* get_head() const {
        return m_list.get_head();
    }

    //! Get the allocated memory in bytes.
    size_type get_memory_in_bytes(bool include_this = true) const {
        const size_type this_bytes = sizeof(*this) * include_this;
        return this_bytes + get_bsts_memory_in_bytes() + m_list.get_memory_in_bytes(false);
    }

    //! Get the allocated memory of BSTs in bytes.
    size_type get_bsts_memory_in_bytes() const {
        return std::accumulate(m_bsts, m_bsts + 256, 0ULL,
                               [&](size_type acc, const bst_type& bst) { return acc + get_bst_memory_in_bytes(bst); });
    }

    void show_memory_statistics() const {
        tfm::printfln("[Memory_FIndex]");
        tfm::reportfln("m_bsts:\t%d", get_bsts_memory_in_bytes());
        tfm::reportfln("m_list:\t%d", m_list.get_memory_in_bytes());
    }

    void show_detailed_statistics() const {
        const auto [num_tlinks, num_fnodes] = compute_num_tlinks();
        tfm::printfln("[Detail_FIndex]");
        tfm::reportfln("num_tlinks:\t%d", num_tlinks);
        tfm::reportfln("num_fnodes:\t%d", num_fnodes);
        tfm::reportfln("ratio_tlinks:\t%g", num_tlinks / double(num_fnodes));
    }

    void clear(lnode_type* y) {
        DEBUG_ABORT_IF_NE(END_MARKER, y->get_data().get_chr());

        m_list.clear();

        fnode_type* x = m_list.get_head();
        x->get_data().set_lnode(y);
        y->get_data().set_fnode(x);

        m_bsts[END_MARKER].insert(m_list.get_head());
    }

    // Update the list from "x1->x2" to "x1->x3->x2"
    // Return x3
    fnode_type* insert_after(fnode_type* x1, lnode_type* y3) {
        fnode_type* x3 = m_list.insert_after(x1);
        x3->get_data().set_lnode(y3);
        y3->get_data().set_fnode(x3);

        DEBUG_ABORT_IF_LT(x3->get_data().get_chr(), x1->get_data().get_chr());

        auto& bst = m_bsts[x3->get_data().get_chr()];
        DEBUG_ABORT_IF(bst.find(x3) != bst.end());
        bst.insert(x3);

        return x3;
    }

    // Search the most backward F-node x such that x->chr() <= chr and x->order() <= px->order().
    fnode_type* predecessor(uchar_type chr, const fnode_type* px) const {
        DEBUG_ABORT_IF_EQ(chr, END_MARKER);
        DEBUG_ABORT_IF(px == get_head());

        // Q: Why to construct BSTs for each character?
        // A: The predecessor query needs to input 'chr', but 'comparer_type' cannot handle it.

        if (m_bsts[chr].empty()) {
            while (--chr != END_MARKER) {
                if (!m_bsts[chr].empty()) {
                    break;
                }
            }
            DEBUG_ABORT_IF(m_bsts[chr].empty());
            return const_cast<fnode_type*>(*m_bsts[chr].crbegin());
        }

        auto itr = m_bsts[chr].upper_bound(px);
        if (itr == m_bsts[chr].end()) {
            return const_cast<fnode_type*>(*m_bsts[chr].crbegin());
        }
        return const_cast<fnode_type*>(*itr)->get_prev();
    }

    // Search the most backward F-node x such that x->chr() <= chr and x->order() <= px->order().
    fnode_type* predecessor_naive(uchar_type chr, const fnode_type* px) const {
        DEBUG_ABORT_IF(px == get_head());

        // 1) Search the head node in the F-interval of the character.
        const fnode_type* x = get_head();
        do {
            x = x->get_next();
            if (chr <= x->get_data().get_chr()) {
                break;
            }
        } while (x != get_head());

        // 2) Search the most frontward node such that px->order() < x->order()
        const loint_type pv = lo_common::get_basic_order(px);

        while (x != get_head()) {
            if (chr < x->get_data().get_chr()) {
                break;
            }
            if (pv < lo_common::get_basic_order(x)) {
                break;
            }
            x = x->get_next();
        }

        return x->get_prev();  // the previous node is the target
    }

    // Search the most backward F-node x such that x->chr() == chr and x->order() <= px->order().
    const fnode_type* exact_predecessor(uchar_type chr, const fnode_type* px) const {
        DEBUG_ABORT_IF_EQ(chr, END_MARKER);
        DEBUG_ABORT_IF_EQ(chr, px->get_data().get_chr());
        DEBUG_ABORT_IF(px == get_head());

        if (m_bsts[chr].empty()) {
            return nullptr;
        }
        if (lo_common::get_basic_order(px) < lo_common::get_basic_order(*m_bsts[chr].begin())) {
            // Then, the predecessor is none
            return nullptr;
        }
        if (lo_common::get_basic_order(px) > lo_common::get_basic_order(*m_bsts[chr].rbegin())) {
            // Then, the predecessor is the tail F-node with chr
            return *m_bsts[chr].rbegin();
        }

        auto itr = m_bsts[chr].upper_bound(px);
        return (*itr)->get_prev();
    }

    // Search the most frontward F-node x such that x->chr() == chr and px->order() <= x->order().
    const fnode_type* exact_successor(uchar_type chr, const fnode_type* px) const {
        DEBUG_ABORT_IF_EQ(chr, END_MARKER);
        DEBUG_ABORT_IF_EQ(chr, px->get_data().get_chr());
        DEBUG_ABORT_IF(px == get_head());

        if (m_bsts[chr].empty()) {
            return nullptr;
        }

        auto itr = m_bsts[chr].lower_bound(px);
        if (itr == m_bsts[chr].end()) {
            return nullptr;
        } else {
            return *itr;
        }
    }

    std::pair<size_type, size_type> compute_num_tlinks() const {
        size_type num_tlinks = 0;
        size_type num_fnodes = 0;

        const fnode_type* fnode = get_head();
        do {
            if (fnode->get_data().get_tnode()) {
                num_tlinks += 1;
            }
            num_fnodes += 1;
            fnode = fnode->get_next();
        } while (fnode != get_head());

        return {num_tlinks, num_fnodes};
    }

    void test_alphabet_order() const {
        const fnode_type* x = get_head();

        // Head F-node does not have END_MARKER.
        ABORT_IF_NE(x->get_data().get_chr(), END_MARKER);

        while (true) {
            const fnode_type* y = x->get_next();
            if (y == get_head()) {
                break;
            }
            // F-node characters are not sorted.
            ABORT_IF_LT(y->get_data().get_chr(), x->get_data().get_chr());
            x = y;
        };
    }
};

}  // namespace rcomp
