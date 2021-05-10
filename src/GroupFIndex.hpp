/**
 * @file GroupFIndex.hpp
 * @author Shunsuke Kanda (kampersanda)
 */
#pragma once

#include <set>

#include "LinkedList.hpp"
#include "abort_if.hpp"

namespace rcomp {

// Very bad manner (should be improved with dynamic polymorphism)
namespace dangerous_statics {

static loint_type q_order;
static size_type  q_second_order;

}  // namespace dangerous_statics

/**
 * A class for self-balancing binary search tree on grouped F-nodes.
 *
 * @tparam t_FData The data type of grouped F-node.
 */
template <class t_FData>
class GroupFIndex {
    static_assert(sizeof(uchar_type) == 1);
    static_assert(END_MARKER == uchar_type(0));

  public:
    using this_type = GroupFIndex<t_FData>;
    using list_type = LinkedList<t_FData>;

    using fdata_type = t_FData;
    using fnode_type = typename LinkedList<fdata_type>::node_type;
    using lnode_type = typename fdata_type::lnode_type;

    using lcursor_type = typename lnode_type::data_type::lcursor_type;
    using fcursor_type = typename fnode_type::data_type::fcursor_type;

    struct comparer_type {
        // Compare the F-nodes in the total order in L
        bool operator()(const fnode_type* x, const fnode_type* y) const {
            DEBUG_ABORT_IF(!x and !y);

            // First, compare with the L-node order
            const loint_type x_order = x ? lo_common::get_basic_order(x) : dangerous_statics::q_order;
            const loint_type y_order = y ? lo_common::get_basic_order(y) : dangerous_statics::q_order;
            if (x_order != y_order) {
                return x_order < y_order;
            }

            // Second, compare with the offset position in the L-node
            const size_type x_second_order = x ? x->get_data().get_second_order() : dangerous_statics::q_second_order;
            const size_type y_second_order = y ? y->get_data().get_second_order() : dangerous_statics::q_second_order;
            return x_second_order < y_second_order;
        }
    };

    using bst_type = std::set<const fnode_type*, comparer_type>;

  private:
    list_type m_list;
    bst_type  m_bsts[256];

  public:
    GroupFIndex() = default;

    virtual ~GroupFIndex() = default;

    //! Copy constructor (deleted)
    GroupFIndex(const GroupFIndex&) = delete;

    //! Copy constructor (deleted)
    GroupFIndex& operator=(const GroupFIndex&) = delete;

    //! Move constructor
    GroupFIndex(GroupFIndex&&) noexcept = default;

    //! Move constructor
    GroupFIndex& operator=(GroupFIndex&&) noexcept = default;

    fnode_type* get_head() const {
        return m_list.get_head();
    }

    // Initialize the index by definiing the head F-node.
    std::pair<fnode_type*, fcursor_type> clear(lnode_type* new_lnode, size_type new_lrnid) {
        m_list.clear();

        fnode_type*  new_fnode = m_list.get_head();
        fcursor_type new_fcrsr = new_fnode->get_data().insert_init(new_lnode, new_lrnid);

        new_lnode->get_data().set_frptr(new_lrnid, new_fnode, new_fcrsr.get_frnid());
        m_bsts[END_MARKER].insert(new_fnode);

        return {new_fnode, new_fcrsr};
    }

    // Insert and set frptr
    std::pair<fnode_type*, fcursor_type> insert_after(fnode_type* fnode, lnode_type* new_lnode, size_type new_lrnid) {
        fnode_type*  new_fnode = m_list.insert_after(fnode);
        fcursor_type new_fcrsr = new_fnode->get_data().insert_init(new_lnode, new_lrnid);

        new_lnode->get_data().set_frptr(new_lrnid, new_fnode, new_fcrsr.get_frnid());
        DEBUG_ABORT_IF_LT(new_fnode->get_data().get_chr(), fnode->get_data().get_chr());

        auto& bst = m_bsts[new_fnode->get_data().get_chr()];
        DEBUG_ABORT_IF(bst.find(new_fnode) != bst.end());

        bst.insert(new_fnode);

        return {new_fnode, new_fcrsr};
    }

    fnode_type* divide_after(fnode_type* befo_fnode) {
        return divide_after(befo_fnode, befo_fnode->get_data().get_num_runs() / 2);
    }
    fnode_type* divide_after(fnode_type* befo_fnode, const fcursor_type fcrsr) {
        DEBUG_ABORT_IF(!befo_fnode->get_data().is_valid_fcursor(fcrsr));
        return divide_after(befo_fnode, fcrsr.get_frpos() + 1);
    }
    fnode_type* divide_after(fnode_type* befo_fnode, const size_type num_befo_runs) {
        fnode_type* aftr_fnode = m_list.insert_after(befo_fnode);
        {
            fdata_type aftr_fdata = befo_fnode->get_data().divide_after(num_befo_runs);
            aftr_fnode->set_data(std::move(aftr_fdata));
        }

        // Reset L-pointers (should be modified)
        {
            fdata_type&  befo_fdata = befo_fnode->get_data();
            fcursor_type befo_fcrsr = befo_fdata.get_first_fcursor();

            while (true) {
                auto [lnode, lrnid] = befo_fdata.get_lrptr(befo_fcrsr);
                DEBUG_ABORT_IF(!lnode->get_data().is_valid_lcursor(lnode->get_data().get_lcursor(lrnid)));

                lnode->get_data().set_frptr(lrnid, befo_fnode, befo_fcrsr.get_frnid());

                if (befo_fdata.is_last_fcursor(befo_fcrsr)) {
                    break;
                }
                befo_fcrsr = befo_fdata.get_next_fcursor(befo_fcrsr);
            }
        }

        {
            fdata_type&  aftr_fdata = aftr_fnode->get_data();
            fcursor_type aftr_fcrsr = aftr_fdata.get_first_fcursor();

            while (true) {
                auto [lnode, lrnid] = aftr_fdata.get_lrptr(aftr_fcrsr);
                DEBUG_ABORT_IF(!lnode->get_data().is_valid_lcursor(lnode->get_data().get_lcursor(lrnid)));

                lnode->get_data().set_frptr(lrnid, aftr_fnode, aftr_fcrsr.get_frnid());

                if (aftr_fdata.is_last_fcursor(aftr_fcrsr)) {
                    break;
                }
                aftr_fcrsr = aftr_fdata.get_next_fcursor(aftr_fcrsr);
            }
        }

        auto& bst = m_bsts[aftr_fnode->get_data().get_chr()];
        DEBUG_ABORT_IF(bst.find(aftr_fnode) != bst.end());
        bst.insert(aftr_fnode);

        return aftr_fnode;
    }

    // Search the most backward F-node x such that x->chr() <= chr and x->order() <= px->order().
    std::pair<fnode_type*, fcursor_type>  //
    predecessor(uchar_type chr, const loint_type q_order, const size_type q_second_order) const {
        DEBUG_ABORT_IF_EQ(chr, END_MARKER);

        // Q: Why to construct BSTs for each character?
        // A: The predecessor query needs to input 'chr', but 'comparer_type' cannot handle it.

        if (m_bsts[chr].empty()) {
            while (--chr != END_MARKER) {
                if (!m_bsts[chr].empty()) {
                    break;
                }
            }
            DEBUG_ABORT_IF(m_bsts[chr].empty());

            fnode_type* fnode = const_cast<fnode_type*>(*m_bsts[chr].crbegin());
            return {fnode, fnode->get_data().get_last_fcursor()};
        }

        dangerous_statics::q_order        = q_order;
        dangerous_statics::q_second_order = q_second_order;

        // Query with the static variables
        auto itr = m_bsts[chr].upper_bound(nullptr);

        if (itr == m_bsts[chr].begin()) {
            // Then, the target position is the most precedent in all groups with the character.
            fnode_type* fnode = const_cast<fnode_type*>(*itr)->get_prev();
            return {fnode, fnode->get_data().get_last_fcursor()};
        }

        if (itr != m_bsts[chr].end()) {
            // Then, the target position is in the previous group.
            fnode_type* fnode = const_cast<fnode_type*>(*itr)->get_prev();
            return {fnode, fnode->get_data().predecessor(q_order, q_second_order)};
        } else {
            // Then, the target position is in the last group.
            fnode_type* fnode = const_cast<fnode_type*>(*m_bsts[chr].crbegin());
            return {fnode, fnode->get_data().predecessor(q_order, q_second_order)};
        }
    }

    // Search the most backword F-node and run-unit whose order <= q_order and lrpos <= q_second_order.
    // Return its pointer and run-ID.
    std::pair<fnode_type*, fcursor_type>  //
    predecessor_naive(const uchar_type chr, const loint_type q_order, const size_type q_second_order) const {
        DEBUG_ABORT_IF_EQ(chr, END_MARKER);

        // 1) Search the first node within the F-intervals of 'chr'.
        const fnode_type* x = get_head();
        do {
            x = x->get_next();
            if (chr <= x->get_data().get_chr()) {
                break;
            }
        } while (x != get_head());

        // 2) Search the upper bound with comparer_type
        while (x != get_head()) {
            if (chr < x->get_data().get_chr()) {
                break;
            }
            if (q_order < lo_common::get_basic_order(x)) {
                break;
            }
            if (q_order == lo_common::get_basic_order(x) and q_second_order < x->get_data().get_second_order()) {
                break;
            }
            x = x->get_next();
        }

        x = x->get_prev();  // the previous node is the target
        DEBUG_ABORT_IF_LT(chr, x->get_data().get_chr());

        if (x->get_data().get_chr() < chr) {
            return {const_cast<fnode_type*>(x), x->get_data().get_last_fcursor()};
        } else {  // x->get_data().get_chr() == chr
            return {const_cast<fnode_type*>(x), x->get_data().predecessor(q_order, q_second_order)};
        }
    }

    size_type get_memory_in_bytes() const {
        return get_bst_memory_in_bytes() + m_list.get_memory_in_bytes();
    }

    // https://ny23.hatenadiary.org/entry/20111129/p1
    size_type get_bst_memory_in_bytes() const {
        size_type mem = 0;
        for (size_type i = 0; i < 256; i++) {
            mem += m_bsts[i].size() * (sizeof(int) +  // Node color (enum; int)
                                       3 * sizeof(void*) +  // Pointers to parent, left, and right
                                       sizeof(const fnode_type*)  // Key
                                      );
        }
        return mem;
    }

    void show_memory_statistics() const {
        tfm::printfln("[Memory_FIndex]");
        tfm::reportfln("m_bsts:\t%d", get_bst_memory_in_bytes());
        tfm::reportfln("m_list:\t%d", m_list.get_memory_in_bytes());

        const auto [size_info, link_info, unit_info, sumexp_hint, lookup_hint] = compute_flist_memory();

        tfm::printfln("[Memory_FList]");
        tfm::reportfln("size_info:\t%d", size_info);
        tfm::reportfln("link_info:\t%d", link_info);
        tfm::reportfln("unit_info:\t%d", unit_info);
        tfm::reportfln("sumexp_hint:\t%d", sumexp_hint);
        tfm::reportfln("lookup_hint:\t%d", lookup_hint);
    }

    std::tuple<size_type, size_type, size_type, size_type, size_type> compute_flist_memory() const {
        size_type size_info   = 0;
        size_type link_info   = 0;
        size_type unit_info   = 0;
        size_type sumexp_hint = 0;
        size_type lookup_hint = 0;

        const fnode_type* fnode = get_head();
        do {
            size_info += fnode->get_data().get_size_info_memory_in_bytes();
            link_info += fnode->get_data().get_link_info_memory_in_bytes();
            unit_info += fnode->get_data().get_unit_info_memory_in_bytes();
            sumexp_hint += fnode->get_data().get_sumexp_hint_memory_in_bytes();
            lookup_hint += fnode->get_data().get_lookup_hint_memory_in_bytes();
            fnode = fnode->get_next();
        } while (fnode != get_head());

        return {size_info, link_info, unit_info, sumexp_hint, lookup_hint};
    }

    void show_detailed_statistics() const {
        const auto [num_tlinks, num_fnodes, size_units, capa_units] = compute_fnode_stats();
        tfm::printfln("[Detail_FIndex]");
        tfm::reportfln("num_tlinks:\t%d", num_tlinks);
        tfm::reportfln("num_fnodes:\t%d", num_fnodes);
        tfm::reportfln("ratio_tlinks:\t%g", num_tlinks / double(num_fnodes));
        tfm::reportfln("size_units:\t%d", size_units);
        tfm::reportfln("capa_units:\t%d", capa_units);
        tfm::reportfln("ave_size_units:\t%g", size_units / double(num_fnodes));
        tfm::reportfln("ave_capa_units:\t%g", capa_units / double(num_fnodes));
    }

    std::tuple<size_type, size_type, size_type, size_type> compute_fnode_stats() const {
        size_type num_tlinks = 0;
        size_type num_fnodes = 0;
        size_type size_units = 0;
        size_type capa_units = 0;

        const fnode_type* fnode = get_head();
        do {
            if (fnode->get_data().get_tnode()) {
                num_tlinks += 1;
            }
            num_fnodes += 1;
            size_units += fnode->get_data().get_num_runs();
            capa_units += fnode->get_data().get_capa_runs();
            fnode = fnode->get_next();
        } while (fnode != get_head());

        return {num_tlinks, num_fnodes, size_units, capa_units};
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
