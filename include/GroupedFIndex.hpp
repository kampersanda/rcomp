/**
 * @file GroupedFIndex.hpp
 */
#pragma once

#include <set>

#include "LinkedList.hpp"
#include "abort_if.hpp"

namespace rcomp {

// Sorry, very bad manner (should be improved with dynamic polymorphism and so on)
namespace dangerous_statics {

static loint_type q_order;
static size_type q_second_order;

}  // namespace dangerous_statics

namespace lo_common {

template <class T>
inline size_type get_second_order(const T* x) {
    return x->get_data().get_second_order();
}

}  // namespace lo_common

/**
 * A class for self-balancing binary search tree on the grouped F-nodes.
 *
 * @tparam t_FData The data type of the grouped F-node.
 * @attention The implementation resorts to static global variables in dangerous_statics.
 */
template <class t_FData>
class GroupedFIndex {
    static_assert(sizeof(uchar_type) == 1);
    static_assert(END_MARKER == uchar_type(0));

  public:
    using this_type = GroupedFIndex<t_FData>;
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
            const size_type x_second_order = x ? lo_common::get_second_order(x) : dangerous_statics::q_second_order;
            const size_type y_second_order = y ? lo_common::get_second_order(y) : dangerous_statics::q_second_order;
            return x_second_order < y_second_order;
        }
    };

    // For predecessor/successor queries
    using bst_type = std::set<const fnode_type*, comparer_type>;

  private:
    list_type m_list;
    bst_type m_bsts[256];

  public:
    //! Default constructor
    GroupedFIndex() = default;

    //! Default destructor
    virtual ~GroupedFIndex() = default;

    //! Copy constructor (deleted)
    GroupedFIndex(const GroupedFIndex&) = delete;

    //! Copy constructor (deleted)
    GroupedFIndex& operator=(const GroupedFIndex&) = delete;

    //! Move constructor
    GroupedFIndex(GroupedFIndex&&) noexcept = default;

    //! Move constructor
    GroupedFIndex& operator=(GroupedFIndex&&) noexcept = default;

    //! Get the head F-node.
    inline fnode_type* get_head() const {
        return m_list.get_head();
    }

    // Initialize the index by definiing the head F-node.
    inline std::tuple<fnode_type*, fcursor_type> clear(lnode_type* new_lnode, size_type new_lrnid) {
        m_list.clear();

        fnode_type* new_fnode = m_list.get_head();
        fcursor_type new_fcrsr = new_fnode->get_data().insert_init(new_lnode, new_lrnid);

        new_lnode->get_data().set_frptr(new_lrnid, new_fnode, new_fcrsr.get_frnid());
        m_bsts[END_MARKER].insert(new_fnode);

        return {new_fnode, new_fcrsr};
    }

    // Insert and set frptr
    std::tuple<fnode_type*, fcursor_type> insert_after(fnode_type* fnode, lnode_type* new_lnode, size_type new_lrnid) {
        fnode_type* new_fnode = m_list.insert_after(fnode);
        fcursor_type new_fcrsr = new_fnode->get_data().insert_init(new_lnode, new_lrnid);

        new_lnode->get_data().set_frptr(new_lrnid, new_fnode, new_fcrsr.get_frnid());
        DEBUG_ABORT_IF_LT(new_fnode->get_data().get_chr(), fnode->get_data().get_chr());

        auto& bst = m_bsts[new_fnode->get_data().get_chr()];
        DEBUG_ABORT_IF(bst.find(new_fnode) != bst.end());

        bst.insert(new_fnode);

        return {new_fnode, new_fcrsr};
    }

    //! Divide the F-node into two nodes
    inline fnode_type* divide_after(fnode_type* div_fnode) {
        return divide_after(div_fnode, div_fnode->get_data().get_num_runs() / 2);
    }

    //! Divide the F-node into two nodes
    inline fnode_type* divide_after(fnode_type* div_fnode, const fcursor_type& fcrsr) {
        DEBUG_ABORT_IF(!div_fnode->get_data().is_valid_fcursor(fcrsr));
        return divide_after(div_fnode, fcrsr.get_frpos() + 1);
    }

    //! Divide the F-node into two nodes such that the first part has 'num_befo' runs.
    inline fnode_type* divide_after(fnode_type* div_fnode, const size_type num_befo) {
        fnode_type* new_fnode = m_list.insert_after(div_fnode);

        {
            fdata_type new_fdata = div_fnode->get_data().divide_after(num_befo);
            new_fnode->set_data(std::move(new_fdata));
        }

        {
            fdata_type& div_fdata = div_fnode->get_data();
            fcursor_type div_fcrsr = div_fdata.get_first_fcursor();

            while (true) {
                auto [lnode, lrnid] = div_fdata.get_lrptr(div_fcrsr);
                DEBUG_ABORT_IF(!lnode->get_data().is_valid_lcursor(lnode->get_data().get_lcursor(lrnid)));

                lnode->get_data().set_frptr(lrnid, div_fnode, div_fcrsr.get_frnid());

                if (div_fdata.is_last_fcursor(div_fcrsr)) {
                    break;
                }
                div_fcrsr = div_fdata.get_next_fcursor(div_fcrsr);
            }
        }

        {
            fdata_type& new_fdata = new_fnode->get_data();
            fcursor_type new_fcrsr = new_fdata.get_first_fcursor();

            while (true) {
                auto [lnode, lrnid] = new_fdata.get_lrptr(new_fcrsr);
                DEBUG_ABORT_IF(!lnode->get_data().is_valid_lcursor(lnode->get_data().get_lcursor(lrnid)));

                lnode->get_data().set_frptr(lrnid, new_fnode, new_fcrsr.get_frnid());

                if (new_fdata.is_last_fcursor(new_fcrsr)) {
                    break;
                }
                new_fcrsr = new_fdata.get_next_fcursor(new_fcrsr);
            }
        }

        auto& bst = m_bsts[new_fnode->get_data().get_chr()];
        DEBUG_ABORT_IF(bst.find(new_fnode) != bst.end());
        bst.insert(new_fnode);

        return new_fnode;
    }

    // Search the most backward F-node x such that x->chr() <= chr and x->order() <= px->order().
    inline std::tuple<fnode_type*, fcursor_type>  //
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

        // Set the query variables
        dangerous_statics::q_order = q_order;
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
    std::tuple<fnode_type*, fcursor_type>  //
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

    //! Search the most backward F-node x such that
    //!  - q_chr == x->chr
    //!  - (x->order, x->second_order) <= (q_order, q_second_order)
    std::tuple<fnode_type*, fcursor_type>  //
    exact_predecessor(uchar_type q_chr, const loint_type q_order, const size_type q_second_order) const {
        DEBUG_ABORT_IF_EQ(q_chr, END_MARKER);

        if (m_bsts[q_chr].empty()) {
            return {nullptr, fcursor_type{}};
        }

        // Set the query variables
        dangerous_statics::q_order = q_order;
        dangerous_statics::q_second_order = q_second_order;

        // Query with the static variables
        auto itr = m_bsts[q_chr].upper_bound(nullptr);

        if (itr == m_bsts[q_chr].begin()) {
            // Then, (q_order, q_second_order) is the smallest
            return {nullptr, fcursor_type{}};
        }

        // Let x = itr->prev and y = itr.
        // (x->order, x->second_order) <= (q_order, q_second_order) < (y->order, y->second_order)

        if (itr != m_bsts[q_chr].end()) {
            fnode_type* fnode = const_cast<fnode_type*>(*itr)->get_prev();
            return {fnode, fnode->get_data().predecessor(q_order, q_second_order)};
        } else {
            fnode_type* fnode = const_cast<fnode_type*>(*m_bsts[q_chr].rbegin());
            return {fnode, fnode->get_data().predecessor(q_order, q_second_order)};
        }
    }

    //! Search the most frontward F-node y such that
    //!  - q_chr == y->chr
    //!  - (q_order, q_second_order) <= (y->order, y->second_order)
    std::tuple<fnode_type*, fcursor_type>  //
    exact_successor(uchar_type q_chr, const loint_type q_order, const size_type q_second_order) const {
        DEBUG_ABORT_IF_EQ(q_chr, END_MARKER);

        if (m_bsts[q_chr].empty()) {
            return {nullptr, fcursor_type{}};
        }

        // Set the query variables
        dangerous_statics::q_order = q_order;
        dangerous_statics::q_second_order = q_second_order;

        // Query with the static variables
        auto itr = m_bsts[q_chr].lower_bound(nullptr);

        // Let x = itr->prev and y = itr.
        // (x->order, x->second_order) < (q_order, q_second_order) <= (y->order, y->second_order)
        // That is, the target run may be the first of y or be in x.

        if (itr == m_bsts[q_chr].begin()) {
            fnode_type* fnode = const_cast<fnode_type*>(*itr);
            return {fnode, fnode->get_data().get_first_fcursor()};
        } else if (itr == m_bsts[q_chr].end()) {
            fnode_type* fnode = const_cast<fnode_type*>(*m_bsts[q_chr].rbegin());
            auto fcrsr_opt = fnode->get_data().successor(q_order, q_second_order);

            if (!fcrsr_opt) {
                return {nullptr, fcursor_type{}};
            }
            return {fnode, fcrsr_opt.value()};
        } else {
            // 1) The first is the target?
            fnode_type* fnode = const_cast<fnode_type*>(*itr);
            const loint_type y_order = lo_common::get_basic_order(fnode);
            const size_type y_second_order = fnode->get_data().get_second_order();
            if (y_order == q_order and y_second_order == q_second_order) {
                return {fnode, fnode->get_data().get_first_fcursor()};
            }
            // 2) The target may be in the prev group.
            fnode = fnode->get_prev();
            auto fcrsr_opt = fnode->get_data().successor(q_order, q_second_order);

            if (fcrsr_opt) {
                return {fnode, fcrsr_opt.value()};
            }

            // 3) The target was not in the prev group.
            fnode = const_cast<fnode_type*>(*itr);
            return {fnode, fnode->get_data().get_first_fcursor()};
        }
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

        const auto [size_info, link_info, unit_info, sumexp_hint, lookup_hint] = compute_flist_memory();

        tfm::printfln("[Memory_FList]");
        tfm::reportfln("size_info:\t%d", size_info);
        tfm::reportfln("link_info:\t%d", link_info);
        tfm::reportfln("unit_info:\t%d", unit_info);
        tfm::reportfln("sumexp_hint:\t%d", sumexp_hint);
        tfm::reportfln("lookup_hint:\t%d", lookup_hint);
    }

    std::tuple<size_type, size_type, size_type, size_type, size_type> compute_flist_memory() const {
        size_type size_info = 0;
        size_type link_info = 0;
        size_type unit_info = 0;
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
