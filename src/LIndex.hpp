/**
 * @file LIndex.hpp
 * @author Shunsuke Kanda (kampersanda)
 */
#pragma once

#include "LinkedList.hpp"
#include "abort_if.hpp"

namespace rcomp {

/**
 * A class for order maintenance data structure on L-nodes.
 *
 * @tparam t_LData The data type of L-node.
 */
template <class t_LData>
class LIndex {
  public:
    using this_type = LIndex<t_LData>;
    using list_type = LinkedList<t_LData>;

    using ldata_type = t_LData;
    using lnode_type = typename LinkedList<t_LData>::node_type;
    using fnode_type = typename ldata_type::fnode_type;

  private:
    list_type m_list;

  public:
    LIndex() = default;

    virtual ~LIndex() = default;

    //! Copy constructor (deleted)
    LIndex(const LIndex&) = delete;

    //! Copy constructor (deleted)
    LIndex& operator=(const LIndex&) = delete;

    //! Move constructor
    LIndex(LIndex&&) noexcept = default;

    //! Move constructor
    LIndex& operator=(LIndex&&) noexcept = default;

    void clear() {
        m_list.clear();

        lnode_type* x = m_list.get_head();
        x->get_data().set_chr(END_MARKER);
        x->get_data().set_exp(1);
        x->get_data().set_order(lo_common::LO_MIN_VALUE);
    }

    lnode_type* get_head() const {
        return m_list.get_head();
    }

    // FROM x1 -> x2 TO x1 -> x3 -> x2
    // RETURN x3
    lnode_type* insert_after(lnode_type* x1, uchar_type chr, size_type exp) {
        relabel(x1);

        lnode_type* x2 = x1->get_next();
        lnode_type* x3 = m_list.insert_after(x1);

        x3->get_data().set_chr(chr);
        x3->get_data().set_exp(exp);
        x3->get_data().set_order(get_middle_order(x1, x2));

        return x3;
    }

    size_type get_memory_in_bytes() const {
        return m_list.get_memory_in_bytes();
    }

    void show_memory_statistics() const {
        tfm::printfln("[Memory_LIndex]");
        tfm::reportfln("m_list:\t%d", m_list.get_memory_in_bytes());
    }

    void show_detailed_statistics() const {
        const auto [num_hlinks, num_lnodes] = compute_num_hlinks();
        tfm::printfln("[Detail_LIndex]");
        tfm::reportfln("num_hlinks:\t%d", num_hlinks);
        tfm::reportfln("num_lnodes:\t%d", num_lnodes);
        tfm::reportfln("ratio_hlinks:\t%g", num_hlinks / double(num_lnodes));
    }

    // not in O(1) time
    std::pair<size_type, size_type> compute_num_hlinks() const {
        size_type num_hlinks = 0;
        size_type num_lnodes = 0;

        const lnode_type* x = get_head()->get_next();
        for (; x != get_head(); x = x->get_next()) {
            if (x->get_data().get_hnode()) {
                num_hlinks += 1;
            }
            num_lnodes += 1;
        }

        return {num_hlinks, num_lnodes};
    }

    void test_order() const {
        const lnode_type* x = get_head();

        // Head L-node order is not LO_MIN_VALUE
        ABORT_IF_NE(lo_common::get_basic_order(x), lo_common::LO_MIN_VALUE);
        // Tail L-node order is not LO_MAX_VALUE
        ABORT_IF_NE(lo_common::get_star_order(x), lo_common::LO_MAX_VALUE);

        do {
            const lnode_type* y = x->get_next();
            // L-node orders are not sorted.
            ABORT_IF_LE(lo_common::get_star_order(y), lo_common::get_basic_order(x));
            x = y;
        } while (x != get_head());
    }

  private:
    void relabel(lnode_type* x) {
        lnode_type* y = x->get_next();

        // TODO: 密です (should be fixed)
        ABORT_IF(lo_common::get_star_order(y) == lo_common::LO_MAX_VALUE && lo_common::get_distance(x, y) == 1);

        loint_type j = 1;
        loint_type w = lo_common::get_distance(x, y);
        while (w <= j * j) {
            y = y->get_next();
            w = lo_common::get_distance(x, y);
            j += 1;
        }

        y = x->get_next();
        for (loint_type k = 1; k < j; k++) {
            loint_type new_v = ((w / j) * k) + x->get_data().get_order();
            y->get_data().set_order(new_v);
            y = y->get_next();
        }
    }

    loint_type get_middle_order(const lnode_type* x, const lnode_type* y) const {
        return lo_common::get_distance(x, y) / 2 + lo_common::get_basic_order(x);
    }
};

}  // namespace rcomp
