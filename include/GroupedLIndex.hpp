/**
 * @file GroupedLIndex.hpp
 */
#pragma once

#include "LinkedList.hpp"
#include "abort_if.hpp"

namespace rcomp {

/**
 * A class for order maintenance data structure on grouped L-nodes.
 *
 * @tparam t_LData The data type of grouped L-node.
 */
template <class t_LData>
class GroupedLIndex {
    static_assert(sizeof(uchar_type) == 1);
    static_assert(END_MARKER == uchar_type(0));

  public:
    using this_type = GroupedLIndex<t_LData>;
    using list_type = LinkedList<t_LData>;

    using ldata_type = t_LData;
    using lnode_type = typename LinkedList<ldata_type>::node_type;
    using fnode_type = typename ldata_type::fnode_type;

    using lcursor_type = typename lnode_type::data_type::lcursor_type;
    using fcursor_type = typename fnode_type::data_type::fcursor_type;

  private:
    list_type m_list;

  public:
    //! Default constructor
    GroupedLIndex() = default;

    //! Default destructor
    virtual ~GroupedLIndex() = default;

    //! Copy constructor (deleted)
    GroupedLIndex(const GroupedLIndex&) = delete;

    //! Copy constructor (deleted)
    GroupedLIndex& operator=(const GroupedLIndex&) = delete;

    //! Move constructor
    GroupedLIndex(GroupedLIndex&&) noexcept = default;

    //! Move constructor
    GroupedLIndex& operator=(GroupedLIndex&&) noexcept = default;

    //! Get the head L-node.
    inline lnode_type* get_head() const {
        return m_list.get_head();
    }

    //!
    inline std::tuple<lnode_type*, lcursor_type> clear() {
        m_list.clear();

        lnode_type* new_lnode = m_list.get_head();
        lcursor_type new_lcurs = new_lnode->get_data().insert_init(END_MARKER, 1);

        new_lnode->get_data().set_order(lo_common::LO_MIN_VALUE);

        return {new_lnode, new_lcurs};
    }

    // FROM x1 -> x2 TO x1 -> x3 -> x2
    // RETURN x3
    std::tuple<lnode_type*, lcursor_type> insert_after(lnode_type* lnode, uchar_type new_chr, size_type new_exp) {
        relabel(lnode);

        lnode_type* nxt_lnode = lnode->get_next();
        lnode_type* new_lnode = m_list.insert_after(lnode);
        lcursor_type new_lcurs = new_lnode->get_data().insert_init(new_chr, new_exp);

        new_lnode->get_data().set_order(get_middle_order(lnode, nxt_lnode));

        return {new_lnode, new_lcurs};
    }

    inline lnode_type* divide_after(lnode_type* div_fnode) {
        return divide_after(div_fnode, div_fnode->get_data().get_num_runs() / 2);
    }

    inline lnode_type* divide_after(lnode_type* div_fnode, const lcursor_type& lcrsr) {
        DEBUG_ABORT_IF(!div_fnode->get_data().is_valid_lcursor(lcrsr));
        return divide_after(div_fnode, lcrsr.get_lrpos() + 1);
    }

    inline lnode_type* divide_after(lnode_type* div_fnode, const size_type num_befo) {
        relabel(div_fnode);

        lnode_type* nxt_lnode = div_fnode->get_next();
        lnode_type* new_lnode = m_list.insert_after(div_fnode);

        {
            ldata_type new_ldata = div_fnode->get_data().divide_after(num_befo);
            new_lnode->set_data(std::move(new_ldata));
        }

        {
            ldata_type& div_ldata = div_fnode->get_data();
            lcursor_type div_lcrsr = div_ldata.get_first_lcursor();

            while (true) {
                auto [fnode, frnid] = div_ldata.get_frptr(div_lcrsr);
                DEBUG_ABORT_IF(!fnode->get_data().is_valid_fcursor(fnode->get_data().get_fcursor(frnid)));

                fnode->get_data().set_lrptr(frnid, div_fnode, div_lcrsr.get_lrnid());

                if (div_ldata.is_last_lcursor(div_lcrsr)) {
                    break;
                }
                div_lcrsr = div_ldata.get_next_lcursor(div_lcrsr);
            }
        }

        {
            ldata_type& new_ldata = new_lnode->get_data();
            lcursor_type new_lcrsr = new_ldata.get_first_lcursor();

            while (true) {
                auto [fnode, frnid] = new_ldata.get_frptr(new_lcrsr);
                DEBUG_ABORT_IF(!fnode->get_data().is_valid_fcursor(fnode->get_data().get_fcursor(frnid)));

                fnode->get_data().set_lrptr(frnid, new_lnode, new_lcrsr.get_lrnid());

                if (new_ldata.is_last_lcursor(new_lcrsr)) {
                    break;
                }
                new_lcrsr = new_ldata.get_next_lcursor(new_lcrsr);
            }
        }

        // Labeling
        new_lnode->get_data().set_order(get_middle_order(div_fnode, nxt_lnode));

        return new_lnode;
    }

    //! Get the allocated memory in bytes.
    size_type get_memory_in_bytes(bool include_this = true) const {
        const size_type this_bytes = sizeof(*this) * include_this;
        return this_bytes + m_list.get_memory_in_bytes(false);
    }

    void show_memory_statistics() const {
        tfm::printfln("[Memory_LIndex]");
        tfm::reportfln("m_list:\t%d", m_list.get_memory_in_bytes());

        const auto [size_info, link_info, unit_info, sumexp_hint, lookup_hint] = compute_llist_memory();

        tfm::printfln("[Memory_LList]");
        tfm::reportfln("size_info:\t%d", size_info);
        tfm::reportfln("link_info:\t%d", link_info);
        tfm::reportfln("unit_info:\t%d", unit_info);
        tfm::reportfln("sumexp_hint:\t%d", sumexp_hint);
        tfm::reportfln("lookup_hint:\t%d", lookup_hint);
    }

    std::tuple<size_type, size_type, size_type, size_type, size_type> compute_llist_memory() const {
        size_type size_info = 0;
        size_type link_info = 0;
        size_type unit_info = 0;
        size_type sumexp_hint = 0;
        size_type lookup_hint = 0;

        const lnode_type* lnode = get_head()->get_next();
        for (; lnode != get_head(); lnode = lnode->get_next()) {
            size_info += lnode->get_data().get_size_info_memory_in_bytes();
            link_info += lnode->get_data().get_link_info_memory_in_bytes();
            unit_info += lnode->get_data().get_unit_info_memory_in_bytes();
            sumexp_hint += lnode->get_data().get_sumexp_hint_memory_in_bytes();
            lookup_hint += lnode->get_data().get_lookup_hint_memory_in_bytes();
        }

        return {size_info, link_info, unit_info, sumexp_hint, lookup_hint};
    }

    void show_detailed_statistics() const {
        const auto [num_hlinks, num_lnodes, size_units, capa_units] = compute_lnode_stats();
        tfm::printfln("[Detail_LIndex]");
        tfm::reportfln("num_hlinks:\t%d", num_hlinks);
        tfm::reportfln("num_lnodes:\t%d", num_lnodes);
        tfm::reportfln("ratio_hlinks:\t%g", num_hlinks / double(num_lnodes));
        tfm::reportfln("size_units:\t%d", size_units);
        tfm::reportfln("capa_units:\t%d", capa_units);
        tfm::reportfln("ave_size_units:\t%g", size_units / double(num_lnodes));
        tfm::reportfln("ave_capa_units:\t%g", capa_units / double(num_lnodes));
    }

    std::tuple<size_type, size_type, size_type, size_type> compute_lnode_stats() const {
        size_type num_hlinks = 0;
        size_type num_lnodes = 0;
        size_type size_units = 0;
        size_type capa_units = 0;

        const lnode_type* lnode = get_head()->get_next();
        for (; lnode != get_head(); lnode = lnode->get_next()) {
            if (lnode->get_data().get_hnode()) {
                num_hlinks += 1;
            }
            num_lnodes += 1;
            size_units += lnode->get_data().get_num_runs();
            capa_units += lnode->get_data().get_capa_runs();
        }

        return {num_hlinks, num_lnodes, size_units, capa_units};
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
