/**
 * @file GroupLIndex.hpp
 * @author Shunsuke Kanda (kampersanda)
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
class GroupLIndex {
    static_assert(sizeof(uchar_type) == 1);
    static_assert(END_MARKER == uchar_type(0));

  public:
    using this_type = GroupLIndex<t_LData>;
    using list_type = LinkedList<t_LData>;

    using ldata_type = t_LData;
    using lnode_type = typename LinkedList<ldata_type>::node_type;
    using fnode_type = typename ldata_type::fnode_type;

    using lcursor_type = typename lnode_type::data_type::lcursor_type;
    using fcursor_type = typename fnode_type::data_type::fcursor_type;

  private:
    list_type m_list;

  public:
    GroupLIndex() = default;

    virtual ~GroupLIndex() = default;

    //! Copy constructor (deleted)
    GroupLIndex(const GroupLIndex&) = delete;

    //! Copy constructor (deleted)
    GroupLIndex& operator=(const GroupLIndex&) = delete;

    //! Move constructor
    GroupLIndex(GroupLIndex&&) noexcept = default;

    //! Move constructor
    GroupLIndex& operator=(GroupLIndex&&) noexcept = default;

    lnode_type* get_head() const {
        return m_list.get_head();
    }

    std::pair<lnode_type*, lcursor_type> clear() {
        m_list.clear();

        lnode_type*  new_lnode = m_list.get_head();
        lcursor_type new_lcurs = new_lnode->get_data().insert_init(END_MARKER, 1);

        new_lnode->get_data().set_order(lo_common::LO_MIN_VALUE);

        return {new_lnode, new_lcurs};
    }

    // FROM x1 -> x2 TO x1 -> x3 -> x2
    // RETURN x3
    std::pair<lnode_type*, lcursor_type> insert_after(lnode_type* befo_lnode, uchar_type new_chr, size_type new_exp) {
        relabel(befo_lnode);

        lnode_type*  next_lnode = befo_lnode->get_next();
        lnode_type*  aftr_lnode = m_list.insert_after(befo_lnode);
        lcursor_type aftr_lcurs = aftr_lnode->get_data().insert_init(new_chr, new_exp);

        aftr_lnode->get_data().set_order(get_middle_order(befo_lnode, next_lnode));

        // lrnidの帰り地いる？
        return {aftr_lnode, aftr_lcurs};
    }

    lnode_type* divide_after(lnode_type* befo_lnode) {
        return divide_after(befo_lnode, befo_lnode->get_data().get_num_runs() / 2);
    }
    lnode_type* divide_after(lnode_type* befo_lnode, const lcursor_type lcrsr) {
        DEBUG_ABORT_IF(!befo_lnode->get_data().is_valid_lcursor(lcrsr));
        return divide_after(befo_lnode, lcrsr.get_lrpos() + 1);
    }
    lnode_type* divide_after(lnode_type* befo_lnode, const size_type num_befo_runs) {
        relabel(befo_lnode);

        lnode_type* next_lnode = befo_lnode->get_next();
        lnode_type* aftr_lnode = m_list.insert_after(befo_lnode);

        {
            ldata_type aftr_ldata = befo_lnode->get_data().divide_after(num_befo_runs);
            aftr_lnode->set_data(std::move(aftr_ldata));
        }

        // Reset F-pointers
        {
            ldata_type&  befo_ldata = befo_lnode->get_data();
            lcursor_type befo_lcrsr = befo_ldata.get_first_lcursor();

            while (true) {
                // DEBUG_PRINT(tfm::printfln("befo_lcrsr: (%d,%d)", befo_lcrsr.get_lrnid(), befo_lcrsr.get_lrpos()));

                auto [fnode, frnid] = befo_ldata.get_frptr(befo_lcrsr);
                DEBUG_ABORT_IF(!fnode->get_data().is_valid_fcursor(fnode->get_data().get_fcursor(frnid)));

                fnode->get_data().set_lrptr(frnid, befo_lnode, befo_lcrsr.get_lrnid());

                if (befo_ldata.is_last_lcursor(befo_lcrsr)) {
                    break;
                }
                befo_lcrsr = befo_ldata.get_next_lcursor(befo_lcrsr);
            }
        }
        {
            ldata_type&  aftr_ldata = aftr_lnode->get_data();
            lcursor_type aftr_lcrsr = aftr_ldata.get_first_lcursor();

            while (true) {
                // DEBUG_PRINT(tfm::printfln("aftr_lcrsr: (%d,%d)", aftr_lcrsr.get_lrnid(), aftr_lcrsr.get_lrpos()));

                auto [fnode, frnid] = aftr_ldata.get_frptr(aftr_lcrsr);
                DEBUG_ABORT_IF(!fnode->get_data().is_valid_fcursor(fnode->get_data().get_fcursor(frnid)));

                fnode->get_data().set_lrptr(frnid, aftr_lnode, aftr_lcrsr.get_lrnid());

                if (aftr_ldata.is_last_lcursor(aftr_lcrsr)) {
                    break;
                }
                aftr_lcrsr = aftr_ldata.get_next_lcursor(aftr_lcrsr);
            }
        }

        // Labeling
        aftr_lnode->get_data().set_order(get_middle_order(befo_lnode, next_lnode));

        return aftr_lnode;
    }

    size_type get_memory_in_bytes() const {
        return m_list.get_memory_in_bytes();
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
        size_type size_info   = 0;
        size_type link_info   = 0;
        size_type unit_info   = 0;
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
