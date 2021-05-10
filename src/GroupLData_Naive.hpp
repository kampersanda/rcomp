/**
 * @file GroupLData_Naive.hpp
 * @author Shunsuke Kanda (kampersanda)
 */
#pragma once

#include <vector>

#include "GroupFData_Naive.hpp"
#include "LinkedList.hpp"

namespace rcomp {

/**
 * A class for naive data type of grouped L-node.
 *
 * @tparam t_GroupBound The upper bound of storable runs.
 * @tparam t_SumExpHint Explicitly keep the sum of exponents?
 * @tparam t_LookupHint Explicitly keep the mapping from run-ids to run-positions?
 */
template <size_type t_GroupBound, bool t_SumExpHint = true, bool t_LookupHint = true>
class GroupLData_Naive {
  public:
    using this_type  = GroupLData_Naive<t_GroupBound, t_SumExpHint, t_LookupHint>;
    using fnode_type = typename LinkedList<GroupFData_Naive<this_type>>::node_type;

    static constexpr auto GROUP_BOUND  = t_GroupBound;
    static constexpr auto SUMEXP_HINT  = t_SumExpHint;
    static constexpr auto LOOKUP_HINT  = t_LookupHint;
    static constexpr auto MAX_NUM_RUNS = GROUP_BOUND + 2;

    // Cursor
    class lcursor_type {
        friend class GroupLData_Naive;

        size_type m_lrnid = MAX_SIZE_INT;
        size_type m_lrpos = MAX_SIZE_INT;

      public:
        lcursor_type() = default;

        size_type get_lrnid() const {
            return m_lrnid;
        }
        size_type get_lrpos() const {
            return m_lrpos;
        }

      private:
        lcursor_type(size_type lrnid, size_type lrpos) : m_lrnid(lrnid), m_lrpos(lrpos) {}
    };

  private:
    struct unit_type {
        uchar_type  chr;  // of run
        size_type   exp;  // of run
        fnode_type* fnode;
        size_type   frnid;  // run-ID in the F-node
        size_type   lrnid;  // run-ID in the L-node (of this)
    };
    std::vector<unit_type> m_units;

    loint_type  m_order  = 0;
    fnode_type* m_hnode  = nullptr;
    offset_type m_hofst  = 0;
    uint8_t     m_weight = 0;

    // hints
    size_type m_sum_exps = 0;
    //
    std::array<size_type, MAX_NUM_RUNS> m_id_to_pos;

  public:
    GroupLData_Naive() = default;

    //! Default destructor
    virtual ~GroupLData_Naive() = default;

    //! Copy constructor (deleted)
    GroupLData_Naive(const GroupLData_Naive&) = delete;

    //! Copy constructor (deleted)
    GroupLData_Naive& operator=(const GroupLData_Naive&) = delete;

    //! Move constructor
    GroupLData_Naive(GroupLData_Naive&&) noexcept = default;

    //! Move constructor
    GroupLData_Naive& operator=(GroupLData_Naive&&) noexcept = default;

    void prefetch() const {}

    bool is_empty() const {
        return m_units.empty();
    }

    // Get the number of bytes
    size_type get_memory_in_bytes() const {
        size_type mem = 0;
        mem += get_size_info_memory_in_bytes();
        mem += get_link_info_memory_in_bytes();
        mem += get_unit_info_memory_in_bytes();
        mem += get_sumexp_hint_memory_in_bytes();
        mem += get_lookup_hint_memory_in_bytes();
        return mem;
    }

    size_type get_size_info_memory_in_bytes() const {
        return 0;
    }
    size_type get_link_info_memory_in_bytes() const {
        return sizeof(m_order) + sizeof(m_hnode) + sizeof(m_hofst) + sizeof(m_weight);
    }
    size_type get_unit_info_memory_in_bytes() const {
        return sizeof(unit_type) * m_units.capacity() + sizeof(m_units);
    }
    size_type get_sumexp_hint_memory_in_bytes() const {
        if constexpr (SUMEXP_HINT) {
            return sizeof(m_sum_exps);
        } else {
            return 0;
        }
    }
    size_type get_lookup_hint_memory_in_bytes() const {
        if constexpr (LOOKUP_HINT) {
            return sizeof(m_id_to_pos);
        } else {
            return 0;
        }
    }

    // Retrun the number of runs in the group
    size_type get_num_runs() const {
        return m_units.size();
    }
    size_type get_capa_runs() const {
        return m_units.capacity();
    }

    // Retrun the sum of exponents in the group
    size_type get_sum_exps() const {
        if constexpr (SUMEXP_HINT) {
            return m_sum_exps;
        } else {
            size_type sum_exps = 0;
            for (const auto& unit : m_units) {
                sum_exps += unit.exp;
            }
            return sum_exps;
        }
    }

    /**
     *  Accessor by L-run ID
     */

    lcursor_type get_lcursor(size_type lrnid) const {
        return {lrnid, get_lrpos(lrnid)};
    }
    lcursor_type get_lcursor_from_position(size_type lrpos) const {
        return {get_lrnid(lrpos), lrpos};
    }

    // Check if lcrsr has a valid id and position.
    bool is_valid_lcursor(const lcursor_type lcrsr) const {
        if (lcrsr.m_lrnid == MAX_SIZE_INT) return false;
        if (lcrsr.m_lrpos == MAX_SIZE_INT) return false;
        return lcrsr.m_lrpos == get_lrpos(lcrsr.m_lrnid);
    }

    lcursor_type get_first_lcursor() const {
        DEBUG_ABORT_IF(m_units.empty());
        return {m_units.front().lrnid, 0};
    }
    lcursor_type get_last_lcursor() const {
        DEBUG_ABORT_IF(m_units.empty());
        return {m_units.back().lrnid, m_units.size() - 1};
    }
    lcursor_type get_prev_lcursor(const lcursor_type lcrsr) const {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        DEBUG_ABORT_IF_LE(lcrsr.get_lrpos(), 0);
        return {m_units[lcrsr.get_lrpos() - 1].lrnid, lcrsr.get_lrpos() - 1};
    }
    lcursor_type get_next_lcursor(const lcursor_type lcrsr) const {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        DEBUG_ABORT_IF_LE(m_units.size() - 1, lcrsr.get_lrpos());
        return {m_units[lcrsr.get_lrpos() + 1].lrnid, lcrsr.get_lrpos() + 1};
    }
    bool is_first_lcursor(const lcursor_type lcrsr) const {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        return lcrsr.get_lrpos() == 0;
    }
    bool is_last_lcursor(const lcursor_type lcrsr) const {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        return lcrsr.get_lrpos() == m_units.size() - 1;
    }

    /**
     *  Getter & Setter for the group member
     */

    // Retrun the order lableled with list-label maintainance
    loint_type get_order() const {
        return m_order;
    }
    // Return the head-linked F-node
    fnode_type* get_hnode() const {
        return m_hnode;
    }
    // Return the offset to the head-linked F-node
    offset_type get_hofst() const {
        return m_hofst;
    }
    std::pair<fnode_type*, offset_type> get_hlink() const {
        return {m_hnode, m_hofst};
    }
    void set_order(loint_type order) {
        m_order = order;
    }
    void set_hnode(fnode_type* hnode) {
        m_hnode = hnode;
    }
    void set_hofst(offset_type hofst) {
        m_hofst = hofst;
    }
    void set_hlink(fnode_type* hnode, offset_type hofst) {
        m_hnode = hnode;
        m_hofst = hofst;
    }
    void reset_hlink() {
        m_hnode = nullptr;
        m_hofst = 0;
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

    /**
     *  Getter & Setter for each run unit
     */

    // Retrun the character of the run-unit with ID 'lrnid'.
    uchar_type get_chr(size_type lrnid) const {
        return get_chr(get_lcursor(lrnid));
    }
    // Retrun the exponent of the run-unit with ID 'lrnid'.
    size_type get_exp(size_type lrnid) const {
        return get_exp(get_lcursor(lrnid));
    }
    //
    void add_exp(size_type lrnid, size_type exp) {
        add_exp(get_lcursor(lrnid), exp);
    }
    // Return the pointer to F's run-unit corresponding to L's run-unit with ID 'lrnid'.
    std::pair<fnode_type*, size_type> get_frptr(size_type lrnid) const {
        return get_frptr(get_lcursor(lrnid));
    }
    //
    void set_frptr(size_type lrnid, fnode_type* new_fnode, size_type new_frnid) {
        set_frptr(get_lcursor(lrnid), new_fnode, new_frnid);
    }

    uchar_type get_chr(lcursor_type lcrsr) const {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        return m_units[lcrsr.get_lrpos()].chr;
    }
    size_type get_exp(lcursor_type lcrsr) const {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        return m_units[lcrsr.get_lrpos()].exp;
    }
    // Increment the exponent of the run-unit with ID 'lrnid' by 'exp'.
    void add_exp(lcursor_type lcrsr, size_type exp) {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        m_units[lcrsr.get_lrpos()].exp += exp;
        if constexpr (SUMEXP_HINT) {
            m_sum_exps += exp;
        }
    }
    std::pair<fnode_type*, size_type> get_frptr(lcursor_type lcrsr) const {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        return {m_units[lcrsr.get_lrpos()].fnode, m_units[lcrsr.get_lrpos()].frnid};
    }
    void set_frptr(lcursor_type lcrsr, fnode_type* new_fnode, size_type new_frnid) {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        m_units[lcrsr.get_lrpos()].fnode = new_fnode;
        m_units[lcrsr.get_lrpos()].frnid = new_frnid;
    }

    /**
     *  Updater
     */

    // Initillay insert a new L-run with (new_chr, new_exp)
    // Return the run-ID
    lcursor_type insert_init(const uchar_type new_chr, const size_type new_exp) {
        ABORT_IF(!is_empty());

        const size_type new_lrnid = make_free_lrnid();
        const size_type new_lrpos = m_units.size();

        m_units.push_back(unit_type{new_chr, new_exp, nullptr, 0, new_lrnid});

        if constexpr (SUMEXP_HINT) {
            m_sum_exps += new_exp;
        }

        if constexpr (LOOKUP_HINT) {
            std::fill_n(m_id_to_pos.data(), MAX_NUM_RUNS, MAX_SIZE_INT);
            m_id_to_pos[new_lrnid] = new_lrpos;
        }

        return {new_lrnid, new_lrpos};  // new run-ID
    }
    // Insert a new L-run with (new_chr, new_exp) before the L-run with lrnid.
    // Return the run-ID
    lcursor_type insert_before(lcursor_type& lcrsr, const uchar_type new_chr, const size_type new_exp) {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));

        const size_type new_lrnid = make_free_lrnid();
        const size_type new_lrpos = lcrsr.get_lrpos();

        m_units.insert(m_units.begin() + new_lrpos, unit_type{new_chr, new_exp, nullptr, 0, new_lrnid});
        lcrsr.m_lrpos += 1;

        if constexpr (SUMEXP_HINT) {
            m_sum_exps += new_exp;
        }

        if constexpr (LOOKUP_HINT) {
            DEBUG_ABORT_IF_NE(m_id_to_pos[new_lrnid], MAX_SIZE_INT);
            m_id_to_pos[new_lrnid] = new_lrpos;
            for (size_type i = new_lrpos + 1; i < m_units.size(); i++) {
                DEBUG_ABORT_IF_NE(m_id_to_pos[m_units[i].lrnid], i - 1);
                m_id_to_pos[m_units[i].lrnid] = i;
            }
        }

        return {new_lrnid, new_lrpos};
    }
    // Insert a new L-run with (new_chr, new_exp) after the L-run with lrnid.
    // Return the run-ID
    lcursor_type insert_after(lcursor_type& lcrsr, const uchar_type new_chr, const size_type new_exp) {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));

        const size_type new_lrnid = make_free_lrnid();
        const size_type new_lrpos = lcrsr.get_lrpos() + 1;

        m_units.insert(m_units.begin() + new_lrpos, unit_type{new_chr, new_exp, nullptr, 0, new_lrnid});

        if constexpr (SUMEXP_HINT) {
            m_sum_exps += new_exp;
        }

        if constexpr (LOOKUP_HINT) {
            DEBUG_ABORT_IF_NE(m_id_to_pos[new_lrnid], MAX_SIZE_INT);
            m_id_to_pos[new_lrnid] = new_lrpos;
            for (size_type i = new_lrpos + 1; i < m_units.size(); i++) {
                DEBUG_ABORT_IF_NE(m_id_to_pos[m_units[i].lrnid], i - 1);
                m_id_to_pos[m_units[i].lrnid] = i;
            }
        }

        return {new_lrnid, new_lrpos};
    }
    // Split the L-run with lrnid into two runs, whose the first L-run is of size new_exp1.
    // Return the ID of the second (i.e., new created) L-run.
    lcursor_type split_after(lcursor_type& lcrsr, const size_type new_exp_a) {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        ABORT_IF_OUT(new_exp_a, 1, m_units[lcrsr.get_lrpos()].exp - 1);

        const size_type new_lrnid = make_free_lrnid();
        const size_type new_lrpos = lcrsr.get_lrpos() + 1;

        unit_type& new_unit_a = m_units[lcrsr.get_lrpos()];
        unit_type  new_unit_b = {new_unit_a.chr, new_unit_a.exp - new_exp_a, nullptr, 0, new_lrnid};

        new_unit_a.exp = new_exp_a;
        m_units.insert(m_units.begin() + lcrsr.get_lrpos() + 1, new_unit_b);

        if constexpr (LOOKUP_HINT) {
            DEBUG_ABORT_IF_NE(m_id_to_pos[new_lrnid], MAX_SIZE_INT);
            m_id_to_pos[new_lrnid] = new_lrpos;
            for (size_type i = new_lrpos + 1; i < m_units.size(); i++) {
                DEBUG_ABORT_IF_NE(m_id_to_pos[m_units[i].lrnid], i - 1);
                m_id_to_pos[m_units[i].lrnid] = i;
            }
        }

        return {new_lrnid, new_lrpos};
    }

    /**
     *  Offset handlers
     */

    // Return the offset of run lrnid.
    offset_type get_offset(size_type lrnid) const {
        return get_offset(get_lcursor(lrnid));
    }
    offset_type get_offset(const lcursor_type lcrsr) const {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        offset_type lofst = 0;
        for (size_type i = 0; i < lcrsr.get_lrpos(); i++) {
            lofst += m_units[i].exp;
        }
        return lofst;
    }
    std::pair<lcursor_type, offset_type> access_with_offset(offset_type lofst) {
        const size_type num_runs = get_num_runs();
        for (size_type lrpos = 0; lrpos < num_runs; lrpos++) {
            if (lofst < offset_type(m_units[lrpos].exp)) {
                return {lcursor_type(m_units[lrpos].lrnid, lrpos), lofst};
            }
            lofst = lofst - m_units[lrpos].exp;
        }
        ABORT_IF(true);
        return {};
    }

    /**
     *  Splitter
     */

    // 後ろ半分のユニットを保存したデータを返し、それをこのグループから消す。
    // L-IDが変わるので注意が必要
    // 返すグループのorderやリンクはみ設定
    this_type divide_after() {
        ABORT_IF_LE(get_num_runs(), 1);
        return divide_after(get_num_runs() / 2);
    }
    this_type divide_after(const size_type num_runs_befo) {
        ABORT_IF_LE(get_num_runs(), num_runs_befo);
        const size_type num_runs_aftr = m_units.size() - num_runs_befo;

        // Clone the rear part
        this_type data_aftr;
        std::copy(m_units.begin() + num_runs_befo, m_units.end(), std::back_inserter(data_aftr.m_units));

        // Remove the rear part
        m_units.resize(num_runs_befo);

        for (size_type i = 0; i < num_runs_befo; i++) {
            m_units[i].lrnid = i;
        }
        for (size_type i = 0; i < num_runs_aftr; i++) {
            data_aftr.m_units[i].lrnid = i;
        }

        /**
         *  Update hints
         */
        if constexpr (SUMEXP_HINT) {
            m_sum_exps = 0;
            for (size_type i = 0; i < num_runs_befo; i++) {
                m_sum_exps += m_units[i].exp;
            }
            data_aftr.m_sum_exps = 0;
            for (size_type i = 0; i < num_runs_aftr; i++) {
                data_aftr.m_sum_exps += data_aftr.m_units[i].exp;
            }
        }

        if constexpr (LOOKUP_HINT) {
            std::fill_n(m_id_to_pos.data(), MAX_NUM_RUNS, MAX_SIZE_INT);
            for (size_type i = 0; i < num_runs_befo; i++) {
                m_id_to_pos[m_units[i].lrnid] = i;
            }
            std::fill_n(data_aftr.m_id_to_pos.data(), MAX_NUM_RUNS, MAX_SIZE_INT);
            for (size_type i = 0; i < num_runs_aftr; i++) {
                data_aftr.m_id_to_pos[data_aftr.m_units[i].lrnid] = i;
            }
        }

        return data_aftr;
    }

  private:
    //
    size_type make_free_lrnid() {
        if (m_units.empty()) {
            return 0;
        }
        size_type max_lrnid = 0;
        for (const auto& unit : m_units) {
            max_lrnid = std::max(max_lrnid, unit.lrnid);
        }
        return max_lrnid + 1;
    }

    /**
     *  ID <-> Pos
     */

    // Retrun the position of the run-unit with ID 'lrnid'.
    size_type get_lrpos(size_type lrnid) const {
        if constexpr (LOOKUP_HINT) {
            DEBUG_ABORT_IF_EQ(m_id_to_pos[lrnid], MAX_SIZE_INT);
            return m_id_to_pos[lrnid];
        } else {
            for (size_type lrpos = 0; lrpos < m_units.size(); lrpos++) {
                if (m_units[lrpos].lrnid == lrnid) {
                    return lrpos;
                }
            }
            return MAX_SIZE_INT;
        }
    }

    // Returns the ID of the run-unit located at 'lrpos'
    size_type get_lrnid(size_type lrpos) const {
        return m_units[lrpos].lrnid;
    }
};

}  // namespace rcomp
