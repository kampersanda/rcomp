/**
 * @file GroupFData_Naive.hpp
 * @author Shunsuke Kanda (kampersanda)
 */
#pragma once

#include <vector>

#include "LinkedList.hpp"

namespace rcomp {

/**
 * A class for naive data type of grouped F-node.
 *
 * @tparam t_LData The data type of grouped L-node.
 */
template <class t_LData>
class GroupFData_Naive {
  public:
    using this_type  = GroupFData_Naive<t_LData>;
    using lnode_type = typename LinkedList<t_LData>::node_type;

    static constexpr size_type GROUP_BOUND = t_LData::GROUP_BOUND;
    static constexpr bool      SUMEXP_HINT = t_LData::SUMEXP_HINT;
    static constexpr bool      LOOKUP_HINT = t_LData::LOOKUP_HINT;

    static constexpr size_type MAX_NUM_RUNS = GROUP_BOUND + 2;

    // Cursor
    class fcursor_type {
        friend class GroupFData_Naive;

        size_type m_frnid = MAX_SIZE_INT;
        size_type m_frpos = MAX_SIZE_INT;

      public:
        fcursor_type() = default;

        size_type get_frnid() const {
            return m_frnid;
        }
        size_type get_frpos() const {
            return m_frpos;
        }

      private:
        fcursor_type(size_type frnid, size_type frpos) : m_frnid(frnid), m_frpos(frpos) {}
    };

  private:
    struct unit_type {
        lnode_type* lnode;
        size_type   lrnid;  // run-ID in the L-node
        size_type   frnid;  // run-ID in the F-node (of this)
    };
    std::vector<unit_type> m_units;

    lnode_type* m_tnode  = nullptr;
    uint8_t     m_weight = 0;

    // hints
    size_type m_sum_exps = 0;
    //
    std::array<size_type, MAX_NUM_RUNS> m_id_to_pos;

  public:
    GroupFData_Naive() = default;

    virtual ~GroupFData_Naive() = default;

    //! Copy constructor (deleted)
    GroupFData_Naive(const GroupFData_Naive&) = delete;

    //! Copy constructor (deleted)
    GroupFData_Naive& operator=(const GroupFData_Naive&) = delete;

    //! Move constructor
    GroupFData_Naive(GroupFData_Naive&&) noexcept = default;

    //! Move constructor
    GroupFData_Naive& operator=(GroupFData_Naive&&) noexcept = default;

    void prefetch() const {}

    bool is_empty() const {
        return m_units.empty();
    }
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
        return sizeof(m_tnode) + sizeof(m_weight);
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
                sum_exps += unit.lnode->get_data().get_exp(unit.lrnid);
            }
            return sum_exps;
        }
    }

    /**
     *  Accessor by L-run ID
     */

    fcursor_type get_fcursor(size_type frnid) const {
        return {frnid, get_frpos(frnid)};
    }

    fcursor_type get_fcursor_from_position(size_type frpos) const {
        return {get_frnid(frpos), frpos};
    }

    bool is_valid_fcursor(const fcursor_type fcrsr) const {
        if (fcrsr.m_frnid == MAX_SIZE_INT) return false;
        if (fcrsr.m_frpos == MAX_SIZE_INT) return false;
        return fcrsr.m_frpos == get_frpos(fcrsr.m_frnid);
    }

    fcursor_type get_first_fcursor() const {
        DEBUG_ABORT_IF(m_units.empty());
        return {m_units.front().frnid, 0};
    }
    fcursor_type get_last_fcursor() const {
        DEBUG_ABORT_IF(m_units.empty());
        return {m_units.back().frnid, m_units.size() - 1};
    }
    fcursor_type get_prev_fcursor(const fcursor_type fcrsr) const {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        DEBUG_ABORT_IF_LE(fcrsr.get_frpos(), 0);
        return {m_units[fcrsr.get_frpos() - 1].frnid, fcrsr.get_frpos() - 1};
    }
    fcursor_type get_next_fcursor(const fcursor_type fcrsr) const {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        DEBUG_ABORT_IF_LE(m_units.size() - 1, fcrsr.get_frpos());
        return {m_units[fcrsr.get_frpos() + 1].frnid, fcrsr.get_frpos() + 1};
    }

    bool is_first_fcursor(const fcursor_type fcrsr) const {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        return fcrsr.get_frpos() == 0;
    }
    bool is_last_fcursor(const fcursor_type fcrsr) const {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        return fcrsr.get_frpos() == m_units.size() - 1;
    }

    /**
     *  Getter & Setter for the group member
     */
    lnode_type* get_tnode() const {
        return m_tnode;
    }
    offset_type get_tofst() const {
        return m_tnode->get_data().get_hofst();
    }
    void set_tnode(lnode_type* tnode) {
        m_tnode = tnode;
    }
    void reset_tlink() {
        m_tnode = nullptr;
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
     * For predecesspr query
     */

    // Retrun the order of the L-node containing the first run-unit of this group.
    loint_type get_order() const {
        DEBUG_ABORT_IF(m_units.empty());
        return m_units.front().lnode->get_data().get_order();
    }
    // Retrun the run-unit position in the L-node for the first run-unit of this group.
    size_type get_second_order() const {  // second order?
        DEBUG_ABORT_IF(m_units.empty());
        return m_units.front().lnode->get_data().get_lcursor(m_units.front().lrnid).get_lrpos();
    }
    // Search the most backword run-unit whose order <= q_order and second_order <= q_second_order.
    // Return its run-ID.
    fcursor_type predecessor(const loint_type q_order, const size_type q_second_order) const {
        size_type i = 0;
        for (; i < m_units.size(); i++) {
            const unit_type& unit  = m_units[i];
            const loint_type order = lo_common::get_basic_order(unit.lnode);
            if (q_order < order) {
                break;
            }
            const size_type second_order = unit.lnode->get_data().get_lcursor(unit.lrnid).get_lrpos();
            if (q_order == order and q_second_order < second_order) {
                break;
            }
        }
        ABORT_IF_EQ(i, 0);
        return {m_units[i - 1].frnid, i - 1};
    }

    /**
     *  Getter & Setter for each run unit
     */

    // Retrun the character of the group.
    uchar_type get_chr() const {
        DEBUG_ABORT_IF(is_empty());
        const auto& unit = m_units.front();
        return unit.lnode->get_data().get_chr(unit.lnode->get_data().get_lcursor(unit.lrnid));
    }
    // Will return the same character for all runs
    uchar_type get_chr(size_type frnid) const {
        return get_chr(get_fcursor(frnid));
    }
    // Retrun the exponent of the run-unit with 'frnid'.
    size_type get_exp(size_type frnid) const {
        return get_exp(get_fcursor(frnid));
    }
    void add_exp(size_type frnid, size_type exp) {
        add_exp(get_fcursor(frnid), exp);
    }
    // Return the pointer to L's run-unit corresponding to F's run-unit with ID 'frnid'.
    std::pair<lnode_type*, size_type> get_lrptr(size_type frnid) const {
        return get_lrptr(get_fcursor(frnid));
    }
    void set_lrptr(size_type frnid, lnode_type* new_lnode, size_type new_lrnid) {
        set_lrptr(get_fcursor(frnid), new_lnode, new_lrnid);
    }

    uchar_type get_chr(const fcursor_type fcrsr) const {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        return m_units[fcrsr.get_frpos()].lnode->get_data().get_chr(m_units[fcrsr.get_frpos()].lrnid);
    }
    // Retrun the exponent of the run-unit with 'frnid'.
    size_type get_exp(const fcursor_type fcrsr) const {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        return m_units[fcrsr.get_frpos()].lnode->get_data().get_exp(m_units[fcrsr.get_frpos()].lrnid);
    }
    void add_exp(const fcursor_type fcrsr, size_type exp) {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        if constexpr (SUMEXP_HINT) {
            m_sum_exps += exp;
        }
    }
    std::pair<lnode_type*, size_type> get_lrptr(const fcursor_type fcrsr) const {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        return {m_units[fcrsr.get_frpos()].lnode, m_units[fcrsr.get_frpos()].lrnid};
    }
    void set_lrptr(const fcursor_type fcrsr, lnode_type* new_lnode, size_type new_lrnid) {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        m_units[fcrsr.get_frpos()].lnode = new_lnode;
        m_units[fcrsr.get_frpos()].lrnid = new_lrnid;
    }

    /**
     *  Updater
     */

    // Initillay insert a new F-run corresponding to L-run (new_lnode, new_lrnid).
    // Return the run-ID
    fcursor_type insert_init(lnode_type* new_lnode, size_type new_lrnid) {
        ABORT_IF(!is_empty());

        const size_type new_frnid = make_free_frnid();
        const size_type new_frpos = m_units.size();

        m_units.push_back(unit_type{new_lnode, new_lrnid, new_frnid});

        if constexpr (SUMEXP_HINT) {
            m_sum_exps += new_lnode->get_data().get_exp(new_lrnid);
        }

        if constexpr (LOOKUP_HINT) {
            std::fill_n(m_id_to_pos.data(), MAX_NUM_RUNS, MAX_SIZE_INT);
            m_id_to_pos[new_frnid] = new_frpos;
        }

        return {new_frnid, new_frpos};  // new run-ID
    }
    // Insert a new F-run corresponding to L-run (new_lnode, new_lrnid) before the F-run with frnid.
    // Return the run-ID
    fcursor_type insert_before(fcursor_type& fcrsr, lnode_type* new_lnode, size_type new_lrnid) {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));

        const size_type new_frnid = make_free_frnid();
        const size_type new_frpos = fcrsr.get_frpos();

        m_units.insert(m_units.begin() + new_frpos, unit_type{new_lnode, new_lrnid, new_frnid});
        fcrsr.m_frpos += 1;

        if constexpr (SUMEXP_HINT) {
            m_sum_exps += new_lnode->get_data().get_exp(new_lrnid);
        }

        if constexpr (LOOKUP_HINT) {
            DEBUG_ABORT_IF_NE(m_id_to_pos[new_frnid], MAX_SIZE_INT);
            m_id_to_pos[new_frnid] = new_frpos;
            for (size_type i = new_frpos + 1; i < m_units.size(); i++) {
                DEBUG_ABORT_IF_NE(m_id_to_pos[m_units[i].frnid], i - 1);
                m_id_to_pos[m_units[i].frnid] = i;
            }
        }

        return {new_frnid, new_frpos};
    }
    // Insert a new F-run corresponding to L-run (new_lnode, new_lrnid) after the F-run with frnid.
    // Return the run-ID
    fcursor_type insert_after(fcursor_type& fcrsr, lnode_type* new_lnode, size_type new_lrnid) {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        // insert F-node
        const size_type new_frnid = make_free_frnid();
        const size_type new_frpos = fcrsr.get_frpos() + 1;

        m_units.insert(m_units.begin() + new_frpos, unit_type{new_lnode, new_lrnid, new_frnid});

        if constexpr (SUMEXP_HINT) {
            m_sum_exps += new_lnode->get_data().get_exp(new_lrnid);
        }

        if constexpr (LOOKUP_HINT) {
            DEBUG_ABORT_IF_NE(m_id_to_pos[new_frnid], MAX_SIZE_INT);
            m_id_to_pos[new_frnid] = new_frpos;
            for (size_type i = new_frpos + 1; i < m_units.size(); i++) {
                DEBUG_ABORT_IF_NE(m_id_to_pos[m_units[i].frnid], i - 1);
                m_id_to_pos[m_units[i].frnid] = i;
            }
        }

        // const size_type new_frnid = make_free_frnid();
        // m_units.insert(m_units.begin() + fcrsr.get_frpos() + 1, unit_type{new_lnode, new_lrnid, new_frnid});
        return {new_frnid, new_frpos};
    }
    // やってることはinsert_afterと同じだけど、sum_of_expsを陽に持つ場合にその能力は発揮されるであろう（ごごご
    fcursor_type split_after(fcursor_type& fcrsr, lnode_type* new_lnode, size_type new_lrnid) {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        // insert F-node
        const size_type new_frnid = make_free_frnid();
        const size_type new_frpos = fcrsr.get_frpos() + 1;

        m_units.insert(m_units.begin() + new_frpos, unit_type{new_lnode, new_lrnid, new_frnid});

        // splitのみなので、m_sum_expsは増えない

        if constexpr (LOOKUP_HINT) {
            DEBUG_ABORT_IF_NE(m_id_to_pos[new_frnid], MAX_SIZE_INT);
            m_id_to_pos[new_frnid] = new_frpos;
            for (size_type i = new_frpos + 1; i < m_units.size(); i++) {
                DEBUG_ABORT_IF_NE(m_id_to_pos[m_units[i].frnid], i - 1);
                m_id_to_pos[m_units[i].frnid] = i;
            }
        }

        // const size_type new_frnid = make_free_frnid();
        // m_units.insert(m_units.begin() + fcrsr.get_frpos() + 1, unit_type{new_lnode, new_lrnid, new_frnid});
        return {new_frnid, new_frpos};
    }

    /**
     *  Offset handlers
     */

    // Return the offset of run frnid.
    offset_type get_offset(size_type frnid) const {
        return get_offset(get_fcursor(frnid));
    }
    offset_type get_offset(const fcursor_type fcrsr) const {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        offset_type fofst = 0;
        for (size_type i = 0; i < fcrsr.get_frpos(); i++) {
            fofst += m_units[i].lnode->get_data().get_exp(m_units[i].lrnid);
        }
        return fofst;
    }
    std::tuple<fcursor_type, offset_type> access_with_offset(offset_type fofst) {
        const size_type num_runs = get_num_runs();
        for (size_type frpos = 0; frpos < num_runs; frpos++) {
            const size_type exp = m_units[frpos].lnode->get_data().get_exp(m_units[frpos].lrnid);
            if (fofst < offset_type(exp)) {
                return {fcursor_type(m_units[frpos].frnid, frpos), fofst};
            }
            fofst = fofst - exp;
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
        const size_type num_runs_aftr = get_num_runs() - num_runs_befo;

        // Clone the rear part (B)
        this_type data_aftr;
        std::copy(m_units.begin() + num_runs_befo, m_units.end(), std::back_inserter(data_aftr.m_units));

        // Update A
        m_units.resize(num_runs_befo);
        for (size_type i = 0; i < num_runs_befo; i++) {
            m_units[i].frnid = i;
        }
        // Update B
        for (size_type i = 0; i < num_runs_aftr; i++) {
            data_aftr.m_units[i].frnid = i;
        }

        if constexpr (SUMEXP_HINT) {
            m_sum_exps = 0;
            for (size_type i = 0; i < num_runs_befo; i++) {
                const unit_type& unit = m_units[i];
                m_sum_exps += unit.lnode->get_data().get_exp(unit.lrnid);
            }
            data_aftr.m_sum_exps = 0;
            for (size_type i = 0; i < num_runs_aftr; i++) {
                const unit_type& unit = data_aftr.m_units[i];
                data_aftr.m_sum_exps += unit.lnode->get_data().get_exp(unit.lrnid);
            }
        }

        if constexpr (LOOKUP_HINT) {
            std::fill_n(m_id_to_pos.data(), MAX_NUM_RUNS, MAX_SIZE_INT);
            for (size_type i = 0; i < num_runs_befo; i++) {
                m_id_to_pos[m_units[i].frnid] = i;
            }
            std::fill_n(data_aftr.m_id_to_pos.data(), MAX_NUM_RUNS, MAX_SIZE_INT);
            for (size_type i = 0; i < num_runs_aftr; i++) {
                data_aftr.m_id_to_pos[data_aftr.m_units[i].frnid] = i;
            }
        }

        return data_aftr;
    }

  private:
    size_type make_free_frnid() {
        if (m_units.empty()) {
            return 0;
        }
        size_type max_frnid = 0;
        for (const auto& unit : m_units) {
            max_frnid = std::max(max_frnid, unit.frnid);
        }
        return max_frnid + 1;
    }

    /**
     *  ID <-> Pos
     */
    // Retrun the position of the run-unit with 'frnid'.
    size_type get_frpos(size_type frnid) const {
        if constexpr (LOOKUP_HINT) {
            DEBUG_ABORT_IF_EQ(m_id_to_pos[frnid], MAX_SIZE_INT);
            return m_id_to_pos[frnid];
        } else {
            for (size_type frpos = 0; frpos < m_units.size(); frpos++) {
                if (m_units[frpos].frnid == frnid) {
                    return frpos;
                }
            }
            return MAX_SIZE_INT;
        }
    }
    // Returns the ID of the run-unit located at 'lrpos'
    size_type get_frnid(offset_type frpos) const {
        return m_units[frpos].frnid;
    }
};

}  // namespace rcomp
