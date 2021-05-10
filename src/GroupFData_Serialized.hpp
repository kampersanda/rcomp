/**
 * @file GroupFData_Serialized.hpp
 * @author Shunsuke Kanda (kampersanda)
 */
#pragma once

#include <memory>
#include <vector>

#include "LinkedList.hpp"
#include "intrinsics.hpp"
#include "utils.hpp"

namespace rcomp {

/**
 * A class for space-efficient data type of grouped F-node.
 *
 * @tparam t_LData The data type of grouped L-node.
 */
template <class t_LData>
class GroupFData_Serialized {
  public:
    using this_type  = GroupFData_Serialized<t_LData>;
    using lnode_type = typename LinkedList<t_LData>::node_type;

    static constexpr auto GROUP_BOUND = t_LData::GROUP_BOUND;
    static constexpr auto SUMEXP_HINT = t_LData::SUMEXP_HINT;
    static constexpr auto LOOKUP_HINT = t_LData::LOOKUP_HINT;
    static constexpr auto GROW_SIZE   = t_LData::GROW_SIZE;

    static constexpr auto MAX_NUM_RUNS = GROUP_BOUND + 2;
    static_assert(MAX_NUM_RUNS <= 256);

    // Cursor
    class fcursor_type {
        friend class GroupFData_Serialized;

        uint8_t m_frnid = 0;
        uint8_t m_frpos = 0;

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
    uint8_t m_exp_bytes = 0;
    uint8_t m_size_runs = 0;
    uint8_t m_capa_runs = 0;

    uchar_type  m_chr    = '\0';
    lnode_type* m_tnode  = nullptr;
    uint8_t     m_weight = 0;

    std::unique_ptr<uint8_t[]> m_units;

    size_type                         m_sum_exps = 0;  // Compressible, but would not be reasonable
    std::array<uint8_t, MAX_NUM_RUNS> m_id_to_pos;  // Should be embeded in m_units?

    struct unit_type {
        lnode_type* lnode;
        size_type   lrnid;
        size_type   frnid;
        size_type   exp;
    };

  public:
    GroupFData_Serialized() = default;

    //! Default destructor
    virtual ~GroupFData_Serialized() = default;

    //! Copy constructor (deleted)
    GroupFData_Serialized(const GroupFData_Serialized&) = delete;

    //! Copy constructor (deleted)
    GroupFData_Serialized& operator=(const GroupFData_Serialized&) = delete;

    //! Move constructor
    GroupFData_Serialized(GroupFData_Serialized&&) noexcept = default;

    //! Move constructor
    GroupFData_Serialized& operator=(GroupFData_Serialized&&) noexcept = default;

    inline void prefetch() const {
        intrinsics::prefetch(m_units.get());
    }

    inline bool is_empty() const {
        return !m_units.get();
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
        return sizeof(m_exp_bytes) + sizeof(m_size_runs) + sizeof(m_capa_runs);
    }
    size_type get_link_info_memory_in_bytes() const {
        return sizeof(m_chr) + sizeof(m_tnode) + sizeof(m_weight);
    }
    size_type get_unit_info_memory_in_bytes() const {
        return get_capa_runs() * get_unit_bytes(m_exp_bytes) + sizeof(m_units);
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

    inline size_type get_num_runs() const {
        return m_size_runs;
    }
    inline size_type get_capa_runs() const {
        return m_capa_runs;
    }

    // Retrun the sum of exponents in the group
    inline size_type get_sum_exps() const {
        if constexpr (SUMEXP_HINT) {
            return m_sum_exps;
        } else {
            size_type sum_exps = 0;
            for (size_type frpos = 0; frpos < get_num_runs(); frpos++) {
                sum_exps += get_exp_from_position(frpos);
            }
            return sum_exps;
        }
    }

    /**
     *  Accessor by L-run ID
     */

    inline fcursor_type get_fcursor(size_type frnid) const {
        return {frnid, get_frpos_from_identifier(frnid)};
    }
    inline fcursor_type get_fcursor_from_position(size_type frpos) const {
        return {get_frnid_from_position(frpos), frpos};
    }

    bool is_valid_fcursor(const fcursor_type fcrsr) const {
        DEBUG_ABORT_IF(is_empty());
        return fcrsr.get_frpos() == get_frpos_from_identifier(fcrsr.get_frnid());
    }

    inline fcursor_type get_first_fcursor() const {
        DEBUG_ABORT_IF(is_empty());
        return {get_frnid_from_position(0), 0};
    }
    inline fcursor_type get_last_fcursor() const {
        DEBUG_ABORT_IF(is_empty());
        return {get_frnid_from_position(get_num_runs() - 1), get_num_runs() - 1};
    }
    inline fcursor_type get_prev_fcursor(const fcursor_type fcrsr) const {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        DEBUG_ABORT_IF_LE(fcrsr.get_frpos(), 0);
        return {get_frnid_from_position(fcrsr.get_frpos() - 1), fcrsr.get_frpos() - 1};
    }
    inline fcursor_type get_next_fcursor(const fcursor_type fcrsr) const {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        DEBUG_ABORT_IF_LE(get_num_runs() - 1, fcrsr.get_frpos());
        return {get_frnid_from_position(fcrsr.get_frpos() + 1), fcrsr.get_frpos() + 1};
    }

    inline bool is_first_fcursor(const fcursor_type fcrsr) const {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        return fcrsr.get_frpos() == 0;
    }
    inline bool is_last_fcursor(const fcursor_type fcrsr) const {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        return fcrsr.get_frpos() == get_num_runs() - 1;
    }

    /**
     *  Getter & Setter for the group member
     */
    inline lnode_type* get_tnode() const {
        return m_tnode;
    }
    inline offset_type get_tofst() const {
        return m_tnode->get_data().get_hofst();
    }
    inline void set_tnode(lnode_type* tnode) {
        m_tnode = tnode;
    }
    inline void reset_tlink() {
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
    inline loint_type get_order() const {
        DEBUG_ABORT_IF(is_empty());
        return get_lnode_from_position(0)->get_data().get_order();
    }
    // Retrun the run-unit position in the L-node for the first run-unit of this group.
    inline size_type get_second_order() const {  // second order?
        DEBUG_ABORT_IF(is_empty());
        return get_lnode_from_position(0)->get_data().get_lcursor(get_lrnid_from_position(0)).get_lrpos();
    }
    // Search the most backword run-unit whose order <= q_order and second_order <= q_second_order.
    // Return its run-ID.
    inline fcursor_type predecessor(const loint_type q_order, const size_type q_second_order) const {
        DEBUG_ABORT_IF(is_empty());

        size_type frpos = 0;
        for (; frpos < get_num_runs(); frpos++) {
            const lnode_type* lnode = get_lnode_from_position(frpos);
            const size_type   lrnid = get_lrnid_from_position(frpos);

            const loint_type order = lo_common::get_basic_order(lnode);
            if (q_order < order) {
                break;
            }

            const size_type second_order = lnode->get_data().get_lcursor(lrnid).get_lrpos();
            if (q_order == order and q_second_order < second_order) {
                break;
            }
        }

        ABORT_IF_EQ(frpos, 0);
        return {get_frnid_from_position(frpos - 1), frpos - 1};
    }

    /**
     *  Getter & Setter for each run unit
     */

    // Retrun the character of the group.
    inline uchar_type get_chr() const {
        DEBUG_ABORT_IF(is_empty());
        DEBUG_ABORT_IF_NE(m_chr, get_chr_from_position(0));
        return m_chr;
    }
    // Will return the same character for all runs
    inline uchar_type get_chr(size_type frnid) const {
        return get_chr(get_fcursor(frnid));
    }
    // Retrun the exponent of the run-unit with 'frnid'.
    inline size_type get_exp(size_type frnid) const {
        return get_exp(get_fcursor(frnid));
    }
    inline void add_exp(size_type frnid, size_type v) {
        add_exp(get_fcursor(frnid), v);
    }
    // Return the pointer to L's run-unit corresponding to F's run-unit with ID 'frnid'.
    inline std::pair<lnode_type*, size_type> get_lrptr(size_type frnid) const {
        return get_lrptr(get_fcursor(frnid));
    }
    inline void set_lrptr(size_type frnid, lnode_type* new_lnode, size_type new_lrnid) {
        set_lrptr(get_fcursor(frnid), new_lnode, new_lrnid);
    }

    inline uchar_type get_chr(const fcursor_type fcrsr) const {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        DEBUG_ABORT_IF_NE(m_chr, get_chr_from_position(fcrsr.get_frpos()));
        return m_chr;
    }
    // Retrun the exponent of the run-unit with 'frnid'.
    inline size_type get_exp(const fcursor_type fcrsr) const {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        return get_exp_from_position(fcrsr.get_frpos());
    }
    inline void add_exp(const fcursor_type fcrsr, size_type v) {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));

        const size_type new_exp = get_exp(fcrsr) + v;

        try_to_realloc(new_exp);
        set_exp_in_position(fcrsr.get_frpos(), new_exp);

        if constexpr (SUMEXP_HINT) {
            m_sum_exps += v;
        }
    }
    inline std::pair<lnode_type*, size_type> get_lrptr(const fcursor_type fcrsr) const {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        return {get_lnode_from_position(fcrsr.get_frpos()), get_lrnid_from_position(fcrsr.get_frpos())};
    }
    inline void set_lrptr(const fcursor_type fcrsr, lnode_type* new_lnode, size_type new_lrnid) {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        set_lnode_in_position(fcrsr.get_frpos(), new_lnode);
        set_lrnid_in_position(fcrsr.get_frpos(), new_lrnid);
    }

    /**
     *  Updater
     */

    // Initillay insert a new F-run corresponding to L-run (new_lnode, new_lrnid).
    // Return the run-ID
    inline fcursor_type insert_init(lnode_type* new_lnode, size_type new_lrnid) {
        ABORT_IF(!is_empty());

        const size_type new_frnid = get_free_lrnid();
        const size_type new_frpos = m_size_runs;

        const size_type new_exp  = new_lnode->get_data().get_exp(new_lrnid);
        const unit_type new_unit = {new_lnode, new_lrnid, new_frnid, new_exp};

        realloc_units(utils::get_nbytes(new_exp), 1, GROW_SIZE);
        set_unit_in_position(new_frpos, new_unit);

        // Set the representative character of this group
        m_chr = new_lnode->get_data().get_chr(new_lrnid);

        if constexpr (SUMEXP_HINT) {
            m_sum_exps += new_exp;
        }
        if constexpr (LOOKUP_HINT) {
            std::fill_n(m_id_to_pos.data(), MAX_NUM_RUNS, UINT8_MAX);
            m_id_to_pos[new_frnid] = uint8_t(new_frpos);
        }

        return {new_frnid, new_frpos};  // new run-ID
    }
    // Insert a new F-run corresponding to L-run (new_lnode, new_lrnid) before the F-run with frnid.
    // Return the run-ID
    inline fcursor_type insert_before(fcursor_type& fcrsr, lnode_type* new_lnode, size_type new_lrnid) {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));

        const size_type new_frnid = get_free_lrnid();
        const size_type new_frpos = fcrsr.get_frpos();

        const size_type new_exp  = new_lnode->get_data().get_exp(new_lrnid);
        const unit_type new_unit = {new_lnode, new_lrnid, new_frnid, new_exp};

        try_to_realloc(new_exp);
        make_free_unit(new_frpos);
        set_unit_in_position(new_frpos, new_unit);

        fcrsr.m_frpos += 1;

        if constexpr (SUMEXP_HINT) {
            m_sum_exps += new_exp;
        }

        if constexpr (LOOKUP_HINT) {
            DEBUG_ABORT_IF_NE(m_id_to_pos[new_frnid], UINT8_MAX);
            m_id_to_pos[new_frnid] = uint8_t(new_frpos);
            for (size_type i = new_frpos + 1; i < get_num_runs(); i++) {
                DEBUG_ABORT_IF_NE(m_id_to_pos[get_frnid_from_position(i)], i - 1);
                m_id_to_pos[get_frnid_from_position(i)] = uint8_t(i);
            }
        }

        return {new_frnid, new_frpos};
    }
    // Insert a new F-run corresponding to L-run (new_lnode, new_lrnid) after the F-run with frnid.
    // Return the run-ID
    inline fcursor_type insert_after(fcursor_type& fcrsr, lnode_type* new_lnode, size_type new_lrnid) {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));

        const size_type new_frnid = get_free_lrnid();
        const size_type new_frpos = fcrsr.get_frpos() + 1;

        const size_type new_exp  = new_lnode->get_data().get_exp(new_lrnid);
        const unit_type new_unit = {new_lnode, new_lrnid, new_frnid, new_exp};

        try_to_realloc(new_exp);
        make_free_unit(new_frpos);
        set_unit_in_position(new_frpos, new_unit);

        if constexpr (SUMEXP_HINT) {
            m_sum_exps += new_exp;
        }

        if constexpr (LOOKUP_HINT) {
            DEBUG_ABORT_IF_NE(m_id_to_pos[new_frnid], UINT8_MAX);
            m_id_to_pos[new_frnid] = uint8_t(new_frpos);
            for (size_type i = new_frpos + 1; i < get_num_runs(); i++) {
                DEBUG_ABORT_IF_NE(m_id_to_pos[get_frnid_from_position(i)], i - 1);
                m_id_to_pos[get_frnid_from_position(i)] = uint8_t(i);
            }
        }

        return {new_frnid, new_frpos};
    }

    inline fcursor_type split_after(fcursor_type& fcrsr, lnode_type* new_lnode, size_type new_lrnid) {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));

        const size_type new_frnid = get_free_lrnid();
        const size_type new_frpos = fcrsr.get_frpos() + 1;

        const size_type new_exp_b = new_lnode->get_data().get_exp(new_lrnid);  // after
        const size_type new_exp_a = get_exp_from_position(fcrsr.get_frpos()) - new_exp_b;  // before

        const unit_type new_unit_b = {new_lnode, new_lrnid, new_frnid, new_exp_b};  // after

        set_exp_in_position(fcrsr.get_frpos(), new_exp_a);

        make_free_unit(new_frpos);
        set_unit_in_position(new_frpos, new_unit_b);

        if constexpr (LOOKUP_HINT) {
            DEBUG_ABORT_IF_NE(m_id_to_pos[new_frnid], UINT8_MAX);
            m_id_to_pos[new_frnid] = uint8_t(new_frpos);
            for (size_type i = new_frpos + 1; i < get_num_runs(); i++) {
                DEBUG_ABORT_IF_NE(m_id_to_pos[get_frnid_from_position(i)], i - 1);
                m_id_to_pos[get_frnid_from_position(i)] = uint8_t(i);
            }
        }

        return {new_frnid, new_frpos};
    }

    /**
     *  Offset handlers
     */

    // Return the offset of run frnid.
    inline offset_type get_offset(size_type frnid) const {
        return get_offset(get_fcursor(frnid));
    }
    inline offset_type get_offset(const fcursor_type fcrsr) const {
        DEBUG_ABORT_IF(!is_valid_fcursor(fcrsr));
        offset_type fofst = 0;
        for (size_type frpos = 0; frpos < fcrsr.get_frpos(); frpos++) {
            fofst += get_exp_from_position(frpos);
        }
        return fofst;
    }
    inline std::tuple<fcursor_type, offset_type> access_with_offset(offset_type fofst) {
        const size_type num_runs = get_num_runs();
        for (size_type frpos = 0; frpos < num_runs; frpos++) {
            const offset_type exp = get_exp_from_position(frpos);
            if (fofst < exp) {
                return {get_fcursor_from_position(frpos), fofst};
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
    inline this_type divide_after() {
        ABORT_IF_LE(get_num_runs(), 1);
        return divide_after(get_num_runs() / 2);
    }
    inline this_type divide_after(const size_type num_runs_befo) {
        ABORT_IF_LE(get_num_runs(), num_runs_befo);
        const size_type num_runs_aftr = get_num_runs() - num_runs_befo;

        // Clone the rear part
        this_type data_aftr;
        {
            const uint8_t* beg_aftr = get_unit_pointer(num_runs_befo);
            const uint8_t* end_aftr = get_unit_pointer(get_num_runs());

            data_aftr.m_exp_bytes = m_exp_bytes;
            data_aftr.m_size_runs = num_runs_aftr;
            data_aftr.m_capa_runs = num_runs_aftr;

            data_aftr.m_chr = m_chr;

            data_aftr.m_units = std::make_unique<uint8_t[]>(std::distance(beg_aftr, end_aftr));
            std::copy(beg_aftr, end_aftr, data_aftr.m_units.get());
        }

        // Resize the front part
        m_size_runs = num_runs_befo;
        m_capa_runs = num_runs_befo;
        {
            const auto new_bytes = get_unit_begin(num_runs_befo, m_exp_bytes);

            auto new_serialized = std::make_unique<uint8_t[]>(new_bytes);
            std::copy_n(m_units.get(), new_bytes, new_serialized.get());
            m_units = std::move(new_serialized);
        }

        for (size_type i = 0; i < num_runs_befo; i++) {
            set_frnid_in_position(i, i);
        }
        // Update B
        for (size_type i = 0; i < num_runs_aftr; i++) {
            data_aftr.set_frnid_in_position(i, i);
        }

        if constexpr (SUMEXP_HINT) {
            m_sum_exps = 0;
            for (size_type i = 0; i < num_runs_befo; i++) {
                m_sum_exps += get_exp_from_position(i);
            }
            data_aftr.m_sum_exps = 0;
            for (size_type i = 0; i < num_runs_aftr; i++) {
                data_aftr.m_sum_exps += data_aftr.get_exp_from_position(i);
            }
        }

        if constexpr (LOOKUP_HINT) {
            std::fill_n(m_id_to_pos.data(), MAX_NUM_RUNS, UINT8_MAX);
            for (size_type i = 0; i < num_runs_befo; i++) {
                m_id_to_pos[get_frnid_from_position(i)] = i;
            }
            std::fill_n(data_aftr.m_id_to_pos.data(), MAX_NUM_RUNS, UINT8_MAX);
            for (size_type i = 0; i < num_runs_aftr; i++) {
                data_aftr.m_id_to_pos[data_aftr.get_frnid_from_position(i)] = i;
            }
        }

        return data_aftr;
    }

  private:
    inline size_type get_free_lrnid() const {
        return m_size_runs;
    }

    inline void realloc_units(const uint8_t new_exp_bytes, const uint8_t new_size_runs, const uint8_t new_capa_runs) {
        auto resize_runs    = std::min(m_size_runs, new_size_runs);
        auto new_serialized = std::make_unique<uint8_t[]>(get_unit_bytes(new_exp_bytes) * new_capa_runs);

        for (size_type frpos = 0; frpos < resize_runs; frpos++) {
            uint8_t* new_ptr = new_serialized.get() + get_unit_begin(frpos, new_exp_bytes);
            set_unit_in_pointer(new_ptr, get_unit_from_position(frpos), new_exp_bytes);
        }

        m_units     = std::move(new_serialized);
        m_exp_bytes = new_exp_bytes;
        m_size_runs = new_size_runs;
        m_capa_runs = new_capa_runs;
    }

    inline void try_to_realloc(const size_type new_exp) {
        const uint8_t new_exp_bytes = utils::get_nbytes(new_exp);
        if (m_exp_bytes < new_exp_bytes) {
            realloc_units(new_exp_bytes, m_size_runs, m_capa_runs);
            DEBUG_ABORT_IF_NE(m_exp_bytes, new_exp_bytes);
        }
    }

    inline void make_free_unit(size_type frpos) {
        if (m_size_runs == m_capa_runs) {
            realloc_units(m_exp_bytes, m_size_runs, std::min(m_capa_runs + GROW_SIZE, MAX_NUM_RUNS));
        }

        uint8_t*  dst   = m_units.get() + get_unit_begin(frpos + 1, m_exp_bytes);
        uint8_t*  src   = m_units.get() + get_unit_begin(frpos, m_exp_bytes);
        size_type bytes = get_unit_bytes(m_exp_bytes) * (m_size_runs - frpos);

        std::memmove(dst, src, bytes);

        m_size_runs += 1;
    }

    inline uint8_t* get_unit_pointer(size_type frpos) {
        return m_units.get() + get_unit_begin(frpos, m_exp_bytes);
    }
    inline const uint8_t* get_unit_pointer(size_type frpos) const {
        return m_units.get() + get_unit_begin(frpos, m_exp_bytes);
    }

    inline size_type get_frpos_from_identifier(size_type frnid) const {
        if constexpr (LOOKUP_HINT) {
            DEBUG_ABORT_IF_EQ(m_id_to_pos[frnid], UINT8_MAX);
            return m_id_to_pos[frnid];
        } else {
            for (size_type frpos = 0; frpos < get_num_runs(); frpos++) {
                if (get_frnid_from_position(frpos) == frnid) {
                    return frpos;
                }
            }
            return MAX_SIZE_INT;
        }
    }

    /**
     *  Accessors for units with position
     */
    inline uchar_type get_chr_from_position(size_type frpos) const {
        const lnode_type* lnode = get_lnode_from_position(frpos);
        const size_type   lrnid = get_lrnid_from_position(frpos);
        return lnode->get_data().get_chr(lrnid);
    }
    inline lnode_type* get_lnode_from_position(size_type frpos) const {
        return get_lnode_from_pointer(get_unit_pointer(frpos), m_exp_bytes);
    }
    inline size_type get_lrnid_from_position(size_type frpos) const {
        return get_lrnid_from_pointer(get_unit_pointer(frpos), m_exp_bytes);
    }
    inline size_type get_frnid_from_position(size_type frpos) const {
        return get_frnid_from_pointer(get_unit_pointer(frpos), m_exp_bytes);
    }
    inline size_type get_exp_from_position(size_type frpos) const {
        return get_exp_from_pointer(get_unit_pointer(frpos), m_exp_bytes);
    }
    inline void set_lnode_in_position(size_type frpos, lnode_type* lnode) {
        set_lnode_in_pointer(get_unit_pointer(frpos), lnode, m_exp_bytes);
    }
    inline void set_lrnid_in_position(size_type frpos, size_type lrnid) {
        set_lrnid_in_pointer(get_unit_pointer(frpos), lrnid, m_exp_bytes);
    }
    inline void set_frnid_in_position(size_type frpos, size_type frnid) {
        set_frnid_in_pointer(get_unit_pointer(frpos), frnid, m_exp_bytes);
    }
    inline void set_exp_in_position(size_type frpos, size_type exp) {
        set_exp_in_pointer(get_unit_pointer(frpos), exp, m_exp_bytes);
    }

    inline unit_type get_unit_from_position(size_type frpos) {
        return get_unit_from_pointer(get_unit_pointer(frpos), m_exp_bytes);
    }
    inline void set_unit_in_position(size_type frpos, const unit_type& unit) {
        set_unit_in_pointer(get_unit_pointer(frpos), unit, m_exp_bytes);
    }

    /**
     *   Static members for handling serialized units
     */

    static inline size_type get_unit_begin(size_type frpos, uint8_t exp_bytes) {
        return get_unit_bytes(exp_bytes) * frpos;
    }
    static inline size_type get_unit_bytes(uint8_t exp_bytes) {
        // lnode + lrnid + frnid + exp
        return sizeof(uintptr_t) + sizeof(uint8_t) + sizeof(uint8_t) + exp_bytes;
    }

    static inline lnode_type* get_lnode_from_pointer(const uint8_t* ptr, uint8_t) {
        return reinterpret_cast<lnode_type*>(utils::get_memory_block<uintptr_t>(ptr));
    }
    static inline size_type get_lrnid_from_pointer(const uint8_t* ptr, uint8_t) {
        const size_type beg = sizeof(uintptr_t);
        return static_cast<size_type>(utils::get_memory_block<uint8_t>(ptr + beg));
    }
    static inline size_type get_frnid_from_pointer(const uint8_t* ptr, uint8_t) {
        const size_type beg = sizeof(uintptr_t) + sizeof(uint8_t);
        return static_cast<size_type>(utils::get_memory_block<uint8_t>(ptr + beg));
    }
    static inline size_type get_exp_from_pointer(const uint8_t* ptr, uint8_t exp_bytes) {
        const size_type beg = sizeof(uintptr_t) + sizeof(uint8_t) + sizeof(uint8_t);
        return utils::get_memory_block<size_type>(ptr + beg, exp_bytes);
    }

    static inline void set_lnode_in_pointer(uint8_t* ptr, lnode_type* lnode, uint8_t) {
        utils::set_memory_block(ptr, reinterpret_cast<uintptr_t>(lnode));
    }
    static inline void set_lrnid_in_pointer(uint8_t* ptr, size_type lrnid, uint8_t) {
        const size_type beg = sizeof(uintptr_t);
        utils::set_memory_block(ptr + beg, static_cast<uint8_t>(lrnid));
    }
    static inline void set_frnid_in_pointer(uint8_t* ptr, size_type frnid, uint8_t) {
        const size_type beg = sizeof(uintptr_t) + sizeof(uint8_t);
        utils::set_memory_block(ptr + beg, static_cast<uint8_t>(frnid));
    }
    static inline void set_exp_in_pointer(uint8_t* ptr, size_type exp, uint8_t exp_bytes) {
        const size_type beg = sizeof(uintptr_t) + sizeof(uint8_t) + sizeof(uint8_t);
        utils::set_memory_block(ptr + beg, exp, exp_bytes);
    }

    // TODO: Optimization
    static inline unit_type get_unit_from_pointer(uint8_t* ptr, uint8_t exp_bytes) {
        return {get_lnode_from_pointer(ptr, exp_bytes),  //
                get_lrnid_from_pointer(ptr, exp_bytes),  //
                get_frnid_from_pointer(ptr, exp_bytes),  //
                get_exp_from_pointer(ptr, exp_bytes)};
    }
    static inline void set_unit_in_pointer(uint8_t* ptr, const unit_type& unit, uint8_t exp_bytes) {
        set_lnode_in_pointer(ptr, unit.lnode, exp_bytes);
        set_lrnid_in_pointer(ptr, unit.lrnid, exp_bytes);
        set_frnid_in_pointer(ptr, unit.frnid, exp_bytes);
        set_exp_in_pointer(ptr, unit.exp, exp_bytes);
    }
};

}  // namespace rcomp
