/**
 * @file GroupedLData_Serialized.hpp
 */
#pragma once

#include <memory>

#include "GroupedFData_Serialized.hpp"
#include "GroupedLData_UnitTraits.hpp"
#include "LinkedList.hpp"
#include "utils.hpp"

namespace rcomp {

/**
 * A class for a space-efficient grouped L-node.
 *
 * @tparam t_GroupedBound The upper bound of #runs stored.
 * @tparam t_WithSAE Store an SA-entry?
 * @tparam t_GrowSize The constant factor for growing the data.
 * @tparam t_SumExpHint Explicitly keep the sum of exponents?
 * @tparam t_LookupHint Explicitly keep the mapping from run-ids to run-positions?
 */
template <size_type t_GroupedBound, bool t_WithSAE, size_type t_GrowSize = t_GroupedBound / 8,  //
          bool t_SumExpHint = true, bool t_LookupHint = true>
class GroupedLData_Serialized {
    static_assert(sizeof(uchar_type) == 1);

  public:
    using this_type = GroupedLData_Serialized<t_GroupedBound, t_WithSAE, t_GrowSize, t_SumExpHint, t_LookupHint>;
    using fnode_type = typename LinkedList<GroupedFData_Serialized<this_type>>::node_type;
    using unit_type = typename GroupedLData_UnitTraits<t_WithSAE, fnode_type>::unit_type;

    static constexpr auto GROUP_BOUND = t_GroupedBound;
    static constexpr bool WITH_SAE = t_WithSAE;
    static constexpr auto GROW_SIZE = std::max<size_type>(1, t_GrowSize);
    static constexpr auto SUMEXP_HINT = t_SumExpHint;
    static constexpr auto LOOKUP_HINT = t_LookupHint;

    static constexpr auto MAX_NUM_RUNS = GROUP_BOUND + 2;  // +2 is for temporally split nodes
    static_assert(MAX_NUM_RUNS <= 256);

    class lcursor_type {
        friend class GroupedLData_Serialized;

        uint8_t m_lrnid = 0;
        uint8_t m_lrpos = 0;

      public:
        lcursor_type() = default;

        inline size_type get_lrnid() const {
            return m_lrnid;
        }
        inline size_type get_lrpos() const {
            return m_lrpos;
        }

        inline bool operator==(const lcursor_type& other) const {
            return m_lrnid == other.m_lrnid;
        }
        inline bool operator!=(const lcursor_type& other) const {
            return m_lrnid != other.m_lrnid;
        }

      private:
        lcursor_type(size_type lrnid, size_type lrpos) : m_lrnid(lrnid), m_lrpos(lrpos) {}
    };

  private:
    // 8-bit members
    uint8_t m_exp_bytes = 0;
    uint8_t m_size_runs = 0;
    uint8_t m_capa_runs = 0;
    uint8_t m_weight = 0;
    std::array<uint8_t, MAX_NUM_RUNS> m_id_to_pos;
    // 64-bit members
    loint_type m_order = 0;
    fnode_type* m_hnode = nullptr;
    offset_type m_hofst = 0;
    size_type m_sum_exps = 0;
    std::unique_ptr<uint8_t[]> m_units;

  public:
    //! Default constructor
    GroupedLData_Serialized() = default;

    //! Default destructor
    virtual ~GroupedLData_Serialized() = default;

    //! Copy constructor (deleted)
    GroupedLData_Serialized(const GroupedLData_Serialized&) = delete;

    //! Copy constructor (deleted)
    GroupedLData_Serialized& operator=(const GroupedLData_Serialized&) = delete;

    //! Move constructor
    GroupedLData_Serialized(GroupedLData_Serialized&&) noexcept = default;

    //! Move constructor
    GroupedLData_Serialized& operator=(GroupedLData_Serialized&&) noexcept = default;

    inline void prefetch() const {
        intrinsics::prefetch(m_units.get());
    }

    bool is_empty() const {
        return !m_units.get();
    }

    //! Get the allocated memory in bytes.
    size_type get_memory_in_bytes(bool include_this = true) const {
        const size_type this_bytes = sizeof(*this) * include_this;
        return this_bytes + get_capa_runs() * get_unit_bytes(m_exp_bytes);
    }

    size_type get_size_info_memory_in_bytes() const {
        return sizeof(m_exp_bytes) + sizeof(m_size_runs) + sizeof(m_capa_runs);
    }
    size_type get_link_info_memory_in_bytes() const {
        return sizeof(m_order) + sizeof(m_hnode) + sizeof(m_hofst) + sizeof(m_weight);
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
            for (size_type lrpos = 0; lrpos < get_num_runs(); lrpos++) {
                sum_exps += get_exp_from_position(lrpos);
            }
            return sum_exps;
        }
    }

    /**
     *  Accessor by L-run ID
     */

    inline lcursor_type get_lcursor(size_type lrnid) const {
        return {lrnid, get_lrpos_from_identifier(lrnid)};
    }
    inline lcursor_type get_lcursor_from_position(size_type lrpos) const {
        return {get_lrnid_from_position(lrpos), lrpos};
    }

    // Check if lcrsr has a valid id and position.
    bool is_valid_lcursor(const lcursor_type lcrsr) const {
        DEBUG_ABORT_IF(is_empty());
        return lcrsr.get_lrpos() == get_lrpos_from_identifier(lcrsr.get_lrnid());
    }

    inline lcursor_type get_first_lcursor() const {
        DEBUG_ABORT_IF(is_empty());
        return {get_lrnid_from_position(0), 0};
    }
    inline lcursor_type get_last_lcursor() const {
        DEBUG_ABORT_IF(is_empty());
        return {get_lrnid_from_position(get_num_runs() - 1), get_num_runs() - 1};
    }
    inline lcursor_type get_prev_lcursor(const lcursor_type lcrsr) const {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        DEBUG_ABORT_IF_LE(lcrsr.get_lrpos(), 0);
        return {get_lrnid_from_position(lcrsr.get_lrpos() - 1), lcrsr.get_lrpos() - 1};
    }
    inline lcursor_type get_next_lcursor(const lcursor_type lcrsr) const {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        DEBUG_ABORT_IF_LE(get_num_runs() - 1, lcrsr.get_lrpos());
        return {get_lrnid_from_position(lcrsr.get_lrpos() + 1), lcrsr.get_lrpos() + 1};
    }
    inline bool is_first_lcursor(const lcursor_type lcrsr) const {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        return lcrsr.get_lrpos() == 0;
    }
    inline bool is_last_lcursor(const lcursor_type lcrsr) const {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        return lcrsr.get_lrpos() == get_num_runs() - 1;
    }

    /**
     *  Getter & Setter for the group member
     */

    // Retrun the order lableled with list-label maintainance
    inline loint_type get_order() const {
        return m_order;
    }
    // Return the head-linked F-node
    inline fnode_type* get_hnode() const {
        return m_hnode;
    }
    // Return the offset to the head-linked F-node
    inline offset_type get_hofst() const {
        return m_hofst;
    }
    inline std::pair<fnode_type*, offset_type> get_hlink() const {
        return {m_hnode, m_hofst};
    }
    inline void set_order(loint_type order) {
        m_order = order;
    }
    inline void set_hnode(fnode_type* hnode) {
        m_hnode = hnode;
    }
    inline void set_hofst(offset_type hofst) {
        m_hofst = hofst;
    }
    inline void set_hlink(fnode_type* hnode, offset_type hofst) {
        m_hnode = hnode;
        m_hofst = hofst;
    }
    inline void reset_hlink() {
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
    inline uchar_type get_chr(size_type lrnid) const {
        return get_chr(get_lcursor(lrnid));
    }
    // Retrun the exponent of the run-unit with ID 'lrnid'.
    inline size_type get_exp(size_type lrnid) const {
        return get_exp(get_lcursor(lrnid));
    }
    //
    inline void add_exp(size_type lrnid, size_type v) {
        add_exp(get_lcursor(lrnid), v);
    }
    // Return the pointer to F's run-unit corresponding to L's run-unit with ID 'lrnid'.
    inline std::pair<fnode_type*, size_type> get_frptr(size_type lrnid) const {
        return get_frptr(get_lcursor(lrnid));
    }
    //
    inline void set_frptr(size_type lrnid, fnode_type* new_fnode, size_type new_frnid) {
        set_frptr(get_lcursor(lrnid), new_fnode, new_frnid);
    }

    // With cursor
    inline uchar_type get_chr(lcursor_type lcrsr) const {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        return get_chr_from_position(lcrsr.get_lrpos());
    }
    inline size_type get_exp(lcursor_type lcrsr) const {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        return get_exp_from_position(lcrsr.get_lrpos());
    }
    // Increment the exponent of the run-unit with ID 'lrnid' by 'exp'.
    inline void add_exp(lcursor_type lcrsr, size_type v) {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));

        const size_type new_exp = get_exp(lcrsr) + v;

        try_to_realloc(new_exp);
        set_exp_in_position(lcrsr.get_lrpos(), new_exp);

        if constexpr (SUMEXP_HINT) {
            m_sum_exps += v;
        }
    }
    inline std::pair<fnode_type*, size_type> get_frptr(lcursor_type lcrsr) const {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        return {get_fnode_from_position(lcrsr.get_lrpos()), get_frnid_from_position(lcrsr.get_lrpos())};
    }
    inline void set_frptr(lcursor_type lcrsr, fnode_type* new_fnode, size_type new_frnid) {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        set_fnode_in_position(lcrsr.get_lrpos(), new_fnode);
        set_frnid_in_position(lcrsr.get_lrpos(), new_frnid);
    }

    inline size_type get_sae(size_type lrnid) const {
        if constexpr (WITH_SAE) {
            return get_sae(get_lcursor(lrnid));
        } else {
            ABORT_IF(true);
        }
    }
    inline size_type get_sae(lcursor_type lcrsr) const {
        if constexpr (WITH_SAE) {
            DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
            return get_sae_from_position(lcrsr.get_lrpos());
        } else {
            ABORT_IF(true);
        }
    }

    inline void set_sae(size_type lrnid, size_type sae) {
        if constexpr (WITH_SAE) {
            set_sae(get_lcursor(lrnid), sae);
        } else {
            ABORT_IF(true);
        }
    }
    inline void set_sae(lcursor_type lcrsr, size_type sae) {
        if constexpr (WITH_SAE) {
            DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
            set_sae_in_position(lcrsr.get_lrpos(), sae);
        } else {
            ABORT_IF(true);
        }
    }

    /**
     *  Updater
     */

    // Initillay insert a new L-run with (new_chr, new_exp)
    // Return the run-ID
    inline lcursor_type insert_init(const uchar_type new_chr, const size_type new_exp) {
        ABORT_IF(!is_empty());

        const size_type new_lrnid = get_free_lrnid();
        const size_type new_lrpos = m_size_runs;
        const unit_type new_unit = {new_chr, new_exp, nullptr, 0, new_lrnid};

        realloc_units(utils::get_nbytes(new_exp), 1, GROW_SIZE);
        set_unit_in_position(new_lrpos, new_unit);

        if constexpr (SUMEXP_HINT) {
            m_sum_exps += new_exp;
        }
        if constexpr (LOOKUP_HINT) {
            std::fill_n(m_id_to_pos.data(), MAX_NUM_RUNS, UINT8_MAX);
            m_id_to_pos[new_lrnid] = uint8_t(new_lrpos);
        }

        return {new_lrnid, new_lrpos};  // new run-ID
    }
    // Insert a new L-run with (new_chr, new_exp) before the L-run with lrnid.
    // Return the run-ID
    inline lcursor_type insert_before(lcursor_type& lcrsr, const uchar_type new_chr, const size_type new_exp) {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));

        const size_type new_lrnid = get_free_lrnid();
        const size_type new_lrpos = lcrsr.get_lrpos();
        const unit_type new_unit = {new_chr, new_exp, nullptr, 0, new_lrnid};

        try_to_realloc(new_exp);
        make_free_unit(new_lrpos);
        set_unit_in_position(new_lrpos, new_unit);

        lcrsr.m_lrpos += 1;

        if constexpr (SUMEXP_HINT) {
            m_sum_exps += new_exp;
        }
        if constexpr (LOOKUP_HINT) {
            m_id_to_pos[new_lrnid] = uint8_t(new_lrpos);
            for (size_type i = new_lrpos + 1; i < get_num_runs(); i++) {
                DEBUG_ABORT_IF_NE(m_id_to_pos[get_lrnid_from_position(i)], i - 1);
                m_id_to_pos[get_lrnid_from_position(i)] = uint8_t(i);
            }
        }

        return {new_lrnid, new_lrpos};
    }
    // Insert a new L-run with (new_chr, new_exp) after the L-run with lrnid.
    // Return the run-ID
    inline lcursor_type insert_after(lcursor_type& lcrsr, const uchar_type new_chr, const size_type new_exp) {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));

        const size_type new_lrnid = get_free_lrnid();
        const size_type new_lrpos = lcrsr.get_lrpos() + 1;

        unit_type new_unit;
        if constexpr (WITH_SAE) {
            new_unit = {new_chr, new_exp, nullptr, 0, new_lrnid, size_type(0)};
        } else {
            new_unit = {new_chr, new_exp, nullptr, 0, new_lrnid};
        }

        try_to_realloc(new_exp);
        make_free_unit(new_lrpos);
        set_unit_in_position(new_lrpos, new_unit);

        if constexpr (SUMEXP_HINT) {
            m_sum_exps += new_exp;
        }
        if constexpr (LOOKUP_HINT) {
            m_id_to_pos[new_lrnid] = uint8_t(new_lrpos);
            for (size_type i = new_lrpos + 1; i < get_num_runs(); i++) {
                DEBUG_ABORT_IF_NE(m_id_to_pos[get_lrnid_from_position(i)], i - 1);
                m_id_to_pos[get_lrnid_from_position(i)] = uint8_t(i);
            }
        }

        return {new_lrnid, new_lrpos};
    }
    // Split the L-run with lrnid into two runs, whose the first L-run is of size new_exp1.
    // Return the ID of the second (i.e., new created) L-run.
    inline lcursor_type split_after(lcursor_type& lcrsr, const size_type new_exp_a) {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        DEBUG_ABORT_IF_OUT(new_exp_a, 1, get_exp(lcrsr) - 1);

        const size_type new_lrnid = get_free_lrnid();
        const size_type new_lrpos = lcrsr.get_lrpos() + 1;

        unit_type new_unit_a = get_unit_from_pointer(get_unit_pointer(lcrsr.get_lrpos()), m_exp_bytes);
        unit_type new_unit_b = {new_unit_a.chr, new_unit_a.exp - new_exp_a, nullptr, 0, new_lrnid};

        set_exp_in_position(lcrsr.get_lrpos(), new_exp_a);

        make_free_unit(new_lrpos);
        set_unit_in_position(new_lrpos, new_unit_b);

        if constexpr (LOOKUP_HINT) {
            m_id_to_pos[new_lrnid] = uint8_t(new_lrpos);
            for (size_type i = new_lrpos + 1; i < get_num_runs(); i++) {
                DEBUG_ABORT_IF_NE(m_id_to_pos[get_lrnid_from_position(i)], i - 1);
                m_id_to_pos[get_lrnid_from_position(i)] = uint8_t(i);
            }
        }

        return {new_lrnid, new_lrpos};
    }

    /**
     *  Offset handlers
     */

    // Return the offset of run lrnid.
    inline offset_type get_offset(size_type lrnid) const {
        return get_offset(get_lcursor(lrnid));
    }
    inline offset_type get_offset(const lcursor_type lcrsr) const {
        DEBUG_ABORT_IF(!is_valid_lcursor(lcrsr));
        offset_type lofst = 0;
        for (size_type lrpos = 0; lrpos < lcrsr.get_lrpos(); lrpos++) {
            lofst += get_exp_from_position(lrpos);
        }
        return lofst;
    }
    inline std::pair<lcursor_type, offset_type> access_with_offset(offset_type lofst) {
        const size_type num_runs = get_num_runs();
        for (size_type lrpos = 0; lrpos < num_runs; lrpos++) {
            const offset_type exp = get_exp_from_position(lrpos);
            if (lofst < exp) {
                return {get_lcursor_from_position(lrpos), lofst};
            }
            lofst = lofst - exp;
        }
        ABORT_IF(true);
        return {};
    }

    /**
     *  Splitter
     */

    //! Return the data of after units and removes them from this group.
    //! Note that L-run's IDs will be changed, and the order and links of the returned group will not set.
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
            set_lrnid_in_position(i, i);
        }
        for (size_type i = 0; i < num_runs_aftr; i++) {
            data_aftr.set_lrnid_in_position(i, i);
        }

        /**
         *  Update hints
         */
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
                m_id_to_pos[get_lrnid_from_position(i)] = uint8_t(i);
            }
            std::fill_n(data_aftr.m_id_to_pos.data(), MAX_NUM_RUNS, UINT8_MAX);
            for (size_type i = 0; i < num_runs_aftr; i++) {
                data_aftr.m_id_to_pos[data_aftr.get_lrnid_from_position(i)] = uint8_t(i);
            }
        }

        return data_aftr;
    }

  private:
    inline size_type get_free_lrnid() const {
        return m_size_runs;
    }

    //
    inline void realloc_units(const uint8_t new_exp_bytes, const uint8_t new_size_runs, const uint8_t new_capa_runs) {
        auto resize_runs = std::min(m_size_runs, new_size_runs);
        auto new_serialized = std::make_unique<uint8_t[]>(get_unit_bytes(new_exp_bytes) * new_capa_runs);

        for (size_type lrpos = 0; lrpos < resize_runs; lrpos++) {
            uint8_t* new_ptr = new_serialized.get() + get_unit_begin(lrpos, new_exp_bytes);
            set_unit_in_pointer(new_ptr, get_unit_from_position(lrpos), new_exp_bytes);
        }

        m_units = std::move(new_serialized);
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

    //
    inline void make_free_unit(size_type lrpos) {
        if (m_size_runs == m_capa_runs) {
            realloc_units(m_exp_bytes, m_size_runs, std::min(m_capa_runs + GROW_SIZE, MAX_NUM_RUNS));
        }

        uint8_t* dst = m_units.get() + get_unit_begin(lrpos + 1, m_exp_bytes);
        uint8_t* src = m_units.get() + get_unit_begin(lrpos, m_exp_bytes);
        size_type bytes = get_unit_bytes(m_exp_bytes) * (m_size_runs - lrpos);

        std::memmove(dst, src, bytes);

        m_size_runs += 1;
    }

    inline uint8_t* get_unit_pointer(size_type lrpos) {
        return m_units.get() + get_unit_begin(lrpos, m_exp_bytes);
    }
    inline const uint8_t* get_unit_pointer(size_type lrpos) const {
        return m_units.get() + get_unit_begin(lrpos, m_exp_bytes);
    }

    inline size_type get_lrpos_from_identifier(size_type lrnid) const {
        if constexpr (LOOKUP_HINT) {
            DEBUG_ABORT_IF_EQ(m_id_to_pos[lrnid], UINT8_MAX);
            return m_id_to_pos[lrnid];
        } else {
            for (size_type lrpos = 0; lrpos < get_num_runs(); lrpos++) {
                if (get_lrnid_from_position(lrpos) == lrnid) {
                    return lrpos;
                }
            }
            return MAX_SIZE_INT;
        }
    }

    /**
     *  Accessors for units with position
     */
    inline uchar_type get_chr_from_position(size_type lrpos) const {
        return get_chr_from_pointer(get_unit_pointer(lrpos), m_exp_bytes);
    }
    inline size_type get_exp_from_position(size_type lrpos) const {
        return get_exp_from_pointer(get_unit_pointer(lrpos), m_exp_bytes);
    }
    inline fnode_type* get_fnode_from_position(size_type lrpos) const {
        return get_fnode_from_pointer(get_unit_pointer(lrpos), m_exp_bytes);
    }
    inline size_type get_frnid_from_position(size_type lrpos) const {
        return get_frnid_from_pointer(get_unit_pointer(lrpos), m_exp_bytes);
    }
    inline size_type get_lrnid_from_position(size_type lrpos) const {
        return get_lrnid_from_pointer(get_unit_pointer(lrpos), m_exp_bytes);
    }
    inline size_type get_sae_from_position(size_type lrpos) const {
        if constexpr (WITH_SAE) {
            return get_sae_from_pointer(get_unit_pointer(lrpos), m_exp_bytes);
        } else {
            ABORT_IF(true);
        }
    }

    inline void set_chr_in_position(size_type lrpos, uchar_type chr) {
        set_chr_in_pointer(get_unit_pointer(lrpos), chr, m_exp_bytes);
    }
    inline void set_exp_in_position(size_type lrpos, size_type exp) {
        set_exp_in_pointer(get_unit_pointer(lrpos), exp, m_exp_bytes);
    }
    inline void set_fnode_in_position(size_type lrpos, fnode_type* fnode) {
        set_fnode_in_pointer(get_unit_pointer(lrpos), fnode, m_exp_bytes);
    }
    inline void set_frnid_in_position(size_type lrpos, size_type frnid) {
        set_frnid_in_pointer(get_unit_pointer(lrpos), frnid, m_exp_bytes);
    }
    inline void set_lrnid_in_position(size_type lrpos, size_type lrnid) {
        set_lrnid_in_pointer(get_unit_pointer(lrpos), lrnid, m_exp_bytes);
    }
    inline void set_sae_in_position(size_type lrpos, size_type sae) {
        if constexpr (WITH_SAE) {
            set_sae_in_pointer(get_unit_pointer(lrpos), sae, m_exp_bytes);
        } else {
            ABORT_IF(true);
        }
    }

    inline unit_type get_unit_from_position(size_type lrpos) {
        return get_unit_from_pointer(get_unit_pointer(lrpos), m_exp_bytes);
    }
    inline void set_unit_in_position(size_type lrpos, const unit_type& unit) {
        set_unit_in_pointer(get_unit_pointer(lrpos), unit, m_exp_bytes);
    }

    /**
     *   Static members for handling serialized units
     */

    static inline size_type get_unit_begin(size_type lrpos, uint8_t exp_bytes) {
        return get_unit_bytes(exp_bytes) * lrpos;
    }
    static inline size_type get_unit_bytes(uint8_t exp_bytes) {
        // chr + exp + fnode + frnid + lrnid + sae
        if constexpr (WITH_SAE) {
            return sizeof(uchar_type) + exp_bytes + sizeof(uintptr_t) + sizeof(uint8_t) + sizeof(uint8_t) +
                   sizeof(size_type);
        } else {
            return sizeof(uchar_type) + exp_bytes + sizeof(uintptr_t) + sizeof(uint8_t) + sizeof(uint8_t);
        }
    }
    static inline uchar_type get_chr_from_pointer(const uint8_t* ptr, uint8_t) {
        return ptr[0];
    }
    static inline size_type get_exp_from_pointer(const uint8_t* ptr, uint8_t exp_bytes) {
        const size_type beg = sizeof(uchar_type);
        return utils::get_memory_block<size_type>(ptr + beg, exp_bytes);
    }
    static inline fnode_type* get_fnode_from_pointer(const uint8_t* ptr, uint8_t exp_bytes) {
        const size_type beg = sizeof(uchar_type) + exp_bytes;
        return reinterpret_cast<fnode_type*>(utils::get_memory_block<uintptr_t>(ptr + beg, sizeof(uintptr_t)));
    }
    static inline size_type get_frnid_from_pointer(const uint8_t* ptr, uint8_t exp_bytes) {
        const size_type beg = sizeof(uchar_type) + exp_bytes + sizeof(uintptr_t);
        return static_cast<size_type>(utils::get_memory_block<uint8_t>(ptr + beg, sizeof(uint8_t)));
    }
    static inline size_type get_lrnid_from_pointer(const uint8_t* ptr, uint8_t exp_bytes) {
        const size_type beg = sizeof(uchar_type) + exp_bytes + sizeof(uintptr_t) + sizeof(uint8_t);
        return static_cast<size_type>(utils::get_memory_block<uint8_t>(ptr + beg, sizeof(uint8_t)));
    }
    static inline size_type get_sae_from_pointer(const uint8_t* ptr, uint8_t exp_bytes) {
        if constexpr (WITH_SAE) {
            const size_type beg = sizeof(uchar_type) + exp_bytes + sizeof(uintptr_t) +  //
                                  sizeof(uint8_t) + sizeof(uint8_t);
            return utils::get_memory_block<size_type>(ptr + beg);
        } else {
            ABORT_IF(true);
        }
    }

    static inline void set_chr_in_pointer(uint8_t* ptr, uchar_type chr, uint8_t) {
        ptr[0] = chr;
    }
    static inline void set_exp_in_pointer(uint8_t* ptr, size_type exp, uint8_t exp_bytes) {
        const size_type beg = sizeof(uchar_type);
        utils::set_memory_block(ptr + beg, exp, exp_bytes);
    }
    static inline void set_fnode_in_pointer(uint8_t* ptr, fnode_type* fnode, uint8_t exp_bytes) {
        const size_type beg = sizeof(uchar_type) + exp_bytes;
        utils::set_memory_block(ptr + beg, reinterpret_cast<uintptr_t>(fnode), sizeof(uintptr_t));
    }
    static inline void set_frnid_in_pointer(uint8_t* ptr, size_type frnid, uint8_t exp_bytes) {
        const size_type beg = sizeof(uchar_type) + exp_bytes + sizeof(uintptr_t);
        utils::set_memory_block(ptr + beg, static_cast<uint8_t>(frnid), sizeof(uint8_t));
    }
    static inline void set_lrnid_in_pointer(uint8_t* ptr, size_type lrnid, uint8_t exp_bytes) {
        const size_type beg = sizeof(uchar_type) + exp_bytes + sizeof(uintptr_t) + sizeof(uint8_t);
        utils::set_memory_block(ptr + beg, static_cast<uint8_t>(lrnid), sizeof(uint8_t));
    }
    static inline void set_sae_in_pointer(uint8_t* ptr, size_type sae, uint8_t exp_bytes) {
        if constexpr (WITH_SAE) {
            const size_type beg = sizeof(uchar_type) + exp_bytes + sizeof(uintptr_t) +  //
                                  sizeof(uint8_t) + sizeof(uint8_t);
            utils::set_memory_block(ptr + beg, sae);
        } else {
            ABORT_IF(true);
        }
    }

    // TODO: Optimization
    static inline unit_type get_unit_from_pointer(uint8_t* ptr, uint8_t exp_bytes) {
        if constexpr (WITH_SAE) {
            return {get_chr_from_pointer(ptr, exp_bytes),  //
                    get_exp_from_pointer(ptr, exp_bytes),  //
                    get_fnode_from_pointer(ptr, exp_bytes),  //
                    get_frnid_from_pointer(ptr, exp_bytes),  //
                    get_lrnid_from_pointer(ptr, exp_bytes),  //
                    get_sae_from_pointer(ptr, exp_bytes)};
        } else {
            return {get_chr_from_pointer(ptr, exp_bytes),  //
                    get_exp_from_pointer(ptr, exp_bytes),  //
                    get_fnode_from_pointer(ptr, exp_bytes),  //
                    get_frnid_from_pointer(ptr, exp_bytes),  //
                    get_lrnid_from_pointer(ptr, exp_bytes)};
        }
    }

    static inline void set_unit_in_pointer(uint8_t* ptr, const unit_type& unit, uint8_t exp_bytes) {
        set_chr_in_pointer(ptr, unit.chr, exp_bytes);
        set_exp_in_pointer(ptr, unit.exp, exp_bytes);
        set_fnode_in_pointer(ptr, unit.fnode, exp_bytes);
        set_frnid_in_pointer(ptr, unit.frnid, exp_bytes);
        set_lrnid_in_pointer(ptr, unit.lrnid, exp_bytes);
        if constexpr (WITH_SAE) {
            set_sae_in_pointer(ptr, unit.sae, exp_bytes);
        }
    }
};

}  // namespace rcomp
