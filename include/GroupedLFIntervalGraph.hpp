/**
 * @file GroupedLFIntervalGraph.hpp
 */
#pragma once

#include <map>

#include "basics.hpp"
#include "defs.hpp"
#include "utils.hpp"

#include "GroupedFIndex.hpp"
#include "GroupedLIndex.hpp"

#include "StaticQueue.hpp"
#include "StaticVector.hpp"

namespace rcomp {

/**
 * An implementation of the grouped LF-interval graph.
 *
 * @tparam t_LData The data type of grouped L-node.
 * @tparam t_FData The data type of grouped F-node.
 * @tparam t_DivBound The threshould for dividing fat nodes.
 *
 * @note The technical terms used are
 *   - L/F-node means the L/F-run group consisting of some runs.
 *   - L/F-runId means the identifier of the run in the group.
 *   - L/F-cursor means the cursor of the run in the group.
 *   - L/F-offset means the offset of the character in the run.
 */
template <class t_LData, class t_FData, size_type t_DivBound = 7>
class GroupedLFIntervalGraph {
  public:
    using this_type = GroupedLFIntervalGraph<t_LData, t_FData, t_DivBound>;
    using lindex_type = GroupedLIndex<t_LData>;
    using findex_type = GroupedFIndex<t_FData>;
    using ldata_type = t_LData;
    using fdata_type = t_FData;
    using lnode_type = typename lindex_type::lnode_type;
    using fnode_type = typename findex_type::fnode_type;
    using lcursor_type = typename lnode_type::data_type::lcursor_type;
    using fcursor_type = typename fnode_type::data_type::fcursor_type;
    using lqueue_type = StaticQueue<lnode_type*, 4>;  // but enough with 3
    using fqueue_type = StaticQueue<fnode_type*, 4>;  // but enough with 3
    using lbuffer_type = StaticVector<lnode_type*, t_DivBound + 2>;
    using fbuffer_type = StaticVector<fnode_type*, t_DivBound + 2>;

    static constexpr size_type GROUP_BOUND = t_LData::GROUP_BOUND;
    static constexpr size_type DIV_BOUND = t_DivBound;

    static_assert(1 < GROUP_BOUND);
    static_assert(7 <= DIV_BOUND, "The division bound must be no less than 7");

    enum class insert_etype : uint8_t {
        // Case C
        MERGE_CURR,
        MERGE_PREV_IN,
        MERGE_PREV_OUT,
        // Case B
        SPLIT_FRONT,
        SPLIT_MIDDLE,
        SPLIT_BACK,
    };

    enum class bwres_etype : uint8_t {
        MATCHED,
        NOT_MATCHED,
        NOT_MATCHED_BY_EM,
        FAILED,
    };

    struct lcover_type {
        lnode_type* lnode;
        fnode_type* fnode;
        offset_type lofst;
        size_type lwght;
    };

    struct fcover_type {
        fnode_type* fnode;
        lnode_type* lnode;
        offset_type fofst;
        size_type fwght;
    };

    struct literator_type {
        lnode_type* lnode = nullptr;
        lcursor_type lcrsr = {};
        offset_type lofst = -1;

        literator_type() = default;

        literator_type(lnode_type* _lnode, size_type _lrnid, offset_type _lofst = -1)
            : lnode(_lnode), lcrsr(_lnode->get_data().get_fcursor(_lrnid)), lofst(_lofst) {}

        literator_type(lnode_type* _lnode, const lcursor_type& _lcrsr, offset_type _lofst = -1)
            : lnode(_lnode), lcrsr(_lcrsr), lofst(_lofst) {}

        literator_type(const std::tuple<lnode_type*, lcursor_type>& _lpair, offset_type _lofst = -1)
            : lnode(std::get<0>(_lpair)), lcrsr(std::get<1>(_lpair)), lofst(_lofst) {}

        inline bool is_valid() const {
            return lnode != nullptr;
        }
        inline bool is_valid_offset() const {
            return lofst >= 0;
        }
    };

    struct fiterator_type {
        fnode_type* fnode = nullptr;
        fcursor_type fcrsr = {};
        offset_type fofst = -1;

        fiterator_type() = default;

        fiterator_type(fnode_type* _fnode, size_type _frnid, offset_type _fofst = -1)
            : fnode(_fnode), fcrsr(_fnode->get_data().get_fcursor(_frnid)), fofst(_fofst) {}

        fiterator_type(fnode_type* _fnode, const fcursor_type& _fcrsr, offset_type _fofst = -1)
            : fnode(_fnode), fcrsr(_fcrsr), fofst(_fofst) {}

        fiterator_type(const std::tuple<fnode_type*, fcursor_type>& _fpair, offset_type _fofst = -1)
            : fnode(std::get<0>(_fpair)), fcrsr(std::get<1>(_fpair)), fofst(_fofst) {}

        inline bool is_valid() const {
            return fnode != nullptr;
        }
        inline bool is_valid_offset() const {
            return fofst >= 0;
        }
    };

  private:
    // The size of the input original text (with $).
    size_type m_num_chars = 0;

    // L/F data structures.
    lindex_type m_lindex;
    findex_type m_findex;

    // $'s iterator of L.
    lnode_type* m_em_lnode = nullptr;
    size_type m_em_lrnid = MAX_SIZE_INT;
    offset_type m_em_lofst = -1;

    // $'s iterator of F.
    fnode_type* m_em_fnode = nullptr;
    size_type m_em_frnid = MAX_SIZE_INT;
    offset_type m_em_fofst = -1;

    // Buffers for fat L/F-nodes.
    lqueue_type m_fat_lnodes;
    fqueue_type m_fat_fnodes;

    // Buffers for updated L/F-nodes.
    lbuffer_type m_lnodes_buffer;
    fbuffer_type m_fnodes_buffer;

#ifdef ENABLE_DEBUG_PRINT
    std::map<intptr_t, char> m_pcmap = {{intptr_t(nullptr), '?'}};
    char m_pcmax = 'A';
#endif

  public:
    //! Default constructor
    GroupedLFIntervalGraph() = default;

    //! Default destructor
    virtual ~GroupedLFIntervalGraph() = default;

    //! Copy constructor (deleted)
    GroupedLFIntervalGraph(const GroupedLFIntervalGraph&) = delete;

    //! Copy constructor (deleted)
    GroupedLFIntervalGraph& operator=(const GroupedLFIntervalGraph&) = delete;

    //! Move constructor
    GroupedLFIntervalGraph(GroupedLFIntervalGraph&&) noexcept = default;

    //! Move constructor
    GroupedLFIntervalGraph& operator=(GroupedLFIntervalGraph&&) noexcept = default;

    /* * * * * * * * * * * * * * * * * * *
     *
     *            Statistics
     *
     * * * * * * * * * * * * * * * * * * */

    //! Check if the data structure is empty.
    inline bool is_empty() const {
        return m_em_lnode == nullptr;
    }

    //! Get the number of stored characters including $.
    inline size_type get_num_chars() const {
        return m_num_chars;
    }

    //! Get the allocated memory in bytes.
    size_type get_memory_in_bytes(bool include_this = true) const {
        const size_type this_bytes = sizeof(*this) * include_this;
        return this_bytes + m_lindex.get_memory_in_bytes(false) + m_findex.get_memory_in_bytes(false);
    }

    //! Print the statistics related to memory.
    void show_memory_statistics() const {
        tfm::printfln("[Memory_Compressor]");
        tfm::reportfln("m_lindex:\t%d", m_lindex.get_memory_in_bytes());
        tfm::reportfln("m_findex:\t%d", m_findex.get_memory_in_bytes());
        m_lindex.show_memory_statistics();
        m_findex.show_memory_statistics();
    }

    //! Print the detailed statistics.
    void show_detailed_statistics() const {
        m_lindex.show_detailed_statistics();
        m_findex.show_detailed_statistics();
    }

    /* * * * * * * * * * * * * * * * * * *
     *
     *        Testers
     *
     * * * * * * * * * * * * * * * * * * */

    //! Test the alphabet order of F-index.
    void test_alphabet_order() const {
        m_findex.test_alphabet_order();
    }

    //! Test the list-labeling order of L-index.
    void test_list_order() const {
        m_lindex.test_order();
    }

    //! Test the LF-mapping of the LF-interval graph.
    void test_lf_mapping() const {
        auto curr_fitr = get_head_fiterator();
        ABORT_IF_NE(END_MARKER, get_chr(curr_fitr));

        set_next_fcursor(curr_fitr);

        while (true) {
            // Check corresponding nodes
            {
                const auto corr_litr = get_corr_lcursor(curr_fitr);
                const auto corr_fitr = get_corr_fcursor(corr_litr);
                ABORT_IF(curr_fitr.fnode != corr_fitr.fnode);
                ABORT_IF_NE(curr_fitr.fcrsr.get_frnid(), corr_fitr.fcrsr.get_frnid())
                ABORT_IF_NE(curr_fitr.fcrsr.get_frpos(), corr_fitr.fcrsr.get_frpos())
            }

            const auto next_fitr = get_next_fcursor(curr_fitr);
            if (!next_fitr.is_valid()) {
                break;
            }

            const uchar_type curr_chr = get_chr(curr_fitr);
            const uchar_type next_chr = get_chr(next_fitr);

            // F-node characters are not sorted in lex.
            ABORT_IF_LT(next_chr, curr_chr);

            if (curr_chr < next_chr) {
                curr_fitr = next_fitr;
                continue;
            }

            const auto curr_litr = get_corr_lcursor(curr_fitr);
            const auto next_litr = get_corr_lcursor(next_fitr);

            const loint_type curr_order = lo_common::get_basic_order(curr_litr.lnode);
            const loint_type next_order = lo_common::get_basic_order(next_litr.lnode);

            if (curr_order != next_order) {
                ABORT_IF_LE(next_order, curr_order);
            } else {
                ABORT_IF_LE(next_litr.lcrsr.get_lrpos(), curr_litr.lcrsr.get_lrpos());
            }

            curr_fitr = next_fitr;
        }
    }

    //! Test the Head/Tail links of L/F-nodes.
    void test_ht_links() const {
        const fnode_type* fnode = get_head_fnode();
        const lnode_type* lnode = get_head_lnode();

        offset_type fpos_b = 0;  // indicate the front of fnode
        offset_type lpos_b = 0;  // indicate the front of lnode

        while (!is_dmmy_lnode(lnode)) {
            const offset_type fpos_e = fpos_b + get_sum_fexps(fnode) - 1;  // indicate the back of fnode
            const offset_type lpos_e = lpos_b + get_sum_lexps(lnode) - 1;  // indicate the back of lnode
            if ((fpos_b <= lpos_b) and (lpos_b <= fpos_e) and (fpos_e <= lpos_e)) {
                ABORT_IF(fnode->get_data().get_tnode() != lnode);
                ABORT_IF(lnode->get_data().get_hnode() != fnode);
                ABORT_IF_NE(fpos_e - lpos_b, lnode->get_data().get_hofst());
                fpos_b = fpos_e + 1;
                lpos_b = lpos_e + 1;
                fnode = fnode->get_next();
                lnode = lnode->get_next();
            } else if (fpos_e < lpos_e) {
                ABORT_IF(fnode->get_data().get_tnode());
                fpos_b = fpos_e + 1;
                fnode = fnode->get_next();
            } else if (fpos_e > lpos_e) {
                ABORT_IF(lnode->get_data().get_hnode());
                lpos_b = lpos_e + 1;
                lnode = lnode->get_next();
            } else {
                ABORT_IF(true);
            }
        }

        while (!is_head_fnode(fnode)) {
            ABORT_IF(fnode->get_data().get_tnode());
            fpos_b = fpos_b + get_sum_fexps(fnode);
            fnode = fnode->get_next();
        }

        ABORT_IF_NE(fpos_b, lpos_b);
        ABORT_IF_NE(fpos_b, offset_type(m_num_chars));
    }

    //! Test the weights of L/F-nodes.
    void test_weights() const {
        for (fnode_type* fnode = get_head_fnode();; fnode = fnode->get_next()) {
            const fcover_type fcover = get_fcover(fnode);
            ABORT_IF_NE(fcover.fwght, fnode->get_data().get_weight());
            if (fnode->get_data().get_weight() >= DIV_BOUND) {
                ABORT_IF(!m_fat_fnodes.is_member(fnode));
            }
            if (fnode == get_tail_fnode()) {
                break;
            }
        }

        for (lnode_type* lnode = get_head_lnode();; lnode = lnode->get_next()) {
            const lcover_type lcover = get_lcover(lnode);
            ABORT_IF_NE(lcover.lwght, lnode->get_data().get_weight());
            if (lnode->get_data().get_weight() >= DIV_BOUND) {
                ABORT_IF(!m_fat_lnodes.is_member(lnode));
            }
            if (lnode == get_tail_lnode()) {
                break;
            }
        }
    }

    /* * * * * * * * * * * * * * * * * * *
     *
     *            Outputters
     *
     * * * * * * * * * * * * * * * * * * */

    /**
     * @brief Decode the original text in the insertion order.
     * @param[in] fn The callback function to get the original text charcater by charcater.
     */
    void decode_text(const std::function<void(uchar_type)>& fn) const {
        if (is_empty()) {
            return;
        }
        for (auto liter = get_head_literator(); !is_em_literator(liter);) {
            fn(get_chr(liter));
            liter = get_overlapped_literator(get_corr_fiterator(liter));
        }
    }

    /**
     * @brief Output the RLBWT text.
     * @param[in] fn The callback function to get the RLBWT text run by run.
     */
    void output_runs(const std::function<void(const run_type&)>& fn) const {
        // Setting END_MARKER is not problematic since the first run is never $.
        run_type prev_rn = {END_MARKER, 0};

        auto update_run = [&](const run_type& curr_rn) {
            if (prev_rn.chr == curr_rn.chr) {
                prev_rn.exp += curr_rn.exp;
                return;
            }
            if (prev_rn.exp != 0) {
                fn(prev_rn);
            }
            prev_rn = curr_rn;
        };

        for (auto liter = get_head_literator(); liter.is_valid();) {
            if (!is_em_lcursor(liter)) {
                update_run(make_run(get_chr(liter), get_lexp(liter)));
            } else {  // $-node
                const offset_type lexp = get_lexp(liter) - 1;  // without $
                if (m_em_lofst == 0) {  // Push front
                    update_run(make_run(END_MARKER, 1));
                    update_run(make_run(get_chr(liter), lexp));
                } else if (m_em_lofst == lexp) {  // Push back
                    update_run(make_run(get_chr(liter), lexp));
                    update_run(make_run(END_MARKER, 1));
                } else {
                    update_run(make_run(get_chr(liter), m_em_lofst));
                    update_run(make_run(END_MARKER, 1));
                    update_run(make_run(get_chr(liter), lexp - m_em_lofst));
                }
            }
            set_next_lcursor(liter);
        }
        fn(prev_rn);
    }

    /* * * * * * * * * * * * * * * * * * *
     *
     *        Initilization
     *
     * * * * * * * * * * * * * * * * * * */

    //! Make the data structure for the RLBWT text of T = (new_chr,$).
    void extend_init(const uchar_type new_chr) {
        // Dummy node for the partner of $'s F-node (for convenience)
        auto [dlr_lnode, dlr_lcrsr] = m_lindex.clear();
        // $'s node at the top of F
        auto dlr_fnode = std::get<0>(m_findex.clear(dlr_lnode, dlr_lcrsr.get_lrnid()));

        // Insert new nodes with new character
        auto [new_lnode, new_lcrsr] = m_lindex.insert_after(dlr_lnode, new_chr, 1);
        auto [new_fnode, new_fcrsr] = m_findex.insert_after(dlr_fnode, new_lnode, new_lcrsr.get_lrnid());

        // 'new_lnode' is always located at the head of L-list.
        // So, the T-link indicates 'dlr_fnode' ($-node) and will be never changed.
        dlr_fnode->get_data().set_tnode(new_lnode);
        new_lnode->get_data().set_hnode(dlr_fnode);
        new_lnode->get_data().set_hofst(0);

        new_fnode->get_data().set_tnode(nullptr);
        dlr_lnode->get_data().set_hnode(nullptr);
        dlr_lnode->get_data().set_hofst(0);

        m_em_lnode = new_lnode;
        m_em_lrnid = new_lcrsr.get_lrnid();
        m_em_lofst = 1;

        m_em_fnode = new_fnode;
        m_em_frnid = new_fcrsr.get_frnid();
        m_em_fofst = 0;

        dlr_fnode->get_data().set_weight(1);  // for new_lnode
        dlr_lnode->get_data().set_weight(0);
        new_fnode->get_data().set_weight(0);
        new_lnode->get_data().set_weight(2);  // for dlr_fnode and new_fnode

        m_num_chars = 2;
    }

    /* * * * * * * * * * * * * * * * * * *
     *
     *        Phase A in Case C/B
     *
     * * * * * * * * * * * * * * * * * * */

    //! Check if the insertion mode is in merge.
    inline bool is_mergable_mode(insert_etype ins_mode) const {
        return insert_etype::MERGE_CURR <= ins_mode and ins_mode <= insert_etype::MERGE_PREV_OUT;
    }

    //! Check if $'s L-node can be merged or split.
    inline insert_etype check_mergable(const uchar_type new_chr) const {
        DEBUG_ABORT_IF(is_dmmy_lnode(m_em_lnode));
        DEBUG_ABORT_IF_OUT(m_em_lofst, 0, get_lexp(m_em_lnode, m_em_lrnid) - 1);  // -1 is for $

        const auto em_liter = get_em_literator();

        // Mergeable with the current node?
        if (new_chr == get_chr(em_liter)) {
            return insert_etype::MERGE_CURR;  // Case C
        }

        if (is_first_char_in_lrun(em_liter)) {
            // Try to check the previous run
            const auto prev_liter = get_prev_lcursor(em_liter);
            if (!is_first_run_in_lnode(em_liter)) {
                // Mergeable with the previous run in the same group?
                if (new_chr == get_chr(prev_liter)) {
                    return insert_etype::MERGE_PREV_IN;  // Case C
                }
            } else {
                // Mergeable with the previous run in the previous group?
                if (new_chr == get_chr(prev_liter)) {
                    return insert_etype::MERGE_PREV_OUT;  // Case C
                }
            }
        }

        // new_chr will not be merged with the next run because $-mark is never placed at the end of some L-run.
        // Only when $-mark is placed at the last of L, it will be placed at the end of the last L-run.
        // But, then new_chr is not be merged with the next run anyway.

        if (is_first_char_in_lrun(em_liter)) {
            return insert_etype::SPLIT_FRONT;  // Case B
        } else if (is_last_char_in_lrun(em_liter)) {  // -1 is for $
            return insert_etype::SPLIT_BACK;  // Case B
        }
        return insert_etype::SPLIT_MIDDLE;  // Case B
    }

    /* * * * * * * * * * * * * * * * * * *
     *
     *       Phase B -- F in Case C
     *
     * * * * * * * * * * * * * * * * * * */

    //! Extend the RLBWT text in Case C.
    void extend_with_merge(const uchar_type, const insert_etype ins_mode) {
        DEBUG_ABORT_IF(!is_mergable_mode(ins_mode));
        DEBUG_PRINT(tfm::printfln("\n==== Case C ====");)

        DEBUG_ABORT_IF(!m_fat_lnodes.is_empty());
        DEBUG_ABORT_IF(!m_fat_fnodes.is_empty());

        m_num_chars += 1;

        // Variables for updated L/F-nodes
        literator_type ins_liter;
        fiterator_type ins_fiter;

        // [NOTE]
        // Depending on 'ins_mode', the L-node with $ can be updated.
        // But Do NOT change the value of 'm_em_lnode' until Phase E.
        // This is because this change will affect 'get_sum_lexps' (in Phase C).

        //
        // Phase B
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase B =="););

        const auto em_liter = get_em_literator();

        if (ins_mode == insert_etype::MERGE_CURR) {
            ins_liter = em_liter;
            ins_fiter = get_corr_fcursor(ins_liter);
            ins_fiter.fofst = m_em_lofst;
        } else {  // ins_mode == MERGE_PREV_IN or MERGE_PREV_OUT
            ins_liter = get_prev_lcursor(em_liter);
            ins_fiter = get_corr_fcursor(ins_liter);
            ins_fiter.fofst = get_fexp(ins_fiter);
        }

        //
        // Phase C (Remove H/T-links of overlapped nodes)
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase C ==");)

        // $-replacement
        if (ins_mode == insert_etype::MERGE_PREV_OUT) {
            const auto em_fiter = get_em_fiterator();
            const bool is_first_em_fofst = is_first_char_in_fnode(em_fiter);
            const bool is_last_em_fofst = is_last_char_in_fnode(em_fiter);

            // Update the weights
            if (is_first_em_fofst) {
                lnode_type* prev_lnode = m_em_lnode->get_prev();
                sub_weight(m_em_lnode, 1);
                add_weight(prev_lnode, 1);
                push_fat_lnode(prev_lnode);
            }

            if (is_last_em_fofst) {
                fnode_type* next_fnode = m_em_fnode->get_next();
                sub_weight(m_em_fnode, 1);
                add_weight(next_fnode, 1);
                push_fat_fnode(next_fnode);
            }

            // Update the linkage
            lnode_type* em_tnode = m_em_fnode->get_data().get_tnode();

            if (is_last_em_fofst) {
                DEBUG_ABORT_IF(em_tnode != m_em_lnode);
                DEBUG_ABORT_IF_NE(em_tnode->get_data().get_hofst(), 0);
                {
                    fnode_type* fnode_befo = m_em_fnode;
                    lnode_type* lnode_befo = m_em_lnode->get_prev();
                    const offset_type new_sum_fexps = fnode_befo->get_data().get_sum_exps();
                    const offset_type new_sum_lexps = lnode_befo->get_data().get_sum_exps() + 1;
                    if (new_sum_lexps <= new_sum_fexps) {
                        DEBUG_ABORT_IF(lnode_befo->get_data().get_hnode());
                        lnode_befo->get_data().set_hlink(fnode_befo, new_sum_lexps - 1);
                        fnode_befo->get_data().set_tnode(lnode_befo);
                    } else {
                        DEBUG_ABORT_IF(!lnode_befo->get_data().get_hnode());
                        fnode_befo->get_data().reset_tlink();
                    }
                }
                {
                    fnode_type* fnode_aftr = m_em_fnode->get_next();
                    lnode_type* lnode_aftr = m_em_lnode;
                    const offset_type new_sum_fexps = fnode_aftr->get_data().get_sum_exps();
                    const offset_type new_sum_lexps = lnode_aftr->get_data().get_sum_exps();
                    if (new_sum_lexps >= new_sum_fexps) {
                        DEBUG_ABORT_IF(fnode_aftr->get_data().get_tnode());
                        lnode_aftr->get_data().set_hlink(fnode_aftr, new_sum_fexps - 1);
                        fnode_aftr->get_data().set_tnode(lnode_aftr);
                    } else {
                        DEBUG_ABORT_IF(!fnode_aftr->get_data().get_tnode());
                        lnode_aftr->get_data().reset_hlink();
                    }
                }
            } else if (em_tnode == m_em_lnode) {
                // The length of L-node with $-mark should be no less than 2.
                // That is, the tofst should be no less than 1.
                DEBUG_ABORT_IF_LT(em_tnode->get_data().get_hofst(), 1);
                em_tnode->get_data().set_hofst(em_tnode->get_data().get_hofst() - 1);
            }
        }

        //
        // Phase D
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase D ==");)

        // ins_lnode->get_data().prefetch();
        // ins_fnode->get_data().prefetch();

        add_exp(ins_liter, 1);
        add_exp(ins_fiter, 1);

        //
        // Phase E
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase E ==");)

        reset_em_lvars();  // Clear because the L-run has already updated.
        set_em_fvars(ins_fiter);  // Cache to the next process
        set_em_lvars(get_overlapped_literator(ins_fiter));  // $-query

        //
        // Phase F
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase F ==");)

        // ?-replacement
        {
            const bool is_last_em_fofst = is_last_char_in_fnode(ins_fiter);  // em_fiter
            const bool is_first_em_lofst = is_first_char_in_lnode(get_em_literator());

            if (is_last_em_fofst and is_first_em_lofst) {
                // ?-2-B
                // Update weights
                add_weight(m_em_fnode, 1);
                sub_weight(m_em_fnode->get_next(), 1);
                push_fat_fnode(m_em_fnode);
                // Update Linkage
                remove_fnode_link(m_em_fnode);
                remove_lnode_link(m_em_lnode);
                m_em_lnode->get_data().set_hlink(m_em_fnode, 0);
                m_em_fnode->get_data().set_tnode(m_em_lnode);
            } else if (m_em_lnode->get_data().get_hnode() == m_em_fnode) {
                // ?-3-A
                const offset_type new_hofst = m_em_lnode->get_data().get_hofst() + 1;
                m_em_lnode->get_data().set_hofst(new_hofst);
            }
        }
    }

    /* * * * * * * * * * * * * * * * * * *
     *
     *       Phase B -- F in Case B
     *
     * * * * * * * * * * * * * * * * * * */

    // Extend the BWT text in Case B.
    void extend_with_split(const uchar_type new_chr, const insert_etype ins_mode) {
        DEBUG_ABORT_IF(is_mergable_mode(ins_mode));

        DEBUG_ABORT_IF(!m_fat_lnodes.is_empty());
        DEBUG_ABORT_IF(!m_fat_fnodes.is_empty());

        DEBUG_PRINT(tfm::printfln("\n==== Case B ====");)

        // Note:
        // In Case B, L-interval is never changed because the split will be done within the group.

        m_num_chars += 1;

        //
        // Phase B
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase B ==");)

        // Variables for the predecessor F-run
        fiterator_type pred_fiter;
        {
            auto pred_liter = get_em_literator();
            if (ins_mode == insert_etype::SPLIT_FRONT) {
                set_prev_lcursor(pred_liter);
            }
            pred_fiter = predecessor(new_chr, pred_liter);
        }

        // In Case B, the three patterns should be considered.
        //  - INSERT_BEFORE: A new F-run is inserted before F-run (ins_fnode, ins_fcrsr) within the group.
        //  - INSERT_AFTER : A new F-run is inserted after F-run (ins_fnode, ins_fcrsr) within the group.
        //  - NEW_CREATE   : A new F-run is inserted before F-group ins_fnode as a new group.
        enum class fnode_insert_etype : uint8_t { INSERT_BEFORE, INSERT_AFTER, NEW_CREATE };
        auto fins_mode = fnode_insert_etype::INSERT_BEFORE;

        // A new F-run will be inserted before/after F-run (ins_fnode, ins_fcrsr).
        // auto [ins_fnode, ins_fcrsr] = get_next_fcursor(pred_fnode, pred_fcrsr);
        auto ins_fiter = get_next_fcursor(pred_fiter);

        if (ins_fiter.is_valid() and new_chr == get_chr(ins_fiter)) {
            fins_mode = fnode_insert_etype::INSERT_BEFORE;
        } else if (new_chr == get_chr(pred_fiter)) {
            ins_fiter = pred_fiter;
            fins_mode = fnode_insert_etype::INSERT_AFTER;
        } else {
            fins_mode = fnode_insert_etype::NEW_CREATE;
        }

        // Variables for the F-intervals to be inserted.
        // i.e., the F-interval is of size 'ins_fsize' starting at the 'ins_fofst'-th position in group 'ins_fnode'.

        //
        // Phase C
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase C ==");)

        // Nothing to do

        //
        // Phase D
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase D ==");)

        // Variables for indicating a new L-run inserted
        auto ins_liter = get_em_literator();
        auto div_fiter = fiterator_type{};

        // Insert a new L-run
        if (ins_mode == insert_etype::SPLIT_FRONT) {
            insert_lrun_before(ins_liter, new_chr, 1);
        } else if (ins_mode == insert_etype::SPLIT_BACK) {
            insert_lrun_after(ins_liter, new_chr, 1);
        } else if (ins_mode == insert_etype::SPLIT_MIDDLE) {
            div_fiter = get_corr_fcursor(ins_liter);
            split_run_after(ins_liter, div_fiter, m_em_lofst);
            insert_lrun_before(ins_liter, new_chr, 1);
        } else {
            ABORT_IF(true);
        }

        // Insert a new F-run
        if (fins_mode == fnode_insert_etype::INSERT_BEFORE) {
            insert_frun_before(ins_fiter, ins_liter);
        } else if (fins_mode == fnode_insert_etype::INSERT_AFTER) {
            insert_frun_after(ins_fiter, ins_liter);
        } else if (fins_mode == fnode_insert_etype::NEW_CREATE) {
            ins_fiter = {m_findex.insert_after(pred_fiter.fnode, ins_liter.lnode, ins_liter.lcrsr.get_lrnid())};
        } else {
            ABORT_IF(true);
        }
        ins_fiter.fofst = 0;

        //
        // Phase E
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase E ==");)

        reset_em_lvars();  // Clear because the L-run has already updated.
        set_em_fvars(ins_fiter);  // Cache to the next process
        set_em_lvars(get_overlapped_literator(ins_fiter));  // $-query

        // Phase F
        DEBUG_PRINT(tfm::printfln("\n== Phase F ==");)

        // ?-replacement
        {
            const bool is_last_em_fofst = is_last_char_in_fnode(get_em_fiterator());
            const bool is_first_em_lofst = is_first_char_in_lnode(get_em_literator());

            if (fins_mode != fnode_insert_etype::NEW_CREATE) {
                // Then, the processing is the same as that of Case C.
                if (is_last_em_fofst and is_first_em_lofst) {
                    // ?-2-B
                    // Update OV
                    m_em_fnode->get_data().add_weight(1);
                    m_em_fnode->get_next()->get_data().sub_weight(1);
                    push_fat_fnode(m_em_fnode);
                    // Update Linkage
                    remove_fnode_link(m_em_fnode);
                    remove_lnode_link(m_em_lnode);
                    m_em_lnode->get_data().set_hlink(m_em_fnode, 0);
                    m_em_fnode->get_data().set_tnode(m_em_lnode);
                } else if (m_em_lnode->get_data().get_hnode() == m_em_fnode) {
                    // ?-3-A
                    const offset_type new_hofst = m_em_lnode->get_data().get_hofst() + 1;
                    m_em_lnode->get_data().set_hofst(new_hofst);
                }
            } else {  // fins_mode == fnode_insert_etype::NEW_CREATE
                if (is_first_em_lofst) {
                    // Update OV
                    m_em_fnode->get_data().set_weight(1);
                    m_em_fnode->get_next()->get_data().sub_weight(1);
                    // Update Linkage
                    remove_lnode_link(m_em_lnode);
                    m_em_lnode->get_data().set_hlink(m_em_fnode, 0);
                    m_em_fnode->get_data().set_tnode(m_em_lnode);
                }
                m_em_lnode->get_data().add_weight(1);  // for new m_em_fnode
                push_fat_lnode(m_em_lnode);
            }
        }

        //
        // Extra Phase: Grouped split
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase X (Extra Phase) ==");)

        // DEBUG_ABORT_IF(!ins_lnode);
        if (ins_liter.lnode->get_data().get_num_runs() >= GROUP_BOUND) {
            fnode_type* updated_fnode = divide_lnode(ins_liter.lnode);
            push_fat_fnode(updated_fnode);
        }
        // DEBUG_ABORT_IF(!ins_fnode);
        if (ins_fiter.fnode->get_data().get_num_runs() >= GROUP_BOUND) {
            lnode_type* updated_lnode = divide_fnode(ins_fiter.fnode);
            push_fat_lnode(updated_lnode);
        }
        if (div_fiter.is_valid() and div_fiter.fnode->get_data().get_num_runs() >= GROUP_BOUND) {
            lnode_type* updated_lnode = divide_fnode(div_fiter.fnode);
            push_fat_lnode(updated_lnode);
        }
    }

    /* * * * * * * * * * * * * * * * * * *
     *
     *      Phase B -- F in Case L/F
     *
     * * * * * * * * * * * * * * * * * * */

    //! Divide the fat nodes in Case L/F.
    bool divide_fat_node() {
        // Case L
        while (!m_fat_lnodes.is_empty()) {
            lnode_type* lnode = m_fat_lnodes.pop();
            DEBUG_ABORT_IF(!lnode);
            if (lnode->get_data().get_weight() >= DIV_BOUND) {
                const lcover_type lcover = get_lcover(lnode);
                DEBUG_ABORT_IF_NE(lcover.lwght, lnode->get_data().get_weight());
                fnode_type* div_fnode = divide_lcover(lcover);
                if (div_fnode and div_fnode->get_data().get_num_runs() >= GROUP_BOUND) {
                    push_fat_lnode(divide_fnode(div_fnode));
                }
                return true;
            }
        }

        // Case F
        while (!m_fat_fnodes.is_empty()) {
            fnode_type* fnode = m_fat_fnodes.pop();
            DEBUG_ABORT_IF(!fnode);
            if (fnode->get_data().get_weight() >= DIV_BOUND) {
                const fcover_type fcover = get_fcover(fnode);
                DEBUG_ABORT_IF_NE(fcover.fwght, fnode->get_data().get_weight());
                lnode_type* div_lnode = divide_fcover(fcover);
                if (div_lnode and div_lnode->get_data().get_num_runs() >= GROUP_BOUND) {
                    push_fat_fnode(divide_lnode(div_lnode));
                }
                return true;
            }
        }

        return false;
    }

    /**
     * @brief Divide a fat L-node (Case L).
     * @param[out] divided The before part of the divided L-nodes if a L-run is split in the fat L-node.
     * @return true if a fat L-node is divided, or false otherwise.
     * @note The division can produce a new fat L-node due to the L-run split.
     */
    inline bool divide_fat_lnode(lnode_type*& divided) {
        while (!m_fat_lnodes.is_empty()) {
            lnode_type* lnode = m_fat_lnodes.pop();
            DEBUG_ABORT_IF(!lnode);
            if (lnode->get_data().get_weight() >= DIV_BOUND) {
                const lcover_type lcover = get_lcover(lnode);
                DEBUG_ABORT_IF_NE(lcover.lwght, lnode->get_data().get_weight());
                fnode_type* div_fnode = divide_lcover(lcover);
                if (div_fnode and div_fnode->get_data().get_num_runs() >= GROUP_BOUND) {
                    push_fat_lnode(divide_fnode(div_fnode));
                }
                divided = div_fnode ? lcover.lnode : nullptr;
                return true;
            }
        }
        return false;
    }

    /**
     * @brief Divide a fat F-node (Case F).
     * @param[out] divided The before part of the divided F-nodes if a F-run is split in the fat F-node.
     * @return true if a fat F-node is divided, or false otherwise.
     * @note The division can produce a new fat F-node due to the F-run split.
     */
    inline bool divide_fat_fnode(fnode_type*& divided) {
        while (!m_fat_fnodes.is_empty()) {
            fnode_type* fnode = m_fat_fnodes.pop();
            DEBUG_ABORT_IF(!fnode);
            if (fnode->get_data().get_weight() >= DIV_BOUND) {
                const fcover_type fcover = get_fcover(fnode);
                DEBUG_ABORT_IF_NE(fcover.fwght, fnode->get_data().get_weight());
                lnode_type* div_lnode = divide_fcover(fcover);
                if (div_lnode and div_lnode->get_data().get_num_runs() >= GROUP_BOUND) {
                    push_fat_fnode(divide_lnode(div_lnode));
                }
                divided = div_lnode ? fcover.fnode : nullptr;
                return true;
            }
        }
        return false;
    }

    /**
     * @brief Divide the fat L-node whose weight is no less than DIV_BOUND.
     * @param[in] lcover The fat L-node and the overlapped information.
     * @return The F-node in which a F-run is split due to the L-run split (as noted).
     * @note The division can split a L-run in the fat L-node, resulting in splitting the corresponding F-run.
     */
    fnode_type* divide_lcover(const lcover_type& lcover) {
        DEBUG_ABORT_IF(!m_lnodes_buffer.is_empty());
        DEBUG_ABORT_IF(!m_fnodes_buffer.is_empty());
        DEBUG_ABORT_IF_LT(lcover.lwght, DIV_BOUND);

        if (!lcover.lnode) {
            return nullptr;
        }

        DEBUG_PRINT(tfm::printfln("\n==== Case L ====");)
        DEBUG_PRINT(tfm::printfln("lcover => %d", get_lcover_str(lcover));)

        //
        // Phase B: Search the L-run to be divided.
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase B ==");)

        // The target L-node to be divided
        lnode_type* div_lnode = lcover.lnode;

        // If div_lnode == m_em_lnode, may need to update the information of $'s L-node.
        // For this, in advance keep the offset of $'s position in the L-node.
        offset_type em_lofst_in_node = -1;
        if (div_lnode == m_em_lnode) {
            // The offset will be correct although get_offset() does not consider the imaginary $-mark.
            em_lofst_in_node = m_em_lnode->get_data().get_offset(m_em_lrnid) + m_em_lofst;
        }

        // Search the L-run to be divided.
        // div_lofst_in_node will indicate the boundary position to divide, corresponding to the last of before part.
        offset_type div_lofst_in_node = lcover.lofst;
        {
            const size_type fnum_befo = lcover.lwght / 2;
            const fnode_type* fnode = lcover.fnode;

            for (size_type i = 1; i < fnum_befo; i++) {
                fnode = fnode->get_next();
                div_lofst_in_node += get_sum_fexps(fnode);
            }

            // SPECIAL CASE: Need to shift div_lofst_in_node?
            if (em_lofst_in_node == div_lofst_in_node) {
                // Then, $-marker will be placed at the end of the divided L-node (before part).
                // To avoid this, shift the division point.
                fnode = fnode->get_next();
                div_lofst_in_node += get_sum_fexps(fnode);
                ABORT_IF_LE(get_sum_lexps(div_lnode), div_lofst_in_node);
            }
        }
        // Here, div_lnode[div_lofst_in_node] is the last character of the before part.

        DEBUG_PRINT(tfm::printfln("div_lofst_in_node=%d", div_lofst_in_node);)
        DEBUG_PRINT(tfm::printfln("em_lofst_in_node=%d", em_lofst_in_node);)

        //
        // Phase C: Remove head links of divided nodes.
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase C ==");)

        get_overlapped_fnodes(div_lnode, 0, get_sum_lexps(div_lnode), m_fnodes_buffer);
        remove_fnode_links(m_fnodes_buffer);

        //
        // Phase D: Divide the L-node.
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase D ==");)
        ABORT_IF_EQ(em_lofst_in_node, div_lofst_in_node);

        // 1) Split the target L-run in the node.
        lcursor_type div_lcrsr_befo = {};
        offset_type div_lofst_befo = -1;

        if (div_lnode == m_em_lnode and em_lofst_in_node < div_lofst_in_node) {
            // Then, $-mark will be placed at the before part, and consider the imaginary $-mark in access_with_offset.
            std::tie(div_lcrsr_befo, div_lofst_befo) = div_lnode->get_data().access_with_offset(div_lofst_in_node - 1);
        } else {
            std::tie(div_lcrsr_befo, div_lofst_befo) = div_lnode->get_data().access_with_offset(div_lofst_in_node);
        }
        const size_type div_lexp_befo = div_lofst_befo + 1;

        DEBUG_PRINT(tfm::printfln("div_lcrsr_befo=%s", get_lcrsr_str(div_lcrsr_befo));)
        DEBUG_PRINT(tfm::printfln("div_lofst_befo=%d, div_lexp_befo=%d", div_lofst_befo, div_lexp_befo);)

        // Need to consider the updated F-node due to the L-run division.
        fnode_type* div_fnode = nullptr;

        // The L-run needs to be split? (Do NOT use get_lexp not to consider $-marker).
        if (div_lnode->get_data().get_exp(div_lcrsr_befo) != div_lexp_befo) {
            std::tie(div_fnode, std::ignore) = split_lrun(div_lnode, div_lcrsr_befo, div_lexp_befo);
        }

        // 2) Divide the L-node. new_lnode is the after part.
        lnode_type* new_lnode = m_lindex.divide_after(div_lnode, div_lcrsr_befo);

        DEBUG_PRINT(tfm::printfln("div_lnode=%s, new_lnode=%s", get_pc(div_lnode), get_pc(new_lnode));)
        DEBUG_PRINT(tfm::printfln("div_fnode=%s", get_pc(div_fnode));)

        //
        // Phase E
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase E ==");)

        if (div_lnode == m_em_lnode) {
            DEBUG_ABORT_IF_LT(em_lofst_in_node, 0);

            if (em_lofst_in_node < div_lofst_in_node) {
                // Then, $-mark will be placed at div_lnode (before part)
                m_em_lnode = div_lnode;
            } else if (div_lofst_in_node < em_lofst_in_node) {
                // Then, $-mark will be placed at new_lnode (after part)
                DEBUG_ABORT_IF_NE(size_type(div_lofst_in_node + 1), div_lnode->get_data().get_sum_exps());
                m_em_lnode = new_lnode;
                em_lofst_in_node -= (div_lofst_in_node + 1);
            } else {
                // Should be avoided in Phase B
                ABORT_IF_EQ(em_lofst_in_node, div_lofst_in_node);
            }

            if (m_em_lnode == get_tail_lnode() and
                size_type(em_lofst_in_node) == m_em_lnode->get_data().get_sum_exps()) {
                // Only in this case, $-mark is placed at the end (i.e. exclusive position)
                // So, should not call access_with_offset
                const lcursor_type em_lcrsr = m_em_lnode->get_data().get_last_lcursor();

                m_em_lrnid = em_lcrsr.get_lrnid();
                m_em_lofst = m_em_lnode->get_data().get_exp(em_lcrsr);
            } else {
                const auto [em_lcrsr, em_lofst] = m_em_lnode->get_data().access_with_offset(em_lofst_in_node);

                m_em_lrnid = em_lcrsr.get_lrnid();
                m_em_lofst = em_lofst;
            }
        }

        //
        // Phase F
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase F ==");)

        // In this call, div_lnode and new_lnode will be linked appropriately.
        reset_fnode_links_in_order(m_fnodes_buffer);
        m_fnodes_buffer.clear();

        {
            lcover_type div_lcover = get_lcover(div_lnode);
            div_lnode->get_data().set_weight(div_lcover.lwght);
        }
        {
            lcover_type new_lcover = get_lcover(new_lnode);
            new_lnode->get_data().set_weight(new_lcover.lwght);
            new_lcover.fnode->get_data().add_weight(1);  // for new_lnode
        }

        return div_fnode;
    }

    /**
     * @brief Devide the fat F-node whose weight is no less than DIV_BOUND.
     * @param[in] fcover The fat F-node and the overlapped information.
     * @return The L-node in which a L-run is split due to the F-run split (as noted).
     * @note The division can split a F-run in the fat F-node, resulting in splitting the corresponding L-run.
     */
    lnode_type* divide_fcover(const fcover_type& fcover) {
        DEBUG_ABORT_IF(!m_lnodes_buffer.is_empty());
        DEBUG_ABORT_IF(!m_fnodes_buffer.is_empty());
        DEBUG_ABORT_IF_LT(fcover.fwght, DIV_BOUND);

        if (!fcover.fnode) {
            return nullptr;
        }

        DEBUG_PRINT(tfm::printfln("\n==== Case F ====");)
        DEBUG_PRINT(tfm::printfln("fcover => %d", get_fcover_str(fcover));)

        //
        // Phase B: Search the F-run to be split.
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase B ==");)

        fnode_type* div_fnode = fcover.fnode;
        offset_type div_fofst_in_node = fcover.fofst;
        {
            const size_type lnum_befo = fcover.fwght / 2;
            const lnode_type* lnode = fcover.lnode;
            for (size_type i = 1; i < lnum_befo; i++) {
                lnode = lnode->get_next();
                div_fofst_in_node += get_sum_lexps(lnode);
            }
        }
        // Here, div_fnode[div_fofst_in_node] is the last character of the front part.

        // If div_fnode == m_em_fnode, need to update m_em_f*.
        // So, in advance keep the offset of $ position in the F-node.
        offset_type em_fofst_in_node = -1;
        if (div_fnode == m_em_fnode) {
            em_fofst_in_node = m_em_fnode->get_data().get_offset(m_em_frnid) + m_em_fofst;
        }

        //
        // Phase C: Remove head links of divided nodes.
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase C ==");)
        get_overlapped_lnodes(div_fnode, 0, get_sum_fexps(div_fnode), m_lnodes_buffer);
        remove_lnode_links(m_lnodes_buffer);

        //
        // Phase D: Update nodes
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase D ==");)

        // 1) Split the target F-run in the node.
        auto [div_fcrsr_befo, div_fofst_befo] = div_fnode->get_data().access_with_offset(div_fofst_in_node);
        const size_type div_fexp_befo = div_fofst_befo + 1;

        // Need to consider OV change depending on F-run split
        lnode_type* div_lnode = nullptr;
        if (get_fexp(div_fnode, div_fcrsr_befo) != offset_type(div_fexp_befo)) {  // the F-run needs to be split?
            std::tie(div_lnode, std::ignore) = split_frun(div_fnode, div_fcrsr_befo, div_fexp_befo);
        }

        // 2) Divide the F-node
        fnode_type* new_fnode = m_findex.divide_after(div_fnode, div_fcrsr_befo);

        //
        // Phase E
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase E ==");)

        if (div_fnode == m_em_fnode) {
            DEBUG_ABORT_IF_LT(em_fofst_in_node, 0);

            const offset_type sum_fexps = div_fnode->get_data().get_sum_exps();

            if (em_fofst_in_node < sum_fexps) {
                m_em_fnode = div_fnode;
            } else {
                m_em_fnode = new_fnode;
                em_fofst_in_node -= sum_fexps;
            }

            const auto [em_fcrsr, em_fofst] = m_em_fnode->get_data().access_with_offset(em_fofst_in_node);

            m_em_frnid = em_fcrsr.get_frnid();
            m_em_fofst = em_fofst;
        }

        //
        // Phase F
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase F ==");)

        // In this call, div_fnode and new_fnode will be linked appropriately.
        reset_lnode_links_in_order(m_lnodes_buffer);
        // Release
        m_lnodes_buffer.clear();

        {
            fcover_type div_fcover = get_fcover(div_fnode);
            div_fnode->get_data().set_weight(div_fcover.fwght);
        }
        {
            fcover_type new_fcover = get_fcover(new_fnode);
            new_fnode->get_data().set_weight(new_fcover.fwght);
            new_fcover.lnode->get_data().add_weight(1);  // for new_fnode
        }

        return div_lnode;
    }

    /**
     * @brief Divide the given L-node into two L-nodes in the middle.
     * @param[in] div_lnode L-node to be divided
     * @return F-node whose weight is updated
     */
    fnode_type* divide_lnode(lnode_type* div_lnode) {
        DEBUG_ABORT_IF(is_dmmy_lnode(div_lnode));
        DEBUG_ABORT_IF_LT(div_lnode->get_data().get_num_runs(), GROUP_BOUND);

        DEBUG_ABORT_IF(!m_lnodes_buffer.is_empty());
        DEBUG_ABORT_IF(!m_fnodes_buffer.is_empty());

        DEBUG_PRINT(tfm::printfln("\n** divide_lnode **");)
        DEBUG_PRINT(tfm::printfln("div_lnode=%s", get_pc(div_lnode));)

        // Keep the original $-mark L-position for Phase E
        const lcursor_type em_lcrsr = m_em_lnode->get_data().get_lcursor(m_em_lrnid);

        // Phase C
        get_overlapped_fnodes(div_lnode, 0, get_sum_lexps(div_lnode), m_fnodes_buffer);
        remove_fnode_links(m_fnodes_buffer);

        DEBUG_PRINT(tfm::printfln("\n== Phase C ==");)

        // Phase D
        DEBUG_PRINT(tfm::printfln("\n== Phase D ==");)

        // div_lnode is divided into (div_lnode, new_lnode)
        lnode_type* new_lnode = m_lindex.divide_after(div_lnode);

        DEBUG_PRINT(tfm::printfln("div_lnode=%s", get_lnode_str(div_lnode));)
        DEBUG_PRINT(tfm::printfln("new_lnode=%s", get_lnode_str(new_lnode));)
        DEBUG_PRINT(tfm::printfln("m_em_lnode=%s, m_em_lrnid=%d", get_pc(m_em_lnode), m_em_lrnid);)

        // Phase E
        DEBUG_PRINT(tfm::printfln("\n== Phase E ==");)

        if (m_em_lnode == div_lnode) {
            const size_type em_lrpos = em_lcrsr.get_lrpos();
            const size_type num_runs = div_lnode->get_data().get_num_runs();
            if (em_lrpos < num_runs) {  // $-run is in the front part?
                m_em_lnode = div_lnode;
                m_em_lrnid = m_em_lnode->get_data().get_lcursor_from_position(em_lrpos).get_lrnid();
            } else {  // $-run is in the rear part?
                m_em_lnode = new_lnode;
                m_em_lrnid = m_em_lnode->get_data().get_lcursor_from_position(em_lrpos - num_runs).get_lrnid();
            }
        }

        //
        // Phase F
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase F ==");)

        // In this call, div_lnode and new_lnode will be linked appropriately.
        reset_fnode_links_in_order(m_fnodes_buffer);
        // Release
        m_fnodes_buffer.clear();

        {
            const lcover_type lcover = get_lcover(div_lnode);
            div_lnode->get_data().set_weight(lcover.lwght);
        }
        {
            const lcover_type lcover = get_lcover(new_lnode);
            new_lnode->get_data().set_weight(lcover.lwght);
        }

        fnode_type* updated_fnode = std::get<0>(get_first_overlapped_fnode(new_lnode));
        updated_fnode->get_data().add_weight(1);  // for new_lnode

        return updated_fnode;
    }

    /**
     * @brief Divide the given F-node into two F-nodes in the middle.
     * @param[in] div_fnode F-node to be divided
     * @return L-node whose weight is updated
     */
    lnode_type* divide_fnode(fnode_type* div_fnode) {
        DEBUG_ABORT_IF(is_head_fnode(div_fnode));
        DEBUG_ABORT_IF_LT(div_fnode->get_data().get_num_runs(), GROUP_BOUND);

        DEBUG_ABORT_IF(!m_lnodes_buffer.is_empty());
        DEBUG_ABORT_IF(!m_fnodes_buffer.is_empty());

        DEBUG_PRINT(tfm::printfln("\n** divide_fnode **");)
        DEBUG_PRINT(tfm::printfln("div_fnode=%s", get_pc(div_fnode));)

        // Keep the original $-mark F-position for Phase E
        const fcursor_type em_fcrsr = m_em_fnode->get_data().get_fcursor(m_em_frnid);

        //
        // Phase C
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase C ==");)

        // TODO: Simplification
        get_overlapped_lnodes(div_fnode, 0, get_sum_fexps(div_fnode), m_lnodes_buffer);
        remove_lnode_links(m_lnodes_buffer);

        //
        // Phase D
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase D ==");)

        // div_fnode is divided into (div_fnode, new_fnode)
        fnode_type* new_fnode = m_findex.divide_after(div_fnode);

        DEBUG_PRINT(tfm::printfln("div_fnode=%s", get_fnode_str(div_fnode));)
        DEBUG_PRINT(tfm::printfln("new_fnode=%s", get_fnode_str(new_fnode));)

        //
        // Phase E
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase E ==");)

        if (m_em_fnode == div_fnode) {
            const size_type em_frpos = em_fcrsr.get_frpos();
            const size_type num_runs = div_fnode->get_data().get_num_runs();
            if (em_frpos < num_runs) {  // $-run is in the before part?
                m_em_fnode = div_fnode;
                m_em_frnid = m_em_fnode->get_data().get_fcursor_from_position(em_frpos).get_frnid();
            } else {  // $-run is in the after part?
                m_em_fnode = new_fnode;
                m_em_frnid = m_em_fnode->get_data().get_fcursor_from_position(em_frpos - num_runs).get_frnid();
            }
        }

        //
        // Phase F
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase F ==");)

        // In this call, div_fnode and new_fnode will be linked appropriately.
        reset_lnode_links_in_order(m_lnodes_buffer);
        // Release
        m_lnodes_buffer.clear();

        {
            const fcover_type fcover = get_fcover(div_fnode);
            div_fnode->get_data().set_weight(fcover.fwght);
        }
        {
            const fcover_type fcover = get_fcover(new_fnode);
            new_fnode->get_data().set_weight(fcover.fwght);
        }

        lnode_type* updated_lnode = std::get<0>(get_first_overlapped_lnode(new_fnode));
        updated_lnode->get_data().add_weight(1);  // for new_fnode

        return updated_lnode;
    }

    /* * * * * * * * * * * * * * * * * * *
     *
     *      Node & Cursor handlers
     *
     * * * * * * * * * * * * * * * * * * */

    //! Check if the given L-node is dummy.
    inline bool is_dmmy_lnode(const lnode_type* lnode) const {
        return lnode == m_lindex.get_head();
    }

    //! Check if the given F-node is head.
    inline bool is_head_fnode(const fnode_type* fnode) const {
        return fnode == m_findex.get_head();
    }

    //! Get the head L-node.
    inline lnode_type* get_head_lnode() const {
        return m_lindex.get_head()->get_next();
    }

    //! Get the tail L-node.
    inline lnode_type* get_tail_lnode() const {
        return m_lindex.get_head()->get_prev();
    }

    //! Get the head F-node.
    inline fnode_type* get_head_fnode() const {
        return m_findex.get_head();
    }

    //! Get the tail F-node.
    inline fnode_type* get_tail_fnode() const {
        return m_findex.get_head()->get_prev();
    }

    //! Check if the given L-run has $'s iterator.
    inline bool is_em_lrun(const lnode_type* lnode, const size_type lrnid) const {
        return (lnode == m_em_lnode) and (lrnid == m_em_lrnid);
    }

    //! Check if the given L-run has $'s iterator.
    inline bool is_em_lrun(const lnode_type* lnode, const lcursor_type& lcrsr) const {
        return (lnode == m_em_lnode) and (lcrsr.get_lrnid() == m_em_lrnid);
    }

    //! Check if the given F-run has $'s iterator.
    inline bool is_em_frun(const fnode_type* fnode, const size_type frnid) const {
        return fnode == m_em_fnode and frnid == m_em_frnid;
    }

    //! Check if the given F-run has $'s iterator.
    inline bool is_em_frun(const fnode_type* fnode, const fcursor_type& fcrsr) const {
        return fnode == m_em_fnode and fcrsr.get_frnid() == m_em_frnid;
    }

    //! Get the exponent of the L-run, considering $-mark.
    inline offset_type get_lexp(const lnode_type* lnode, const size_type lrnid) const {
        const offset_type lexp = lnode->get_data().get_exp(lrnid);
        return is_em_lrun(lnode, lrnid) ? lexp + 1 : lexp;
    }

    //! Get the exponent of the L-run, considering $-mark.
    inline offset_type get_lexp(const lnode_type* lnode, const lcursor_type& lcrsr) const {
        const offset_type lexp = lnode->get_data().get_exp(lcrsr);
        return is_em_lrun(lnode, lcrsr) ? lexp + 1 : lexp;
    }

    //! Get the exponent of the F-run.
    inline offset_type get_fexp(const fnode_type* fnode, const size_type frnid) const {
        return fnode->get_data().get_exp(frnid);
    }

    //! Get the exponent of the F-run.
    inline offset_type get_fexp(const fnode_type* fnode, const fcursor_type& fcrsr) const {
        return fnode->get_data().get_exp(fcrsr);
    }

    //! Get the sum of exponents of the L-group, considering $-marker.
    inline offset_type get_sum_lexps(const lnode_type* lnode) const {
        const offset_type lsize = lnode->get_data().get_sum_exps();
        return lnode == m_em_lnode ? lsize + 1 : lsize;
    }

    //! Get the sum of exponents of the F-group.
    inline offset_type get_sum_fexps(const fnode_type* fnode) const {
        return fnode->get_data().get_sum_exps();
    }

    //! Get the L-run corresponding to the given F-run.
    inline std::tuple<lnode_type*, lcursor_type> get_corr_lcursor(const fnode_type* fnode,
                                                                  const fcursor_type& fcrsr) const {
        const auto [lnode, lrnid] = fnode->get_data().get_lrptr(fcrsr);
        return {lnode, lnode->get_data().get_lcursor(lrnid)};
    }

    //! Get the F-run corresponding to the given L-run.
    inline std::tuple<fnode_type*, fcursor_type> get_corr_fcursor(const lnode_type* lnode,
                                                                  const lcursor_type& lcrsr) const {
        const auto [fnode, frnid] = lnode->get_data().get_frptr(lcrsr);
        return {fnode, fnode->get_data().get_fcursor(frnid)};
    }

    //! Add the weight of the given L-node.
    inline void add_weight(lnode_type* lnode, size_type wght) const {
        lnode->get_data().add_weight(wght);
    }

    //! Sub the weight of the given L-node.
    inline void sub_weight(lnode_type* lnode, size_type wght) const {
        lnode->get_data().sub_weight(wght);
    }

    //! Add the weight of the given F-node.
    inline void add_weight(fnode_type* fnode, size_type wght) const {
        fnode->get_data().add_weight(wght);
    }

    //! Sub the weight of the given F-node.
    inline void sub_weight(fnode_type* fnode, size_type wght) const {
        fnode->get_data().sub_weight(wght);
    }

    //! The distance between the begin of the given F-node and that of the linked L-node.
    inline offset_type get_distance_for_begins(const fnode_type* fnode) const {
        DEBUG_ABORT_IF(!fnode->get_data().get_tnode());
        DEBUG_ABORT_IF_LT(0, fnode->get_data().get_tofst() - get_sum_fexps(fnode) + 1);
        return -1 * (fnode->get_data().get_tofst() - get_sum_fexps(fnode) + 1);
    }

    //! The distance between the end of the given F-node and that of the linked L-node.
    inline offset_type get_distance_for_ends(const fnode_type* fnode) const {
        DEBUG_ABORT_IF(!fnode->get_data().get_tnode());
        const lnode_type* lnode = fnode->get_data().get_tnode();
        DEBUG_ABORT_IF_LT(0, lnode->get_data().get_hofst() - get_sum_lexps(lnode) + 1);
        return -1 * (lnode->get_data().get_hofst() - get_sum_lexps(lnode) + 1);
    }

    /* * * * * * * * * * * * * * * * * * *
     *
     *          Iterator tools
     *
     * * * * * * * * * * * * * * * * * * */

    //! Get the character that the iterator indicates.
    inline uchar_type get_chr(const literator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        return iter.lnode->get_data().get_chr(iter.lcrsr);
    }

    //! Get the character that the iterator indicates.
    inline uchar_type get_chr(const fiterator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        return iter.fnode->get_data().get_chr(iter.fcrsr);
    }

    //! Increment the exponent that the iterator indicates.
    inline void add_exp(const literator_type& iter, size_type exp) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        iter.lnode->get_data().add_exp(iter.lcrsr, exp);
    }

    //! Increment the exponent that the iterator indicates.
    inline void add_exp(const fiterator_type& iter, size_type exp) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        iter.fnode->get_data().add_exp(iter.fcrsr, exp);
    }

    //! Get the SA-entry that the iterator indicates.
    inline size_type get_sae(const literator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        return iter.lnode->get_data().get_sae(iter.lcrsr);
    }

    //! Set the SA-entry that the iterator indicates.
    inline void set_sae(const literator_type& iter, size_type sae) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        iter.lnode->get_data().set_sae(iter.lcrsr, sae);
    }

    //! Reset $'s member variables of L.
    inline void reset_em_lvars() {
        m_em_lnode = nullptr;
        m_em_lrnid = MAX_SIZE_INT;
        m_em_lofst = -1;
    }

    //! Reset $'s member variables of F.
    inline void reset_em_fvars() {
        m_em_fnode = nullptr;
        m_em_frnid = MAX_SIZE_INT;
        m_em_fofst = -1;
    }

    //! Set $'s member variables of L.
    inline void set_em_lvars(const literator_type& iter) {
        m_em_lnode = iter.lnode;
        m_em_lrnid = iter.lcrsr.get_lrnid();
        m_em_lofst = iter.lofst;
    }

    //! Set $'s member variables of F.
    inline void set_em_fvars(const fiterator_type& iter) {
        m_em_fnode = iter.fnode;
        m_em_frnid = iter.fcrsr.get_frnid();
        m_em_fofst = iter.fofst;
    }

    //! Get $'s iterator of L.
    inline literator_type get_em_literator() const {
        return {m_em_lnode, m_em_lnode->get_data().get_lcursor(m_em_lrnid), m_em_lofst};
    }

    //! Get $'s iterator of F.
    inline fiterator_type get_em_fiterator() const {
        return {m_em_fnode, m_em_fnode->get_data().get_fcursor(m_em_frnid), m_em_fofst};
    }

    //! Check if the cursor indicates $'s run.
    inline bool is_em_lcursor(const literator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        return (iter.lnode == m_em_lnode) and (iter.lcrsr.get_lrnid() == m_em_lrnid);
    }

    //! Check if the cursor indicates $'s run.
    inline bool is_em_fcursor(const fiterator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        return iter.fnode == m_em_fnode and iter.fcrsr.get_frnid() == m_em_frnid;
    }

    //! Check if the iterator indicates $'s marker.
    inline bool is_em_literator(const literator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        DEBUG_ABORT_IF(!iter.is_valid_offset());
        return is_em_lcursor(iter) and iter.lofst == m_em_lofst;
    }

    //! Check if the iterator indicates $'s marker.
    inline bool is_em_fiterator(const fiterator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        DEBUG_ABORT_IF(!iter.is_valid_offset());
        return is_em_fcursor(iter) and iter.fofst == m_em_fofst;
    }

    //!
    inline literator_type get_last_lcursor(lnode_type* lnode) const {
        return {lnode, lnode->get_data().get_last_lcursor()};
    }

    //!
    inline fiterator_type get_last_fcursor(fnode_type* fnode) const {
        return {fnode, fnode->get_data().get_last_fcursor()};
    }

    //! Get the L-iterator of the head.
    inline literator_type get_head_literator() const {
        lnode_type* lnode = get_head_lnode();
        return {lnode, lnode->get_data().get_first_lcursor(), 0};
    }

    //! Get the F-iterator of the head.
    inline fiterator_type get_head_fiterator() const {
        fnode_type* fnode = get_head_fnode();
        return {fnode, fnode->get_data().get_first_fcursor(), 0};
    }

    //! Get the L-iterator of the tail.
    inline literator_type get_tail_literator() const {
        // auto [lnode, lcrsr] = get_tail_lcursor();
        lnode_type* lnode = get_tail_lnode();
        lcursor_type lcrsr = lnode->get_data().get_last_lcursor();
        return {lnode, lcrsr, get_lexp(lnode, lcrsr) - 1};
    }

    //! Get the F-iterator of the tail.
    inline fiterator_type get_tail_fiterator() const {
        // auto [fnode, fcrsr] = get_tail_fcursor();
        fnode_type* fnode = get_tail_fnode();
        fcursor_type fcrsr = fnode->get_data().get_last_fcursor();
        return {fnode, fcrsr, get_fexp(fnode, fcrsr) - 1};
    }

    //! Check if the L-cursor is the head in L.
    inline bool is_head_lcursor(const literator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        const lnode_type* lhead = get_head_lnode();
        if (lhead != iter.lnode) {
            return false;
        }
        return lhead->get_data().is_first_lcursor(iter.lcrsr);
    }

    //! Check if the F-cursor is the head in F.
    inline bool is_head_fcursor(const fiterator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        const fnode_type* fhead = get_head_fnode();
        if (fhead != iter.fnode) {
            return false;
        }
        return fhead->get_data().is_first_fcursor(iter.fcrsr);
    }

    //! Get the previous L-cursor.
    inline literator_type get_prev_lcursor(const literator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        auto tmp = iter;
        set_prev_lcursor(tmp);
        return tmp;
    }

    //! Get the next L-cursor.
    inline literator_type get_next_lcursor(const literator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        auto tmp = iter;
        set_next_lcursor(tmp);
        return tmp;
    }

    //! Get the previous F-cursor.
    inline fiterator_type get_prev_fcursor(const fiterator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        auto tmp = iter;
        set_prev_fcursor(tmp);
        return tmp;
    }

    //! Get the next F-cursor.
    inline fiterator_type get_next_fcursor(const fiterator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        auto tmp = iter;
        set_next_fcursor(tmp);
        return tmp;
    }

    //! Get the previous L-cursor.
    inline void set_prev_lcursor(literator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        if (iter.lnode->get_data().is_first_lcursor(iter.lcrsr)) {
            if (iter.lnode == get_head_lnode()) {
                iter = literator_type{};  // set invalid
            } else {
                iter.lnode = iter.lnode->get_prev();
                iter.lcrsr = iter.lnode->get_data().get_last_lcursor();
            }
        } else {
            iter.lcrsr = iter.lnode->get_data().get_prev_lcursor(iter.lcrsr);
        }
    }

    //! Get the next L-cursor.
    inline void set_next_lcursor(literator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        if (iter.lnode->get_data().is_last_lcursor(iter.lcrsr)) {
            if (iter.lnode == get_tail_lnode()) {
                iter = literator_type{};  // set invalid
            } else {
                iter.lnode = iter.lnode->get_next();
                iter.lcrsr = iter.lnode->get_data().get_first_lcursor();
            }
        } else {
            iter.lcrsr = iter.lnode->get_data().get_next_lcursor(iter.lcrsr);
        }
    }

    //! Get the previous F-cursor.
    inline void set_prev_fcursor(fiterator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        if (iter.fnode->get_data().is_first_fcursor(iter.fcrsr)) {
            if (iter.fnode == get_head_fnode()) {
                iter = fiterator_type{};  // set invalid
            } else {
                iter.fnode = iter.fnode->get_prev();
                iter.fcrsr = iter.fnode->get_data().get_last_fcursor();
            }
        } else {
            iter.fcrsr = iter.fnode->get_data().get_prev_fcursor(iter.fcrsr);
        }
    }

    //! Get the next F-cursor.
    inline void set_next_fcursor(fiterator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        if (iter.fnode->get_data().is_last_fcursor(iter.fcrsr)) {
            if (iter.fnode == get_tail_fnode()) {
                iter = fiterator_type{};  // set invalid
            } else {
                iter.fnode = iter.fnode->get_next();
                iter.fcrsr = iter.fnode->get_data().get_first_fcursor();
            }
        } else {
            iter.fcrsr = iter.fnode->get_data().get_next_fcursor(iter.fcrsr);
        }
    }

    inline void set_next_literator(literator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        if (is_last_char_in_lrun(iter)) {
            set_next_lcursor(iter);
            iter.lofst = 0;
        } else {
            iter.lofst += 1;
        }
    }

    inline void set_next_fiterator(fiterator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        if (is_last_char_in_run(iter)) {
            set_next_fcursor(iter);
            iter.lofst = 0;
        } else {
            iter.lofst += 1;
        }
    }

    //! Check if the given offset indicates the first char of the given L-run, considering $-marker.
    inline bool is_first_char_in_lrun(const literator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        return iter.lofst == 0;
    }

    //! Check if the given offset indicates the last char of the given L-run, considering $-marker.
    inline bool is_last_char_in_lrun(const literator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        return get_lexp(iter.lnode, iter.lcrsr) == iter.lofst + 1;
    }

    //! Check if the given offset indicates the first char of the given F-run.
    inline bool is_first_char_in_frun(const fiterator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        return iter.fofst == 0;
    }

    //! Check if the given offset indicates the last char of the given F-run.
    inline bool is_last_char_in_frun(const fiterator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        return get_fexp(iter.fnode, iter.fcrsr) == iter.fofst + 1;
    }

    //! Check if the given offset indicates the first char of the given L-node, considering $-marker.
    inline bool is_first_char_in_lnode(const literator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        return iter.lofst == 0 and iter.lnode->get_data().is_first_lcursor(iter.lcrsr);
    }

    //! Check if the given offset indicates the last char of the given L-node, considering $-marker.
    inline bool is_last_char_in_lnode(const literator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        return iter.lnode->get_data().is_last_lcursor(iter.lcrsr) and
               get_lexp(iter.lnode, iter.lcrsr) == iter.lofst + 1;
    }

    //! Check if the given offset indicates the first char of the given F-node.
    inline bool is_first_char_in_fnode(const fiterator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        return iter.fofst == 0 and iter.fnode->get_data().is_first_fcursor(iter.fcrsr);
    }

    //! Check if the given offset indicates the last char of the given F-node.
    inline bool is_last_char_in_fnode(const fiterator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        return iter.fnode->get_data().is_last_fcursor(iter.fcrsr) and
               get_fexp(iter.fnode, iter.fcrsr) == iter.fofst + 1;
    }

    inline bool is_first_run_in_lnode(const literator_type& iter) const {
        return iter.lnode->get_data().is_first_lcursor(iter.lcrsr);
    }

    inline bool is_first_run_in_fnode(const fiterator_type& iter) const {
        return iter.fnode->get_data().is_first_fcursor(iter.fcrsr);
    }

    //! Get the L-run corresponding to the given F-run.
    inline literator_type get_corr_lcursor(const fiterator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        const auto [lnode, lrnid] = iter.fnode->get_data().get_lrptr(iter.fcrsr);
        return {lnode, lnode->get_data().get_lcursor(lrnid)};
    }

    //! Get the F-run corresponding to the given L-run.
    inline fiterator_type get_corr_fcursor(const literator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        const auto [fnode, frnid] = iter.lnode->get_data().get_frptr(iter.lcrsr);
        return {fnode, fnode->get_data().get_fcursor(frnid)};
    }

    //! Get the L-run corresponding to the given F-run.
    // inline literator_type get_corr_literator(const fiterator_type& fiter) const {
    //     DEBUG_ABORT_IF(!fiter.is_valid());
    //     const auto [lnode, lrnid] = fiter.fnode->get_data().get_lrptr(fiter.fcrsr);
    //     return {lnode, lnode->get_data().get_lcursor(lrnid), fiter.fofst};
    // }

    //! Get the F-run corresponding to the given L-run.
    inline fiterator_type get_corr_fiterator(const literator_type& liter) const {
        DEBUG_ABORT_IF(!liter.is_valid());
        auto [fnode, frnid] = liter.lnode->get_data().get_frptr(liter.lcrsr);
        if (is_em_lcursor(liter)) {
            if (m_em_lofst < liter.lofst) {
                return {fnode, fnode->get_data().get_fcursor(frnid), liter.lofst - 1};
            } else if (m_em_lofst == liter.lofst) {
                return get_head_fiterator();
            }
        }
        return {fnode, fnode->get_data().get_fcursor(frnid), liter.lofst};
    }

    //! Get the exponent of the given L-run, considering $-mark.
    inline offset_type get_lexp(const literator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        const offset_type lexp = iter.lnode->get_data().get_exp(iter.lcrsr);
        return is_em_lrun(iter.lnode, iter.lcrsr) ? lexp + 1 : lexp;
    }

    //! Get the exponent of the given F-run.
    inline offset_type get_fexp(const fiterator_type& iter) const {
        DEBUG_ABORT_IF(!iter.is_valid());
        return iter.fnode->get_data().get_exp(iter.fcrsr);
    }

    //! Insert a new L-run before the given L-run, and update the cursor of the iterator.
    inline void insert_lrun_before(literator_type& liter, uchar_type chr, size_type exp) const {
        DEBUG_ABORT_IF(!liter.is_valid());
        liter.lcrsr = liter.lnode->get_data().insert_before(liter.lcrsr, chr, exp);
    }

    //! Insert a new F-run before the given F-run, and update the cursor of the iterator.
    inline void insert_frun_before(fiterator_type& fiter, const literator_type& liter) const {
        DEBUG_ABORT_IF(!liter.is_valid());
        DEBUG_ABORT_IF(!fiter.is_valid());
        fiter.fcrsr = fiter.fnode->get_data().insert_before(fiter.fcrsr, liter.lnode, liter.lcrsr.get_lrnid());
        liter.lnode->get_data().set_frptr(liter.lcrsr, fiter.fnode, fiter.fcrsr.get_frnid());
    }

    //! Insert a new L-run after the given L-run, and update the cursor of the iterator.
    inline void insert_lrun_after(literator_type& liter, uchar_type chr, size_type exp) const {
        DEBUG_ABORT_IF(!liter.is_valid());
        liter.lcrsr = liter.lnode->get_data().insert_after(liter.lcrsr, chr, exp);
    }

    //! Insert a new F-run after the given F-run, and update the cursor of the iterator.
    inline void insert_frun_after(fiterator_type& fiter, const literator_type& liter) const {
        DEBUG_ABORT_IF(!liter.is_valid());
        DEBUG_ABORT_IF(!fiter.is_valid());
        fiter.fcrsr = fiter.fnode->get_data().insert_after(fiter.fcrsr, liter.lnode, liter.lcrsr.get_lrnid());
        liter.lnode->get_data().set_frptr(liter.lcrsr, fiter.fnode, fiter.fcrsr.get_frnid());
    }

    inline void split_run_after(literator_type& liter, fiterator_type& fiter, size_type exp_befo) const {
        DEBUG_ABORT_IF(!liter.is_valid());
        DEBUG_ABORT_IF(!fiter.is_valid());
        liter.lcrsr = liter.lnode->get_data().split_after(liter.lcrsr, exp_befo);
        fiter.fcrsr = fiter.fnode->get_data().split_after(fiter.fcrsr, liter.lnode, liter.lcrsr.get_lrnid());
        liter.lnode->get_data().set_frptr(liter.lcrsr, fiter.fnode, fiter.fcrsr.get_frnid());
    }

    //! Check if the given L-addresses are equivalent.
    inline bool is_equal_literator(const literator_type& a, const literator_type& b) const {
        return (a.lnode == b.lnode) and (a.lcrsr.get_lrnid() == b.lcrsr.get_lrnid()) and (a.lofst == b.lofst);
    }

    //! Check if the given F-addresses are equivalent.
    inline bool is_equal_fiterator(const fiterator_type& a, const fiterator_type& b) const {
        return (a.fnode == b.fnode) and (a.fcrsr.get_frnid() == b.fcrsr.get_frnid()) and (a.fofst == b.fofst);
    }

    inline bool compare_order(const literator_type& litr1, const literator_type& litr2) const {
        // First, compare with the L-node order
        const loint_type order1 = lo_common::get_basic_order(litr1.lnode);
        const loint_type order2 = lo_common::get_basic_order(litr2.lnode);
        if (order1 != order2) {
            return order1 < order2;
        }
        // Second, compare with the offset position in the L-node
        if (litr1.lcrsr.get_lrpos() != litr2.lcrsr.get_lrpos()) {
            return litr1.lcrsr.get_lrpos() < litr2.lcrsr.get_lrpos();
        }
        return litr1.lofst < litr2.lofst;
    }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
     *
     *      Backward searchers (in log time)
     *
     * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    //! Search the most backword F-node and run in the node such that
    //!  - F-run's character is no greater than 'q_chr',
    //!  - F-run's order is no greater than 'lo_common::get_basic_order(q_lnode)', and
    //!  - F-run's position in the corresponding L-node is no greater than 'q_lnode->get_data().get_lrpos(q_lrnid)'.
    //! Return its pointer and run-ID.
    inline fiterator_type predecessor(uchar_type chr, literator_type liter) const {
        const loint_type q_order = lo_common::get_basic_order(liter.lnode);
        const size_type q_second_order = liter.lcrsr.get_lrpos();

        // TODO: If the group contains the same character, this can be used to omit binary search.
        return {m_findex.predecessor(chr, q_order, q_second_order)};
    }

    inline bwres_etype predecessor_in_bwsearch(uchar_type chr, const literator_type& liter,
                                               fiterator_type& fiter) const {
        DEBUG_ABORT_IF_EQ(chr, END_MARKER);  // Search for END_MARKER is not supported.
        DEBUG_ABORT_IF(m_em_lnode == nullptr);
        DEBUG_ABORT_IF(m_em_fnode == nullptr);

        if (chr == get_chr(liter)) {
            // The LF-mapped F-node is the linked one.
            fiter = get_corr_fcursor(liter);

            // SPECIAL CASE: The address indicates $-marker
            if (is_em_literator(liter)) {
                if (is_first_char_in_lrun(get_em_literator())) {
                    // Visit the previous L-node instead of the decrement
                    set_prev_fcursor(fiter);
                    if (chr != get_chr(fiter)) {
                        return bwres_etype::FAILED;
                    }
                    fiter.fofst = get_fexp(fiter) - 1;
                } else {
                    fiter.fofst = liter.lofst - 1;  // due to the virtual $-marker
                }
                return bwres_etype::NOT_MATCHED;
            }
            // SPECIAL CASE: The L-node has $-marker (but, not indicating)
            else if (is_em_lcursor(liter) and liter.lofst != m_em_lofst) {
                if (liter.lofst < m_em_lofst) {
                    fiter.fofst = liter.lofst;
                } else {
                    fiter.fofst = liter.lofst - 1;  // due to the virtual $-marker
                }
                return bwres_etype::MATCHED;
            }
            // NORMAL CASE
            fiter.fofst = liter.lofst;
            return bwres_etype::MATCHED;
        } else {
            const loint_type q_order = lo_common::get_basic_order(liter.lnode);
            const size_type q_second_order = liter.lcrsr.get_lrpos();

            auto [fnode, fcrsr] = m_findex.exact_predecessor(chr, q_order, q_second_order);
            if (fnode == nullptr) {
                return bwres_etype::FAILED;
            }
            fiter = {fnode, fcrsr, get_fexp(fnode, fcrsr) - 1};
            return bwres_etype::NOT_MATCHED;
        }
    }

    inline bwres_etype successor_in_bwsearch(uchar_type chr, const literator_type& liter, fiterator_type& fiter) const {
        DEBUG_ABORT_IF_EQ(chr, END_MARKER);  // Search for END_MARKER is not supported.
        DEBUG_ABORT_IF(m_em_lnode == nullptr);
        DEBUG_ABORT_IF(m_em_fnode == nullptr);

        if (chr == get_chr(liter)) {
            // The LF-mapped F-run is the linked one.
            fiter = get_corr_fcursor(liter);

            // SPECIAL CASE: The address indicates $-marker
            if (is_em_literator(liter)) {
                if (is_last_char_in_lrun(get_em_literator())) {
                    // This case will arise only when $-marker is placed at the tail of L.
                    return bwres_etype::FAILED;
                }
                fiter.fofst = liter.lofst;
                return bwres_etype::NOT_MATCHED_BY_EM;
            }
            // SPECIAL CASE: The L-node has $-marker (but, not indicating)
            // liter.lofst != m_em_lofst may not need.
            else if (is_em_lcursor(liter) and liter.lofst != m_em_lofst) {
                if (liter.lofst < m_em_lofst) {
                    fiter.fofst = liter.lofst;
                } else {
                    fiter.fofst = liter.lofst - 1;  // due to the virtual $-marker
                }
                return bwres_etype::MATCHED;
            }
            // NORMAL CASE
            fiter.fofst = liter.lofst;
            return bwres_etype::MATCHED;
        } else {
            const loint_type q_order = lo_common::get_basic_order(liter.lnode);
            const size_type q_second_order = liter.lcrsr.get_lrpos();

            auto [fnode, fcrsr] = m_findex.exact_successor(chr, q_order, q_second_order);
            if (fnode == nullptr) {
                return bwres_etype::FAILED;
            }

            fiter = {fnode, fcrsr, 0};
            return bwres_etype::NOT_MATCHED;
        }
    }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
     *
     *      Getters for linkages (in constant time)
     *
     * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    inline literator_type get_overlapped_literator(const fiterator_type& fiter) const {
        DEBUG_ABORT_IF(!fiter.fnode);
        DEBUG_ABORT_IF(is_head_fnode(fiter.fnode));

        auto [q_fnode, q_fcrsr, q_fofst] = fiter;

        // Exception case: \Lambda_t indicates the tail character of F.
        {
            auto [tail_fnode, tail_fcrsr, tail_fexp] = get_tail_fiterator();
            if (q_fnode == tail_fnode and q_fcrsr == tail_fcrsr and tail_fexp == q_fofst) {
                auto [tail_lnode, tail_lcrsr, tail_lexp] = get_tail_literator();
                if (m_em_lnode != nullptr) {
                    return {tail_lnode, tail_lcrsr, tail_lexp};
                }
                return {tail_lnode, tail_lcrsr, tail_lexp + 1};
            }
        }

        auto [lnode, fofst] = get_first_overlapped_lnode(q_fnode);
        const offset_type q_fofst_in_node = q_fnode->get_data().get_offset(q_fcrsr) + q_fofst;

        while (fofst < q_fofst_in_node) {
            lnode = lnode->get_next();
            fofst += get_sum_lexps(lnode);
        }

        lcursor_type lcrsr = {};
        offset_type lofst = (get_sum_lexps(lnode) - 1) - (fofst - q_fofst_in_node);

        if (lnode == m_em_lnode) {
            const offset_type em_lofst_in_node = m_em_lnode->get_data().get_offset(m_em_lrnid) + m_em_lofst;
            if (em_lofst_in_node < lofst) {
                // Then, lofst will not indicate $ since em_lofst_in_node < lofst.
                std::tie(lcrsr, lofst) = lnode->get_data().access_with_offset(lofst - 1);
                if (is_em_lrun(lnode, lcrsr) and m_em_lofst <= lofst) {
                    lofst += 1;  // for $
                }
            } else {
                std::tie(lcrsr, lofst) = lnode->get_data().access_with_offset(lofst);
            }
        } else {
            std::tie(lcrsr, lofst) = lnode->get_data().access_with_offset(lofst);
        }

        return {lnode, lcrsr, lofst};
    }

    inline fiterator_type get_overlapped_fiterator(const literator_type& liter) const {
        DEBUG_ABORT_IF(!liter.lnode);
        DEBUG_ABORT_IF(is_dmmy_lnode(liter.lnode));

        auto [q_lnode, q_lcrsr, q_lofst] = liter;
        auto [fnode, lofst] = get_first_overlapped_fnode(q_lnode);
        offset_type q_lofst_in_node = q_lnode->get_data().get_offset(q_lcrsr) + q_lofst;

        // Solution for that get_offset() does not consider $-marker
        if (q_lnode == m_em_lnode) {
            const lcursor_type em_lcrsr = m_em_lnode->get_data().get_lcursor(m_em_lrnid);
            if (em_lcrsr.get_lrpos() < q_lcrsr.get_lrpos()) {
                q_lofst_in_node += 1;
            }
        }

        while (lofst < q_lofst_in_node) {
            fnode = fnode->get_next();
            lofst += get_sum_fexps(fnode);
        }

        fcursor_type fcrsr = {};
        offset_type fofst = (get_sum_fexps(fnode) - 1) - (lofst - q_lofst_in_node);
        std::tie(fcrsr, fofst) = fnode->get_data().access_with_offset(fofst);

        return {fnode, fcrsr, fofst};
    }

    //!  Return (lnode, fofst) such that fnode[fofst] corresponds to the back of lnode.
    //!
    //!    e.g.,
    //!               [F]      [L]
    //!            v1  a
    //!                a        b  u1
    //!            v2  a        a
    //!                a        c  u2 *      The output is (u2, 1)
    //!       (in) v3  b        c            since v3[1] corresponds to u2.back().
    //!                b  --->  b
    //!                b        a  u3
    //!                         b
    //!
    inline std::tuple<lnode_type*, offset_type> get_first_overlapped_lnode(const fnode_type* fnode) const {
        if (is_head_fnode(fnode)) {
            lnode_type* lnode = get_head_lnode();
            return {lnode, get_sum_lexps(lnode) - 1};
        }

        // Distance from the head of the given F-node
        offset_type fofst = 0;

        // Move to the closest F-node with a tail link
        // (This will be done since the head nodes are always linked)
        fnode = fnode->get_prev();  // start at the previous F-node
        while (!fnode->get_data().get_tnode()) {
            fofst -= get_sum_fexps(fnode);
            fnode = fnode->get_prev();
        }
        fofst -= get_sum_fexps(fnode);

        // Move to the tail-linked L-node
        lnode_type* lnode = fnode->get_data().get_tnode();
        fofst += get_distance_for_begins(fnode);

        // Move to the first overlapped L-node with the given F-node
        do {
            fofst += get_sum_lexps(lnode);
            lnode = lnode->get_next();
        } while (fofst <= 0);

        return {lnode->get_prev(), fofst - 1};
    }

    //!  Return (fnode, lofst) such that lnode[lofst] corresponds to the back of fnode.
    //!
    //!     e.g.,
    //!               [F]      [L]
    //!            v1  a
    //!                a        b  u1
    //!          * v2  a        a
    //!                a        c  u2 (in)   The output is (v2, 2)
    //!                a        c            since u2[2] corresponds to v2.back.
    //!                a  <---  b
    //!            v3  b        a  u3
    //!                b        b
    //!
    inline std::tuple<fnode_type*, offset_type> get_first_overlapped_fnode(const lnode_type* lnode) const {
        DEBUG_ABORT_IF(is_dmmy_lnode(lnode));

        // Distance from the head of the given L-node
        offset_type lofst = 0;

        // Move to the closest L-node with a head link
        // (This will be always done since the head L-node has the head link to the head F-node)
        while (!lnode->get_data().get_hnode()) {
            lnode = lnode->get_prev();
            lofst -= get_sum_lexps(lnode);
        }

        // Move to the tail-linked L-node
        fnode_type* fnode = lnode->get_data().get_hnode();
        lofst -= get_distance_for_begins(fnode);

        // Move to the first overlapped F-node with the given L-node
        while (true) {
            lofst += get_sum_fexps(fnode);
            if (lofst > 0) {  // overlapped?
                break;
            }
            fnode = fnode->get_next();
        }

        return {fnode, lofst - 1};
    }

    //!  Return the L-nodes overlapped with the F-interval of size 'fsize' starting at the 'fofst'-th of 'fnode'.
    //!
    //!     e.g.,
    //!                     [F]      [L]
    //!                  v1  a
    //!                      a        b  u1
    //!             (in) v2  a        a
    //!     fofst = 1 -->    a +    + c  u2 *       The output is [u2, u3].
    //!                  v3  b |    | c
    //!                      b |    | b
    //!     fsize = 4 -->    b +    | a  u3 *
    //!                             + b
    //!
    inline void get_overlapped_lnodes(const fnode_type* fnode, offset_type fofst, const offset_type fsize,
                                      lbuffer_type& buf) {
        DEBUG_ABORT_IF(is_head_fnode(fnode));

        if (fsize <= 0) {
            return;
        }

        const offset_type fofst_b = fofst;
        const offset_type fofst_e = fofst + fsize - 1;

        lnode_type* lnode = nullptr;
        std::tie(lnode, fofst) = get_first_overlapped_lnode(fnode);

        // non-overlapped?
        while (fofst < fofst_b) {
            lnode = lnode->get_next();
            fofst += get_sum_lexps(lnode);
        }

        // overlapped?
        while (!is_dmmy_lnode(lnode)) {
            buf.push_back(lnode);
            if (fofst_e <= fofst) {
                break;
            }
            lnode = lnode->get_next();
            fofst += get_sum_lexps(lnode);
        }
    }

    //!  Return the F-nodes overlapped with the L-interval of size 'lsize' starting at the 'lofst'-th of 'lnode'.
    //!
    //!     e.g.,
    //!                 [F]      [L]
    //!              v1  a                          The output is [v2, v3].
    //!                  a        b  u1 (in)
    //!            * v2  a +      a
    //!                  a |    + c  u2  <-- lofst = 2
    //!            * v3  b |    | c
    //!                  b |    + b      <-- lsize = 3
    //!                  b +      a  u3
    //!                           b
    //!
    inline void get_overlapped_fnodes(const lnode_type* lnode, offset_type lofst, const offset_type lsize,
                                      fbuffer_type& buf) {
        DEBUG_ABORT_IF(is_dmmy_lnode(lnode));

        if (lsize <= 0) {
            return;
        }

        const offset_type lofst_b = lofst;
        const offset_type lofst_e = lofst + lsize - 1;

        fnode_type* fnode = nullptr;
        std::tie(fnode, lofst) = get_first_overlapped_fnode(lnode);

        // non-overlapped?
        while (lofst < lofst_b) {
            fnode = fnode->get_next();
            lofst += get_sum_fexps(fnode);
        }

        // overlapped?
        while (true) {
            buf.push_back(fnode);
            if (lofst_e <= lofst) {
                break;
            }
            fnode = fnode->get_next();
            lofst += get_sum_fexps(fnode);
            if (is_head_fnode(fnode)) {  // cycled?
                break;
            }
        }
    }

    //!  Return lcover_type(lnode, fnode, lofst, weight) such that
    //!   - 'lnode' is the input L-node,
    //!   - 'fnode' is the first F-node whose head is covered by 'lnode',
    //!   - 'lofst' means that 'lnode[lofst]' corresponds to the back of 'fnode', and
    //!   - 'weight' means that 'weight's F-nodes from 'fnode' are covered by 'lnode'.
    //!
    //!     e.g.,
    //!                     [F]      [L]
    //!                  v1  a
    //!                      a        b  u1 (in)
    //!                * v2  a        b
    //!                      a  <---  b           The output is (u1, v2, 2, 3)
    //!                * v3  b        b           since u1[2] corresponds to v2.back.
    //!                      b        b
    //!                * v4  c        b
    //!                      c
    //!
    inline lcover_type get_lcover(const lnode_type* lnode) const {
        auto [fnode, lofst_e] = get_first_overlapped_fnode(lnode);
        offset_type lofst_b = lofst_e - get_sum_fexps(fnode) + 1;

        // Is the first overlapped F-node not covered? (i.e., the head is not overlapped)
        if (lofst_b < 0) {
            fnode = fnode->get_next();
            lofst_b = lofst_e + 1;
            lofst_e = lofst_b + get_sum_fexps(fnode) - 1;
        }

        lcover_type ret = {const_cast<lnode_type*>(lnode), fnode, lofst_e, 0};
        const offset_type lexp = get_sum_lexps(lnode);

        while (lofst_b < lexp) {
            ret.lwght += 1;
            lofst_b += get_sum_fexps(fnode);
            fnode = fnode->get_next();
            if (is_head_fnode(fnode)) {  // cycled?
                break;
            }
        }
        return ret;
    }

    //!  Return fcover_type(fnode, lnode, fofst, weight) such that
    //!   - 'fnode' is the input F-node,
    //!   - 'lnode' is the first L-node whose head is covered by 'fnode',
    //!   - 'fofst' means that 'fnode[fofst]' corresponds to the back of 'lnode', and
    //!   - 'weight' means that 'weight's L-nodes from 'lnode' are covered by 'fnode'.
    //!
    //!     e.g.,
    //!                     [F]      [L]
    //!                               b  u1
    //!             (in) v1  a        b
    //!                      a        a  u2 *
    //!                      a        a         The output is (v1, u2, 3, 3)
    //!                      a  --->  a           since v1[3] corresponds to u2.back.
    //!                      a        b  u3 *
    //!                      a        c  u4 *
    //!                      a        c
    //!
    inline fcover_type get_fcover(const fnode_type* fnode) const {
        auto [lnode, fofst_e] = get_first_overlapped_lnode(fnode);
        offset_type fofst_b = fofst_e - get_sum_lexps(lnode) + 1;

        // Is the first overlapped L-node not covered? (i.e., the head is not overlapped)
        if (fofst_b < 0) {
            lnode = lnode->get_next();
            fofst_b = fofst_e + 1;
            fofst_e = fofst_b + get_sum_lexps(lnode) - 1;
        }

        fcover_type ret = {const_cast<fnode_type*>(fnode), lnode, fofst_e, 0};
        const offset_type fexp = get_sum_fexps(fnode);

        while (fofst_b < fexp) {
            ret.fwght += 1;
            fofst_b += get_sum_lexps(lnode);
            lnode = lnode->get_next();
            if (is_dmmy_lnode(lnode)) {  // cycled?
                break;
            }
        }
        return ret;
    }

    /* * * * * * * * * * * * * * * * * * *
     *
     *         Update operations
     *
     * * * * * * * * * * * * * * * * * * */

    inline void push_fat_lnode(lnode_type* lnode) {
        if (lnode and lnode->get_data().get_weight() >= DIV_BOUND) {
            m_fat_lnodes.push(lnode);
        }
    }
    inline void push_fat_fnode(fnode_type* fnode) {
        if (fnode and fnode->get_data().get_weight() >= DIV_BOUND) {
            m_fat_fnodes.push(fnode);
        }
    }

    //! Split L-run (lnode, lcrsr_befo) such that the before part has lexp_befo characters.
    //! Note that, if L-run (lnode, lcrsr_befo) has $-mark, lexp_befo should NOT consider $.
    inline std::pair<fnode_type*, fcursor_type> split_lrun(lnode_type* lnode, lcursor_type lcrsr_befo,
                                                           size_type lexp_befo) {
        DEBUG_ABORT_IF_EQ(lexp_befo, get_lexp(lnode, lcrsr_befo));  // not need to split?

        auto [fnode, fcrsr_befo] = get_corr_fcursor(lnode, lcrsr_befo);

        lcursor_type lcrsr_aftr = lnode->get_data().split_after(lcrsr_befo, lexp_befo);
        fcursor_type fcrsr_aftr = fnode->get_data().split_after(fcrsr_befo, lnode, lcrsr_aftr.get_lrnid());
        lnode->get_data().set_frptr(lcrsr_aftr, fnode, fcrsr_aftr.get_frnid());

        // split_after(liter, fiter, lexp_befo);

        if (is_em_lrun(lnode, lcrsr_befo) and lexp_befo <= size_type(m_em_lofst)) {
            // Then, em_lrnid is updated to the after L-run
            m_em_lrnid = lcrsr_aftr.get_lrnid();
            m_em_lofst -= lexp_befo;
        }
        if (is_em_frun(fnode, fcrsr_befo) and lexp_befo <= size_type(m_em_fofst)) {
            // Then, em_frnid is updated to the after F-run
            m_em_frnid = fcrsr_aftr.get_frnid();
            m_em_fofst -= lexp_befo;
        }

        return {fnode, fcrsr_befo};
    }

    //! Split F-run (fnode, fcrsr_befo) such that the front part has fexp_befo characters.
    inline std::pair<lnode_type*, lcursor_type> split_frun(fnode_type* fnode, fcursor_type fcrsr_befo,
                                                           size_type fexp_befo) {
        DEBUG_ABORT_IF_EQ(fexp_befo, get_fexp(fnode, fcrsr_befo));  // not need to split?

        auto [lnode, lcrsr_befo] = get_corr_lcursor(fnode, fcrsr_befo);

        lcursor_type lcrsr_aftr = lnode->get_data().split_after(lcrsr_befo, fexp_befo);
        fcursor_type fcrsr_aftr = fnode->get_data().split_after(fcrsr_befo, lnode, lcrsr_aftr.get_lrnid());
        lnode->get_data().set_frptr(lcrsr_aftr, fnode, fcrsr_aftr.get_frnid());

        if (is_em_lrun(lnode, lcrsr_befo) and fexp_befo <= size_type(m_em_lofst)) {
            // Then, em_lrnid is updated to the after L-run
            m_em_lrnid = lcrsr_aftr.get_lrnid();
            m_em_lofst -= fexp_befo;
        }
        if (is_em_frun(fnode, fcrsr_befo) and fexp_befo <= size_type(m_em_fofst)) {
            // Then, em_frnid is updated to the after F-run
            m_em_frnid = fcrsr_aftr.get_frnid();
            m_em_fofst -= fexp_befo;
        }

        return {lnode, lcrsr_befo};
    }

    // Remove the head link of 'lnode' and the corresponding tail link.
    inline void remove_lnode_link(lnode_type* lnode) {
        fnode_type* fnode = lnode->get_data().get_hnode();
        if (!fnode) {
            return;
        }

        // Keep F-head's link because it is never changed (also the hofst).
        // Note: DON'T remove the link because FL-pair queries assume that the heads are linked.
        if (is_head_fnode(fnode)) {
            return;
        }

        // Remove links
        lnode->get_data().reset_hlink();
        fnode->get_data().reset_tlink();
    }

    // Remove the head links in 'lnodes' and the corresponding tail links.
    inline void remove_lnode_links(lbuffer_type& lnodes) {
        for (size_type i = 0; i < lnodes.get_size(); i++) {
            remove_lnode_link(lnodes[i]);
        }
    }

    // Remove the tail link of 'fnode' and the corresponding head link.
    inline void remove_fnode_link(fnode_type* fnode) {
        // Keep F-head's link because it is never changed (also the hofst).
        // Note: DON'T remove the link because FL-pair queries assume that the heads are linked.
        if (is_head_fnode(fnode)) {
            return;
        }

        lnode_type* lnode = fnode->get_data().get_tnode();
        if (!lnode) {
            return;
        }

        // Remove links
        lnode->get_data().reset_hlink();
        fnode->get_data().reset_tlink();
    }

    // Remove the tail links in 'fnodes' and the corresponding head links.
    inline void remove_fnode_links(fbuffer_type& fnodes) {
        for (size_type i = 0; i < fnodes.get_size(); i++) {
            remove_fnode_link(fnodes[i]);
        }
    }

    inline void reset_lnode_link(lnode_type* lnode) {
        DEBUG_ABORT_IF(is_dmmy_lnode(lnode));

        if (lnode->get_data().get_hnode()) {
            return;
        }

        // Consider 'lnode' as the pivot, the L-interval is [lofst_lb, lofst_le].
        const offset_type lofst_lb = 0;
        const offset_type lofst_le = get_sum_lexps(lnode) - 1;

        // Consider 'lnode' as the pivot, the F-interval is [lofst_fb, lofst_fe].
        auto [fnode, lofst_fe] = get_first_overlapped_fnode(lnode);
        offset_type lofst_fb = lofst_fe - get_sum_fexps(fnode) + 1;

        if ((lofst_fb <= lofst_lb) and (lofst_lb <= lofst_fe) and (lofst_fe <= lofst_le)) {
            DEBUG_ABORT_IF(fnode->get_data().get_tnode());
            lnode->get_data().set_hlink(fnode, lofst_fe - lofst_lb);
            fnode->get_data().set_tnode(lnode);
        }
    }

    inline void reset_lnode_links(lbuffer_type& lnodes) {
        for (lnode_type* lnode : lnodes) {
            reset_lnode_link(lnode);
        }
    }

    inline void reset_lnode_links_in_order(const lbuffer_type& lnodes) {
        if (lnodes.is_empty()) {
            return;
        }

#ifndef NDEBUG
        for (size_type i = 1; i < lnodes.get_size(); i++) {
            DEBUG_ABORT_IF(lnodes[i - 1]->get_next() != lnodes[i]);
        }
#endif

        lnode_type* lnode = lnodes.front();

        // Consider 'lnode' as the pivot, the L-interval is [lofst_lb, lofst_le].
        offset_type sum_lexps = get_sum_lexps(lnode);
        offset_type lofst_lb = 0;
        offset_type lofst_le = sum_lexps - 1;

        // Consider 'lnode' as the pivot, the F-interval is [lofst_fb, lofst_fe].
        auto [fnode, lofst_fe] = get_first_overlapped_fnode(lnode);
        offset_type sum_fexps = get_sum_fexps(fnode);
        offset_type lofst_fb = lofst_fe - sum_fexps + 1;

        const lnode_type* lnode_end = lnodes.back()->get_next();

        while (true) {
            bool lskip = false;
            bool fskip = false;

            if ((lofst_fb <= lofst_lb) and (lofst_lb <= lofst_fe) and (lofst_fe <= lofst_le)) {
                lnode->get_data().set_hlink(fnode, lofst_fe - lofst_lb);
                fnode->get_data().set_tnode(lnode);
                lskip = true;
                fskip = true;
            } else if (lofst_fe < lofst_le) {
                fskip = true;
            } else if (lofst_fe > lofst_le) {
                lskip = true;
            } else {
                ABORT_IF(true);
            }

            if (lskip) {
                lofst_fb = lofst_fb - sum_lexps;
                lofst_fe = lofst_fe - sum_lexps;

                lnode = lnode->get_next();
                if (lnode == lnode_end) {
                    break;
                }

                sum_lexps = get_sum_lexps(lnode);
                lofst_lb = 0;
                lofst_le = sum_lexps - 1;
            }

            if (fskip) {
                fnode = fnode->get_next();

                sum_fexps = get_sum_fexps(fnode);
                lofst_fe = lofst_fe + sum_fexps;
                lofst_fb = lofst_fe - sum_fexps + 1;
            }
        }
    }

    inline void reset_fnode_link(fnode_type* fnode) {
        if (fnode->get_data().get_tnode()) {
            return;
        }

        // Consider 'fnode' as the pivot, the F-interval is [fofst_fb, fofst_fe].
        const offset_type fofst_fb = 0;
        const offset_type fofst_fe = get_sum_fexps(fnode) - 1;

        // Consider 'fnode' as the pivot, the L-interval is [fofst_lb, fofst_le].
        auto [lnode, fofst_le] = get_first_overlapped_lnode(fnode);
        offset_type fofst_lb = fofst_le - get_sum_lexps(lnode) + 1;

        while (fofst_lb <= fofst_fe) {
            if ((fofst_fb <= fofst_lb) and (fofst_lb <= fofst_fe) and (fofst_fe <= fofst_le)) {
                DEBUG_ABORT_IF(lnode->get_data().get_hnode());
                lnode->get_data().set_hlink(fnode, fofst_fe - fofst_lb);
                fnode->get_data().set_tnode(lnode);
                break;
            }
            lnode = lnode->get_next();
            fofst_lb = fofst_le + 1;
            fofst_le = fofst_le + get_sum_lexps(lnode);
        }
    }

    inline void reset_fnode_links(fbuffer_type& fnodes) {
        for (fnode_type* fnode : fnodes) {
            reset_fnode_link(fnode);
        }
    }

    inline void reset_fnode_links_in_order(const fbuffer_type& fnodes) {
        if (fnodes.is_empty()) {
            return;
        }

#ifndef NDEBUG
        for (size_type i = 1; i < fnodes.get_size(); i++) {
            DEBUG_ABORT_IF(fnodes[i - 1]->get_next() != fnodes[i]);
        }
#endif

        fnode_type* fnode = fnodes.front();

        if (is_head_fnode(fnode)) {
            ABORT_IF_EQ(fnodes.get_size(), 1);
            fnode = fnode->get_next();
        }

        // Consider 'fnode' as the pivot, the F-interval is [fofst_fb, fofst_fe].
        offset_type sum_fexps = get_sum_fexps(fnode);
        offset_type fofst_fb = 0;
        offset_type fofst_fe = sum_fexps - 1;

        // Consider 'lnode' as the pivot, the F-interval is [lofst_fb, lofst_fe].
        auto [lnode, fofst_le] = get_first_overlapped_lnode(fnode);
        offset_type sum_lexps = get_sum_lexps(lnode);
        offset_type fofst_lb = fofst_le - sum_lexps + 1;

        const fnode_type* fnode_end = fnodes.back()->get_next();

        while (true) {
            bool fskip = false;
            bool lskip = false;

            if ((fofst_fb <= fofst_lb) and (fofst_lb <= fofst_fe) and (fofst_fe <= fofst_le)) {
                lnode->get_data().set_hlink(fnode, fofst_fe - fofst_lb);
                fnode->get_data().set_tnode(lnode);
                lskip = true;
                fskip = true;
            } else if (fofst_fe < fofst_le) {
                fskip = true;
            } else if (fofst_fe > fofst_le) {
                lskip = true;
            } else {
                lskip = true;
                fskip = true;
            }

            if (fskip) {
                fofst_lb = fofst_lb - sum_fexps;
                fofst_le = fofst_le - sum_fexps;

                fnode = fnode->get_next();
                if (fnode == fnode_end) {
                    break;
                }

                sum_fexps = get_sum_fexps(fnode);
                fofst_fb = 0;
                fofst_fe = sum_fexps - 1;
            }

            if (lskip) {
                lnode = lnode->get_next();

                sum_lexps = get_sum_lexps(lnode);
                fofst_le = fofst_le + sum_lexps;
                fofst_lb = fofst_le - sum_lexps + 1;
            }
        }
    }

#ifdef ENABLE_DEBUG_PRINT
    void debug_print() {
        tfm::printfln("");
        tfm::printfln("**** GroupedLFIntervalGraph::debug_print() ****");
        tfm::printfln("");

        tfm::printfln("=== $ ===");
        tfm::printfln("m_em_lnode=%s, m_em_lrnid=%d, m_em_lofst=%d", get_pc(m_em_lnode), m_em_lrnid, m_em_lofst);
        tfm::printfln("m_em_fnode=%s, m_em_frnid=%d, m_em_fofst=%d", get_pc(m_em_fnode), m_em_frnid, m_em_fofst);

        {
            tfm::printfln("\n=== F ===");
            const fnode_type* fnode = get_head_fnode();
            while (true) {
                tfm::printfln("%s", get_fnode_str(fnode));
                if (fnode == get_tail_fnode()) {
                    break;
                }
                fnode = fnode->get_next();
            }
        }

        {
            tfm::printfln("\n=== L ===");
            const lnode_type* lnode = get_head_lnode();
            while (true) {
                tfm::printfln("%s", get_lnode_str(lnode));
                if (is_dmmy_lnode(lnode)) {
                    break;
                }
                lnode = lnode->get_next();
            }
        }
        tfm::printfln("");
    }

    void update_pcmap(const void* p) {
        if (m_pcmap.find(intptr_t(p)) == m_pcmap.end()) {
            m_pcmap.insert(std::make_pair(intptr_t(p), m_pcmax++));
        }
    }
    std::string get_pc(const void* p) {
        update_pcmap(p);
        return tfm::format("[%c]", m_pcmap.find(intptr_t(p))->second);
    }
    std::string get_pc(const void* p) const {
        return tfm::format("[%c]", m_pcmap.find(intptr_t(p))->second);
    }

    std::string get_lnode_str(const lnode_type* p) {
        if (!p) return tfm::format("(%s)", get_pc(p));
        const auto& data = p->get_data();
        std::string str = tfm::format("%s => H=%s, O=%d, W=%d",  //
                                      get_pc(p), get_pc(data.get_hnode()), data.get_hofst(), data.get_weight());
        for (size_type i = 0; i < data.get_num_runs(); i++) {
            const lcursor_type crsr = data.get_lcursor_from_position(i);
            str += tfm::format("\n\t%s", get_lrun_str(p, crsr));
        }
        return str;
    }
    std::string get_fnode_str(const fnode_type* p) {
        if (!p) return tfm::format("(%s)", get_pc(p));
        const auto& data = p->get_data();
        std::string str = tfm::format("%s => T=%s, W=%d", get_pc(p), get_pc(data.get_tnode()), data.get_weight());
        for (size_type i = 0; i < data.get_num_runs(); i++) {
            const fcursor_type crsr = data.get_fcursor_from_position(i);
            str += tfm::format("\n\t%s", get_frun_str(p, crsr));
        }
        return str;
    }

    std::string get_litr_str(const literator_type itr) const {
        return tfm::format("run=%s, ofst=%d", get_lrun_str(itr.lnode, itr.lcrsr), itr.lofst);
    }
    std::string get_fitr_str(const fiterator_type itr) const {
        return tfm::format("run=%s, ofst=%d", get_frun_str(itr.fnode, itr.fcrsr), itr.fofst);
    }

    std::string get_lrun_str(const lnode_type* p, const lcursor_type crsr) const {
        if (!p) return "()";
        const auto& data = p->get_data();
        return tfm::format("(R=%s, F=(%s,%d), L=%d, SA=%d)",  //
                           make_run(data.get_chr(crsr), get_lexp(p, crsr)),  //
                           get_pc(data.get_frptr(crsr).first), data.get_frptr(crsr).second,  //
                           crsr.get_lrnid(),  //
                           data.get_sae(crsr));
    }
    std::string get_lrun_str(const lnode_type* p, const lcursor_type crsr) {
        if (!p) return "()";
        const auto& data = p->get_data();
        return tfm::format("(R=%s, F=(%s,%d), L=%d, SA=%d)",  //
                           make_run(data.get_chr(crsr), get_lexp(p, crsr)),  //
                           get_pc(data.get_frptr(crsr).first), data.get_frptr(crsr).second,  //
                           crsr.get_lrnid(),  //
                           data.get_sae(crsr));
    }
    std::string get_frun_str(const fnode_type* p, const fcursor_type crsr) const {
        if (!p) return "()";
        const auto& data = p->get_data();
        return tfm::format("(R=%s, L=(%s,%d), F=%d)", make_run(data.get_chr(crsr), get_fexp(p, crsr)),
                           get_pc(data.get_lrptr(crsr).first), data.get_lrptr(crsr).second, crsr.get_frnid());
    }
    std::string get_frun_str(const fnode_type* p, const fcursor_type crsr) {
        if (!p) return "()";
        const auto& data = p->get_data();
        return tfm::format("(R=%s, L=(%s,%d), F=%d)", make_run(data.get_chr(crsr), get_fexp(p, crsr)),
                           get_pc(data.get_lrptr(crsr).first), data.get_lrptr(crsr).second, crsr.get_frnid());
    }

    std::string get_lcrsr_str(const lcursor_type crsr) const {
        return tfm::format("(%d,%d)", crsr.get_lrnid(), crsr.get_lrpos());
    }
    std::string get_lcrsr_str(const lcursor_type crsr) {
        return tfm::format("(%d,%d)", crsr.get_lrnid(), crsr.get_lrpos());
    }
    std::string get_fcrsr_str(const fcursor_type crsr) const {
        return tfm::format("(%d,%d)", crsr.get_frnid(), crsr.get_frpos());
    }
    std::string get_fcrsr_str(const fcursor_type crsr) {
        return tfm::format("(%d,%d)", crsr.get_frnid(), crsr.get_frpos());
    }

    std::string get_lcover_str(const lcover_type& c) {
        return tfm::format("L=%s, F=%s, O=%d, N=%d", get_pc(c.lnode), get_pc(c.fnode), c.lofst, c.lwght);
    }
    std::string get_fcover_str(const fcover_type& c) {
        return tfm::format("F=%s, L=%s, O=%d, N=%d", get_pc(c.fnode), get_pc(c.lnode), c.fofst, c.fwght);
    }
#endif
};

}  // namespace rcomp
