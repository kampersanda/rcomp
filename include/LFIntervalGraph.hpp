/**
 * @file LFIntervalGraph.hpp
 */
#pragma once

#include <map>

#include "basics.hpp"
#include "defs.hpp"
#include "utils.hpp"

#include "FIndex.hpp"
#include "LIndex.hpp"

#include "StaticQueue.hpp"
#include "StaticVector.hpp"

namespace rcomp {

/**
 * An implementation of the straightforward LF-interval graph.
 *
 * @tparam t_LData The data type of L-node.
 * @tparam t_FData The data type of F-node.
 * @tparam t_DivBound The threshould for dividing fat nodes.
 *
 * @note The technical terms used are:
 *   - L/F-node means the L/F-run.
 *   - L/F-ofst means the offset of the character in the run.
 *   - L/F-addr means the 2-tuple of (node, ofst).
 *   - L/F-wght means the number of overlapped nodes (for the head characters).
 */
template <class t_LData, class t_FData, size_type t_DivBound = 7>
class LFIntervalGraph {
  public:
    using this_type = LFIntervalGraph<t_LData, t_FData, t_DivBound>;
    using lindex_type = LIndex<t_LData>;
    using findex_type = FIndex<t_FData>;
    using lnode_type = typename lindex_type::lnode_type;
    using fnode_type = typename findex_type::fnode_type;
    using ldata_type = t_LData;
    using fdata_type = t_FData;
    using lqueue_type = StaticQueue<lnode_type*, 4>;
    using fqueue_type = StaticQueue<fnode_type*, 4>;
    using lbuffer_type = StaticVector<lnode_type*, t_DivBound + 2>;
    using fbuffer_type = StaticVector<fnode_type*, t_DivBound + 2>;

    static constexpr size_type DIV_BOUND = t_DivBound;
    static_assert(7 <= DIV_BOUND, "The division bound must be no less than 7");

    enum class insert_modes : uint8_t {
        // Case C
        MERGE_PREV,
        MERGE_CURR,
        // Case B
        SPLIT_FRONT,
        SPLIT_MIDDLE,
        SPLIT_BACK,
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

  private:
    // The input original text size (with $).
    size_type m_num_chars = 0;

    // L/F data structures.
    lindex_type m_lindex;
    findex_type m_findex;

    // The L-address of $'s character.
    lnode_type* m_em_lnode = nullptr;
    offset_type m_em_lofst = -1;

    // The F-address overlapping the L-address.
    fnode_type* m_em_fnode = nullptr;
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
    LFIntervalGraph() = default;

    //! Default destructor
    virtual ~LFIntervalGraph() = default;

    //! Copy constructor (deleted)
    LFIntervalGraph(const LFIntervalGraph&) = delete;

    //! Copy constructor (deleted)
    LFIntervalGraph& operator=(const LFIntervalGraph&) = delete;

    //! Move constructor
    LFIntervalGraph(LFIntervalGraph&&) noexcept = default;

    //! Move constructor
    LFIntervalGraph& operator=(LFIntervalGraph&&) noexcept = default;

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
    inline size_type get_memory_in_bytes(bool include_this = true) const {
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
     *             Testers
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
        const fnode_type* fnode = m_findex.get_head();

        // Head F-node is not linked to dummy L-node.
        ABORT_IF(!is_dmmy_lnode(fnode->get_data().get_lnode()));

        fnode = fnode->get_next();
        while (true) {
            const fnode_type* fnext = fnode->get_next();
            if (is_head_fnode(fnext)) {
                break;
            }

            // F-node characters are not sorted in lex.
            ABORT_IF_LT(fnext->get_data().get_chr(), fnode->get_data().get_chr());

            if (fnode->get_data().get_chr() < fnext->get_data().get_chr()) {
                fnode = fnext;
                continue;
            }

            // L-node orders are not sorted for the F-interval.
            ABORT_IF_LE(lo_common::get_basic_order(fnext), lo_common::get_basic_order(fnode));
            fnode = fnext;
        }
    }

    //! Test the Head/Tail links of L/F-nodes.
    void test_ht_links() const {
        const fnode_type* fnode = get_head_fnode();
        const lnode_type* lnode = get_head_lnode();

        offset_type fpos_b = 0;  // indicate the front of fnode
        offset_type lpos_b = 0;  // indicate the front of lnode

        while (!is_dmmy_lnode(lnode)) {
            const offset_type fpos_e = fpos_b + get_fexp(fnode) - 1;  // indicate the back of fnode
            const offset_type lpos_e = lpos_b + get_lexp(lnode) - 1;  // indicate the back of lnode
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
            fpos_b = fpos_b + get_fexp(fnode);
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

        auto [lnode, lofst] = get_head_laddress();

        while (!(lnode == m_em_lnode and lofst == m_em_lofst)) {
            fn(lnode->get_data().get_chr());

            auto fnode = lnode->get_data().get_fnode();
            auto fofst = lofst;

            if (lnode == m_em_lnode and m_em_lofst < lofst) {
                fofst -= 1;  // due to $
            }

            std::tie(lnode, lofst) = get_overlapped_laddress(fnode, fofst);
        }
    }

    /**
     * @brief Output the RLBWT text.
     * @param[in] fn The callback function to get the RLBWT text run by run.
     */
    void output_runs(const std::function<void(const run_type&)>& fn) const {
        if (is_empty()) {
            return;
        }

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

        for (auto lnode = get_head_lnode();; lnode = lnode->get_next()) {
            if (lnode != m_em_lnode) {
                update_run(make_run(lnode->get_data().get_chr(), get_lexp(lnode)));
            } else {
                // $-node
                const offset_type this_exp = get_lexp(lnode) - 1;  // without $

                if (m_em_lofst == 0) {
                    // Push front
                    update_run(make_run(END_MARKER, 1));
                    update_run(make_run(lnode->get_data().get_chr(), this_exp));
                } else if (m_em_lofst == this_exp) {
                    // Push back
                    update_run(make_run(lnode->get_data().get_chr(), this_exp));
                    update_run(make_run(END_MARKER, 1));
                } else {
                    const offset_type new_exp1 = m_em_lofst;
                    const offset_type new_exp2 = this_exp - m_em_lofst;
                    update_run(make_run(lnode->get_data().get_chr(), new_exp1));
                    update_run(make_run(END_MARKER, 1));
                    update_run(make_run(lnode->get_data().get_chr(), new_exp2));
                }
            }
            if (lnode == get_tail_lnode()) {
                break;
            }
        }

        fn(prev_rn);
    }

    /* * * * * * * * * * * * * * * * * * *
     *
     *            Initilaizer
     *
     * * * * * * * * * * * * * * * * * * */

    //! Make the data structure for the RLBWT text of T = (new_chr,$).
    auto extend_init(uchar_type new_chr) -> std::tuple<lnode_type*, lnode_type*> {
        // Dummy node for the partner of $'s F-node (for convenience)
        m_lindex.clear();
        lnode_type* ldmmy = m_lindex.get_head();

        // $'s node at the top of F
        m_findex.clear(ldmmy);
        fnode_type* fhead = m_findex.get_head();

        lnode_type* lnode = m_lindex.insert_after(ldmmy, new_chr, 1);
        fnode_type* fnode = m_findex.insert_after(fhead, lnode);

        // 'lnode' is always located at the head of L-list.
        // So, the tail link indicates 'fhead' ($-node) and will be never changed.
        fhead->get_data().set_tnode(lnode);
        lnode->get_data().set_hnode(fhead);
        lnode->get_data().set_hofst(0);

        fnode->get_data().set_tnode(nullptr);
        ldmmy->get_data().set_hnode(nullptr);
        ldmmy->get_data().set_hofst(0);

        // The last character on L always becomes $
        m_em_lnode = lnode;
        m_em_lofst = 1;

        m_em_fnode = fnode;
        m_em_fofst = 0;

        fhead->get_data().set_weight(1);  // for lnode
        ldmmy->get_data().set_weight(0);
        fnode->get_data().set_weight(0);
        lnode->get_data().set_weight(2);  // for fhead and fnode

        m_num_chars = 2;

        return std::make_tuple(lnode, ldmmy);
    }

    /* * * * * * * * * * * * * * * * * * *
     *
     *        Phase A in Case C/B
     *
     * * * * * * * * * * * * * * * * * * */

    //! Check if the insertion mode is in merge.
    inline bool is_mergable_mode(insert_modes ins_mode) const {
        return ins_mode == insert_modes::MERGE_PREV || ins_mode == insert_modes::MERGE_CURR;
    }

    //! Check if $'s L-node can be merged or split.
    inline insert_modes check_mergable(const uchar_type new_chr) const {
        DEBUG_ABORT_IF(is_dmmy_lnode(m_em_lnode));
        DEBUG_ABORT_IF_OUT(m_em_lofst, 0, get_lexp(m_em_lnode) - 1);  // -1 is for $

        // Mergeable with the current node?
        if (new_chr == m_em_lnode->get_data().get_chr()) {
            return insert_modes::MERGE_CURR;  // Case C
        }

        if (m_em_lofst == 0) {
            // Mergeable with the previous node?
            if (new_chr == m_em_lnode->get_prev()->get_data().get_chr()) {
                return insert_modes::MERGE_PREV;  // Case C
            }
        }

        // new_chr will not be merged with the next run because $-mark is never placed at the end of some L-run.
        // Only when $-mark is placed at the last of L, it will be placed at the end of the last L-run.
        // But, then new_chr is not be merged with the next run anyway.

        if (m_em_lofst == 0) {
            return insert_modes::SPLIT_FRONT;
        } else if (m_em_lofst == get_lexp(m_em_lnode) - 1) {  // -1 is for $
            ABORT_IF(m_em_lnode != get_tail_lnode());
            return insert_modes::SPLIT_BACK;
        } else {
            return insert_modes::SPLIT_MIDDLE;
        }
    }

    /* * * * * * * * * * * * * * * * * * *
     *
     *       Phase B -- F in Case C
     *
     * * * * * * * * * * * * * * * * * * */

    //! Extend the RLBWT text in Case C.
    auto extend_with_merge(const uchar_type, const insert_modes ins_mode)
        -> std::tuple<lnode_type*, fnode_type*, offset_type> {
        DEBUG_ABORT_IF(!is_mergable_mode(ins_mode));
        DEBUG_PRINT(tfm::printfln("\n==== Case C ====");)

        DEBUG_ABORT_IF(!m_fat_lnodes.is_empty());
        DEBUG_ABORT_IF(!m_fat_fnodes.is_empty());

        m_num_chars += 1;

        lnode_type* ins_lnode = nullptr;
        fnode_type* ins_fnode = nullptr;
        offset_type ins_fofst = 0;

        // NOTE:
        // Depending on 'ins_mode', the L-node with $ can be updated, but Do NOT change the value of 'm_em_lnode' until
        // Phase E. This is because this change will affect 'get_lexp'.

        //
        // Phase B: Set the L/F-nodes to be updated
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase B ==");)

        if (ins_mode == insert_modes::MERGE_CURR) {
            // Then, the new char will be inserted into 'm_em_lnode'
            DEBUG_ABORT_IF_OUT(m_em_lofst, 0, get_lexp(m_em_lnode) - 1);
            ins_lnode = m_em_lnode;
            ins_fnode = ins_lnode->get_data().get_fnode();
            ins_fofst = m_em_lofst;
        } else if (ins_mode == insert_modes::MERGE_PREV) {
            // Then, the new char will be inserted into m_em_lnode->get_prev()
            DEBUG_ABORT_IF_NE(m_em_lofst, 0);
            ins_lnode = m_em_lnode->get_prev();
            ins_fnode = ins_lnode->get_data().get_fnode();
            ins_fofst = get_fexp(ins_fnode);
        } else {
            ABORT_IF(true);
        }

        //
        // Phase C: Update H/T-links and -weights for merging $'s L-node
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase C ==");)

        // $-replacement
        if (ins_mode == insert_modes::MERGE_PREV) {
            lnode_type* em_tnode = m_em_fnode->get_data().get_tnode();

            const bool is_first_em_fofst = is_first_faddress(m_em_fnode, m_em_fofst);
            const bool is_last_em_fofst = is_last_faddress(m_em_fnode, m_em_fofst);

            // Weights
            if (is_first_em_fofst) {
                lnode_type* prev_lnode = m_em_lnode->get_prev();
                m_em_lnode->get_data().sub_weight(1);
                prev_lnode->get_data().add_weight(1);
                push_fat_lnode(prev_lnode);
            }
            if (is_last_em_fofst) {
                fnode_type* next_fnode = m_em_fnode->get_next();
                m_em_fnode->get_data().sub_weight(1);
                next_fnode->get_data().add_weight(1);
                push_fat_fnode(next_fnode);
            }

            // Linkage
            if (is_last_em_fofst) {
                DEBUG_ABORT_IF(em_tnode != m_em_lnode);
                DEBUG_ABORT_IF_NE(em_tnode->get_data().get_hofst(), 0);
                {
                    fnode_type* fnode_befo = m_em_fnode;
                    lnode_type* lnode_befo = m_em_lnode->get_prev();
                    const offset_type new_fexp = fnode_befo->get_data().get_exp();
                    const offset_type new_lexp = lnode_befo->get_data().get_exp() + 1;
                    if (new_lexp <= new_fexp) {
                        DEBUG_ABORT_IF(lnode_befo->get_data().get_hnode());
                        lnode_befo->get_data().set_hlink(fnode_befo, new_lexp - 1);
                        fnode_befo->get_data().set_tnode(lnode_befo);
                    } else {
                        DEBUG_ABORT_IF(!lnode_befo->get_data().get_hnode());
                        fnode_befo->get_data().reset_tlink();
                    }
                }
                {
                    fnode_type* fnode_aftr = m_em_fnode->get_next();
                    lnode_type* lnode_aftr = m_em_lnode;
                    const offset_type new_fexp = fnode_aftr->get_data().get_exp();
                    const offset_type new_lexp = lnode_aftr->get_data().get_exp();
                    if (new_lexp >= new_fexp) {
                        DEBUG_ABORT_IF(fnode_aftr->get_data().get_tnode());
                        lnode_aftr->get_data().set_hlink(fnode_aftr, new_fexp - 1);
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
        // Phase D: Merge $'s L-node
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase D ==");)

        ins_lnode->get_data().add_exp(1);

        DEBUG_PRINT(tfm::printfln("ins_lnode=%s", get_lnode_str(ins_lnode));)

        //
        // Phase E: Compute the new $'s position
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase E ==");)

        // Clear since the L-run has already updated.
        m_em_lnode = nullptr;
        m_em_lofst = 0;

        // Cache to the next process
        m_em_fnode = ins_fnode;
        m_em_fofst = ins_fofst;

        std::tie(m_em_lnode, m_em_lofst) = get_overlapped_laddress(ins_fnode, ins_fofst);

        DEBUG_PRINT(tfm::printfln("m_em_lnode=%s", get_lnode_str(m_em_lnode));)
        DEBUG_PRINT(tfm::printfln("m_em_lofst=%d", m_em_lofst);)

        //
        // Phase F: Update H/T-links and -weights for putting $-mark (?-replacement)
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase F ==");)

        {
            const bool is_last_em_fofst = is_last_faddress(m_em_fnode, m_em_fofst);
            const bool is_first_em_lofst = is_first_laddress(m_em_lnode, m_em_lofst);

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
        }

        return std::make_tuple(ins_lnode, ins_fnode, ins_fofst);
    }

    /* * * * * * * * * * * * * * * * * * *
     *
     *       Phase B -- F in Case B
     *
     * * * * * * * * * * * * * * * * * * */

    //! Divide $'s L-node in advance.
    insert_modes divide_if_need(insert_modes ins_mode) {
        DEBUG_ABORT_IF(!m_fat_lnodes.is_empty());
        DEBUG_ABORT_IF(!m_fat_fnodes.is_empty());

        if (ins_mode == insert_modes::SPLIT_MIDDLE) {
            divide_lnode(m_em_lnode, m_em_lofst);
            ins_mode = insert_modes::SPLIT_FRONT;
        }
        return ins_mode;
    }

    //! Extend the BWT text in Case B.
    void extend_with_split(const uchar_type new_chr, insert_modes ins_mode) {
        DEBUG_ABORT_IF(is_mergable_mode(ins_mode));
        DEBUG_ABORT_IF(ins_mode == insert_modes::SPLIT_MIDDLE);
        DEBUG_PRINT(tfm::printfln("\n==== Case B ====");)

        m_num_chars += 1;

        //
        // Phase B: Search the L/F-nodes to be updated
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase B ==");)

        lnode_type* lpred = ins_mode == insert_modes::SPLIT_FRONT ? m_em_lnode->get_prev() : m_em_lnode;
        fnode_type* fpred = m_findex.predecessor(new_chr, lpred->get_data().get_fnode());

        DEBUG_PRINT(tfm::printfln("lpred=%s, fpred=%s", get_pc(lpred), get_pc(fpred));)

        //
        // Phase C: Remove H/T-links of overlapped nodes
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase C =="););

        lnode_type* ins_lnode = m_em_lnode;
        fnode_type* ins_fnode = get_tail_fnode() != fpred ? fpred->get_next() : fpred;
        get_overlapped_fnodes(ins_lnode, 0, get_lexp(ins_lnode), m_fnodes_buffer);
        get_overlapped_lnodes(ins_fnode, 0, get_fexp(ins_fnode), m_lnodes_buffer);
        remove_fnode_links(m_fnodes_buffer);
        remove_lnode_links(m_lnodes_buffer);

        DEBUG_PRINT(tfm::printf("m_lnodes_buffer:");)
        DEBUG_PRINT(for (const auto* p : m_lnodes_buffer) tfm::printf(" %s", get_pc(p));)
        DEBUG_PRINT(tfm::printf("\n");)
        DEBUG_PRINT(tfm::printf("m_fnodes_buffer:");)
        DEBUG_PRINT(for (const auto* p : m_fnodes_buffer) tfm::printf(" %s", get_pc(p));)
        DEBUG_PRINT(tfm::printf("\n");)

        //
        // Phase D: Insert the new character
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase D ==");)

        lnode_type* new_lnode = m_lindex.insert_after(lpred, new_chr, 1);
        fnode_type* new_fnode = m_findex.insert_after(fpred, new_lnode);

        DEBUG_PRINT(tfm::printfln("new_lnode =%s", get_lnode_str(new_lnode));)
        DEBUG_PRINT(tfm::printfln("new_fnode =%s", get_fnode_str(new_fnode));)

        //
        // Phase E: Compute new $'s position
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase E ==");)

        // Clear since the L-run has already updated.
        m_em_lnode = nullptr;
        m_em_lofst = 0;

        // Cache to the next process
        m_em_fnode = new_fnode;
        m_em_fofst = 0;

        std::tie(m_em_lnode, m_em_lofst) = get_overlapped_laddress(new_fnode, 0);

        DEBUG_PRINT(tfm::printfln("m_em_lnode=%s", get_lnode_str(m_em_lnode));)
        DEBUG_PRINT(tfm::printfln("m_em_lofst=%d", m_em_lofst);)

        //
        // Phase F: Reset H/T-links of overlapped nodes
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase F ==");)

        reset_lnode_links_and_weights(m_lnodes_buffer);
        reset_fnode_links_and_weights(m_fnodes_buffer);
        reset_lnode_link_and_weight(ins_lnode);
        reset_lnode_link_and_weight(new_lnode);
        reset_fnode_link_and_weight(ins_fnode);
        reset_fnode_link_and_weight(new_fnode);

        // Release
        m_fnodes_buffer.clear();
        m_lnodes_buffer.clear();
    }

    /* * * * * * * * * * * * * * * * * * *
     *
     *      Phase B -- F in Case L/F
     *
     * * * * * * * * * * * * * * * * * * */

    //! Divide a fat node in Case L/F.
    bool divide_fat_node() {
        // Case L
        while (!m_fat_lnodes.is_empty()) {
            lnode_type* lnode = m_fat_lnodes.pop();
            DEBUG_ABORT_IF(!lnode);
            if (lnode->get_data().get_weight() >= DIV_BOUND) {
                const lcover_type lcover = get_lcover(lnode);
                DEBUG_ABORT_IF_NE(lcover.lwght, lnode->get_data().get_weight());
                divide_lcover(lcover);
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
                divide_fcover(fcover);
                return true;
            }
        }
        return false;
    }

    //! Divide a fat L-node in Case L.
    auto divide_fat_lnode() -> std::tuple<lnode_type*, lnode_type*> {
        while (!m_fat_lnodes.is_empty()) {
            lnode_type* lnode = m_fat_lnodes.pop();
            DEBUG_ABORT_IF(!lnode);
            if (lnode->get_data().get_weight() >= DIV_BOUND) {
                const lcover_type lcover = get_lcover(lnode);
                DEBUG_ABORT_IF_NE(lcover.lwght, lnode->get_data().get_weight());
                return divide_lcover(lcover);
            }
        }
        return {nullptr, nullptr};
    }

    //! Divide a fat F-node in Case F.
    auto divide_fat_fnode() -> std::tuple<lnode_type*, lnode_type*> {
        while (!m_fat_fnodes.is_empty()) {
            fnode_type* fnode = m_fat_fnodes.pop();
            DEBUG_ABORT_IF(!fnode);
            if (fnode->get_data().get_weight() >= DIV_BOUND) {
                const fcover_type fcover = get_fcover(fnode);
                DEBUG_ABORT_IF_NE(fcover.fwght, fnode->get_data().get_weight());
                return divide_fcover(fcover);
            }
        }
        return {nullptr, nullptr};
    }

    //! Return the first part of the divided L-node, u, and the L-node corresponding to LF^{-1}(u).
    auto divide_lcover(const lcover_type& lcover) -> std::tuple<lnode_type*, lnode_type*> {
        DEBUG_PRINT(tfm::printfln("\n==== Case L ====");)
        DEBUG_PRINT(tfm::printfln("lcover => %d", get_lcover_str(lcover));)

        // Phase B
        lnode_type* div_lnode = lcover.lnode;  // to be split
        offset_type div_lofst = lcover.lofst;

        const fnode_type* lfm_fnode = lcover.fnode;
        {
            const size_type fnum1 = lcover.lwght / 2;
            for (size_type i = 1; i < fnum1; i++) {
                lfm_fnode = lfm_fnode->get_next();
                div_lofst += get_fexp(lfm_fnode);
            }
        }
        // Here, div_lnode[div_lofst] is the last character of the front part

        // Phase C--F
        if (div_lnode == m_em_lnode and m_em_lofst <= div_lofst) {
            divide_lnode(div_lnode, div_lofst);
        } else {
            divide_lnode(div_lnode, div_lofst + 1);
        }
        return std::make_tuple(div_lnode, lfm_fnode->get_data().get_lnode());
    }

    //! Return the first part of the divided L-node, u, and the L-node corresponding to LF(u).
    auto divide_fcover(const fcover_type& fcover) -> std::tuple<lnode_type*, lnode_type*> {
        DEBUG_PRINT(tfm::printfln("\n==== Case F ====");)
        DEBUG_PRINT(tfm::printfln("fcover => %d", get_fcover_str(fcover));)

        // Phase B
        fnode_type* div_fnode = fcover.fnode;  // to be split
        lnode_type* div_lnode = div_fnode->get_data().get_lnode();
        offset_type div_fofst = fcover.fofst;

        lnode_type* lfm_lnode = fcover.lnode;  // the last partner with div_fnode
        {
            const size_type lnum1 = fcover.fwght / 2;
            for (size_type i = 1; i < lnum1; i++) {
                lfm_lnode = lfm_lnode->get_next();
                div_fofst += get_lexp(lfm_lnode);
            }
        }
        // Here, div_fnode[div_fofst] is the last character of the front part

        // Phase C--F
        divide_lnode(div_lnode, div_fofst + 1);
        return std::make_tuple(div_lnode, lfm_lnode);
    }

    //! Divide the L-node with $-mark and the corresponding F-node.
    //! As the result, $-mark will be placed at the first of the after L-node.
    lnode_type* divide_lnode(lnode_type* div_lnode, const offset_type new_exp1) {
        DEBUG_ABORT_IF(!m_lnodes_buffer.is_empty());
        DEBUG_ABORT_IF(!m_fnodes_buffer.is_empty());

        fnode_type* div_fnode = div_lnode->get_data().get_fnode();

        const offset_type div_lexp = get_lexp(div_lnode);
        const offset_type div_fexp = get_fexp(div_fnode);
        const offset_type new_exp2 = div_fexp - new_exp1;  // Note div_lexp can include that of $-mark.

        ABORT_IF_LE(new_exp1, 0);
        ABORT_IF_LE(new_exp2, 0);

        //
        //  Phase C
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase C =="););

        get_overlapped_fnodes(div_lnode, 0, div_lexp, m_fnodes_buffer);
        get_overlapped_lnodes(div_fnode, 0, div_fexp, m_lnodes_buffer);
        remove_fnode_links(m_fnodes_buffer);
        remove_lnode_links(m_lnodes_buffer);

        DEBUG_PRINT(tfm::printf("m_lnodes_buffer:");)
        DEBUG_PRINT(for (const auto* p : m_lnodes_buffer) tfm::printf(" %s", get_pc(p));)
        DEBUG_PRINT(tfm::printf("\n");)
        DEBUG_PRINT(tfm::printf("m_fnodes_buffer:");)
        DEBUG_PRINT(for (const auto* p : m_fnodes_buffer) tfm::printf(" %s", get_pc(p));)
        DEBUG_PRINT(tfm::printf("\n");)

        //
        // Phase D
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase D ==");)

        div_lnode->get_data().set_exp(new_exp1);

        lnode_type* new_lnode = m_lindex.insert_after(div_lnode, div_lnode->get_data().get_chr(), new_exp2);
        fnode_type* new_fnode = m_findex.insert_after(div_fnode, new_lnode);

        DEBUG_PRINT(tfm::printfln("new_lnode =%s", get_lnode_str(new_lnode));)
        DEBUG_PRINT(tfm::printfln("new_fnode =%s", get_fnode_str(new_fnode));)

        //
        // Phase E
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase E ==");)

        if (div_lnode == m_em_lnode and new_exp1 <= m_em_lofst) {
            m_em_lnode = new_lnode;
            m_em_lofst = m_em_lofst - new_exp1;
        }
        if (div_fnode == m_em_fnode and new_exp1 <= m_em_fofst) {
            m_em_fnode = new_fnode;
            m_em_fofst = m_em_fofst - new_exp1;
        }

        //
        // Phase F
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase F ==");)

        reset_fnode_links_and_weights(m_fnodes_buffer);
        reset_lnode_links_and_weights(m_lnodes_buffer);
        reset_lnode_link_and_weight(div_lnode);
        reset_lnode_link_and_weight(new_lnode);
        reset_fnode_link_and_weight(div_fnode);
        reset_fnode_link_and_weight(new_fnode);

        // Release
        m_fnodes_buffer.clear();
        m_lnodes_buffer.clear();

        return new_lnode;
    }

    /* * * * * * * * * * * * * * * * * * * * * * *
     *
     *       Basic handlers for LF-nodes
     *
     * * * * * * * * * * * * * * * * * * * * * * */

    //! Get the L-address indicating $-marker.
    inline std::tuple<lnode_type*, offset_type> get_em_laddress() const {
        return {m_em_lnode, m_em_lofst};
    }

    //! Get the F-address overlapping the L-address.
    inline std::tuple<fnode_type*, offset_type> get_em_faddress() const {
        return {m_em_fnode, m_em_fofst};
    }

    //! Get the L-address of the head character.
    inline std::tuple<lnode_type*, offset_type> get_head_laddress() const {
        return {get_head_lnode(), 0};
    }

    //! Get the L-address of the tail character.
    inline std::tuple<lnode_type*, offset_type> get_tail_laddress() const {
        return {get_tail_lnode(), get_lexp(get_tail_lnode()) - 1};
    }

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

    //! Get the exponent of the given L-run considering $-marker.
    inline offset_type get_lexp(const lnode_type* lnode) const {
        const offset_type lexp = lnode->get_data().get_exp();
        return lnode == m_em_lnode ? lexp + 1 : lexp;
    }

    //! Get the exponent of the given F-run.
    inline offset_type get_fexp(const fnode_type* fnode) const {
        return fnode->get_data().get_exp();
    }

    //! Check if the given offset indicates the first char of the given L-node.
    inline bool is_first_laddress(const lnode_type*, offset_type lofst) const {
        return lofst == 0;
    }

    //! Check if the given offset indicates the last char of the given L-node, considering $-marker.
    inline bool is_last_laddress(const lnode_type* lnode, offset_type lofst) const {
        return get_lexp(lnode) == lofst + 1;
    }

    //! Check if the given offset indicates the first char of the given F-node.
    inline bool is_first_faddress(const fnode_type*, offset_type fofst) const {
        return fofst == 0;
    }

    //! Check if the given offset indicates the last char of the given F-node.
    inline bool is_last_faddress(const fnode_type* fnode, offset_type fofst) const {
        return get_fexp(fnode) == fofst + 1;
    }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
     *
     *      Backward searchers (in log time)
     *
     * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    //! Predeccessor query in backward search.
    inline std::tuple<const fnode_type*, offset_type, bool>  //
    predecessor_in_bwsearch(uchar_type q_chr, const lnode_type* q_lnode, offset_type q_lofst) const {
        DEBUG_ABORT_IF_EQ(q_chr, END_MARKER);  // Currently, search for END_MARKER is not supported.
        DEBUG_ABORT_IF(m_em_lnode == nullptr);
        DEBUG_ABORT_IF(m_em_fnode == nullptr);

        const fnode_type* lfm_fnode = nullptr;
        offset_type lfm_fofst = -1;
        bool matched = false;

        if (q_chr == q_lnode->get_data().get_chr()) {
            // The LF-mapped F-node is the linked one.
            lfm_fnode = q_lnode->get_data().get_fnode();

            // SPECIAL CASE: The cursor indicates $-marker
            if (q_lnode == m_em_lnode and q_lofst == m_em_lofst) {
                if (is_first_laddress(m_em_lnode, m_em_lofst)) {
                    // Visit the previous L-node instead of the decrement
                    lfm_fnode = lfm_fnode->get_prev();
                    if (q_chr != lfm_fnode->get_data().get_chr()) {
                        return {nullptr, -1, false};
                    }
                    lfm_fofst = get_fexp(lfm_fnode) - 1;
                } else {
                    // decrement due to the virtual $-marker
                    lfm_fofst = q_lofst - 1;
                }
                matched = false;  // NOT MATCHED!!!
            }
            // SPECIAL CASE: The L-node has $-marker (but, not indicating)
            else if (q_lnode == m_em_lnode and q_lofst != m_em_lofst) {
                if (q_lofst < m_em_lofst) {
                    lfm_fofst = q_lofst;
                } else {
                    // decrement due to the virtual $-marker
                    lfm_fofst = q_lofst - 1;
                }
                matched = true;
            }
            // NORMAL CASE
            else {
                lfm_fofst = q_lofst;
                matched = true;
            }
        } else {
            lfm_fnode = m_findex.exact_predecessor(q_chr, q_lnode->get_data().get_fnode());
            if (lfm_fnode == nullptr) {
                return {nullptr, -1, false};
            }
            lfm_fofst = get_fexp(lfm_fnode) - 1;
            matched = false;
        }

        return {lfm_fnode, lfm_fofst, matched};
    }

    //! Successor query in backward search.
    inline std::tuple<const fnode_type*, offset_type, bool>  //
    successor_in_bwsearch(uchar_type q_chr, const lnode_type* q_lnode, offset_type q_lofst) const {
        // Currently, search for END_MARKER is not supported.
        DEBUG_ABORT_IF_EQ(q_chr, END_MARKER);
        DEBUG_ABORT_IF(m_em_lnode == nullptr);
        DEBUG_ABORT_IF(m_em_fnode == nullptr);

        const fnode_type* lfm_fnode = nullptr;
        offset_type lfm_fofst = -1;
        bool matched = false;

        if (q_chr == q_lnode->get_data().get_chr()) {
            // The LF-mapped F-node is the linked one.
            lfm_fnode = q_lnode->get_data().get_fnode();

            // SPECIAL CASE: The cursor indicates $-marker
            if (q_lnode == m_em_lnode and q_lofst == m_em_lofst) {
                if (is_last_laddress(m_em_lnode, m_em_lofst)) {
                    // This case will arise only when $-marker is placed at the tail of L.
                    // So, there is no the successor.
                    return {nullptr, -1, false};
                } else {
                    lfm_fofst = q_lofst;
                }
                matched = false;  // NOT MATCHED!!!
            }
            // SPECIAL CASE: The L-node has $-marker (but, not indicating)
            else if (q_lnode == m_em_lnode and q_lofst != m_em_lofst) {
                if (q_lofst < m_em_lofst) {
                    lfm_fofst = q_lofst;
                } else {
                    // decrement due to the virtual $-marker
                    lfm_fofst = q_lofst - 1;
                }
                matched = true;
            }
            // NORMAL CASE
            else {
                lfm_fofst = q_lofst;
                matched = true;
            }
        } else {
            lfm_fnode = m_findex.exact_successor(q_chr, q_lnode->get_data().get_fnode());
            if (lfm_fnode == nullptr) {
                return {nullptr, -1, false};
            }
            lfm_fofst = 0;
            matched = false;
        }

        return {lfm_fnode, lfm_fofst, matched};
    }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
     *
     *      Getters for linkages (in constant time)
     *
     * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    //! Get the distance between the begin positions of the given F-node and the tail-linked L-node.
    inline offset_type get_distance_for_begins(const fnode_type* fnode) const {
        DEBUG_ABORT_IF(!fnode->get_data().get_tnode());
        DEBUG_ABORT_IF_LT(0, fnode->get_data().get_tofst() - get_fexp(fnode) + 1);
        return -1 * (fnode->get_data().get_tofst() - get_fexp(fnode) + 1);
    }

    //! Return (lnode, lofst) such that q_fnode[q_fofst] overlaps lnode[lofst].
    inline std::tuple<lnode_type*, offset_type>  //
    get_overlapped_laddress(const fnode_type* q_fnode, offset_type q_fofst) const {
        DEBUG_ABORT_IF(!q_fnode);
        DEBUG_ABORT_IF(is_head_fnode(q_fnode));

        // Exception case: \Lambda_t indicates the tail character of F?
        if (q_fnode == get_tail_fnode() and get_fexp(q_fnode) == q_fofst + 1) {
            lnode_type* ltail = get_tail_lnode();
            if (m_em_lnode != nullptr) {
                return {ltail, get_lexp(ltail) - 1};
            } else {
                // It will be called during update
                // Then, lofst to be returned indicates exclusive position (only in this case).
                return {ltail, get_lexp(ltail)};
            }
        }

        auto [lnode, fofst] = get_first_overlapped_lnode(q_fnode);

        while (fofst < q_fofst) {
            lnode = lnode->get_next();
            fofst += get_lexp(lnode);
        }
        return {lnode, (get_lexp(lnode) - 1) - (fofst - q_fofst)};
    }

    //! Return (fnode, fofst) such that q_lnode[q_lofst] overlaps fnode[fofst].
    inline std::tuple<fnode_type*, offset_type>  //
    get_overlapped_faddress(const lnode_type* q_lnode, offset_type q_lofst) const {
        DEBUG_ABORT_IF(!q_lnode);
        DEBUG_ABORT_IF(is_dmmy_lnode(q_lnode));

        auto [fnode, lofst] = get_first_overlapped_fnode(q_lnode);

        while (lofst < q_lofst) {
            fnode = fnode->get_next();
            lofst += get_fexp(fnode);
        }
        return {fnode, (get_fexp(fnode) - 1) - (lofst - q_lofst)};
    }

    //! Return (lnode, fofst) such that fnode[fofst] overlaps the back of lnode.
    //!
    //!    e.g.,
    //!               [F]      [L]
    //!            v1  a
    //!                a        b  u1
    //!            v2  a        b
    //!                a        c  u2 *      The output is (u2, 1)
    //!       (in) v3  b        c            since v3[1] overlaps u2.back().
    //!                b  --->  c
    //!                b        b  u3
    //!                         b
    //!
    std::pair<lnode_type*, offset_type> get_first_overlapped_lnode(const fnode_type* fnode) const {
        if (is_head_fnode(fnode)) {
            lnode_type* lnode = get_head_lnode();
            return {lnode, get_lexp(lnode) - 1};
        }

        // Distance from the head of the given F-node
        offset_type fofst = 0;

        // Move to the closest F-node with a tail link
        // (This will be done since the head nodes are always linked)
        fnode = fnode->get_prev();  // start at the previous F-node
        while (!fnode->get_data().get_tnode()) {
            fofst -= get_fexp(fnode);
            fnode = fnode->get_prev();
        }
        fofst -= get_fexp(fnode);

        // Move to the tail-linked L-node
        lnode_type* lnode = fnode->get_data().get_tnode();
        fofst += get_distance_for_begins(fnode);

        // Move to the first overlapped L-node with the given F-node
        do {
            fofst += get_lexp(lnode);
            lnode = lnode->get_next();
        } while (fofst <= 0);

        return {lnode->get_prev(), fofst - 1};
    }

    //! Return (fnode, lofst) such that lnode[lofst] overlaps the back of fnode.
    //!
    //!     e.g.,
    //!               [F]      [L]
    //!            v1  a
    //!                a        b  u1
    //!          * v2  a        b
    //!                a        c  u2 (in)   The output is (v2, 2)
    //!                a        c            since u2[2] overlaps v2.back.
    //!                a  <---  c
    //!            v3  b        b  u3
    //!                b        b
    //!
    std::pair<fnode_type*, offset_type> get_first_overlapped_fnode(const lnode_type* lnode) const {
        DEBUG_ABORT_IF(is_dmmy_lnode(lnode));

        // Distance from the head of the given L-node
        offset_type lofst = 0;

        // Move to the closest L-node with a head link
        // (This will be always done since the head L-node has the head link to the head F-node)
        while (!lnode->get_data().get_hnode()) {
            lnode = lnode->get_prev();
            lofst -= get_lexp(lnode);
        }

        // Move to the tail-linked L-node
        fnode_type* fnode = lnode->get_data().get_hnode();
        lofst -= get_distance_for_begins(fnode);

        // Move to the first overlapped F-node with the given L-node
        while (true) {
            lofst += get_fexp(fnode);
            if (lofst > 0) {  // overlapped?
                break;
            }
            fnode = fnode->get_next();
        }

        return {fnode, lofst - 1};
    }

    //! Return the L-nodes overlapped with the F-interval of size 'fsize' starting at the 'fofst'-th of 'fnode'.
    //!
    //!     e.g.,
    //!                     [F]      [L]
    //!                  v1  a
    //!                      a        b  u1
    //!             (in) v2  a        b
    //!     fofst = 1 -->    a +    + c  u2 *       The output is [u2, u3].
    //!                  v3  b |    | c
    //!                      b |    | c
    //!     fsize = 4 -->    b +    | b  u3 *
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
            fofst += get_lexp(lnode);
        }

        // overlapped?
        while (!is_dmmy_lnode(lnode)) {
            buf.push_back(lnode);
            if (fofst_e <= fofst) {
                break;
            }
            lnode = lnode->get_next();
            fofst += get_lexp(lnode);
        }
    }

    //! Return the F-nodes overlapped with the L-interval of size 'lsize' starting at the 'lofst'-th of 'lnode'.
    //!
    //!     e.g.,
    //!                 [F]      [L]
    //!              v1  a                          The output is [v2, v3].
    //!                  a        b  u1 (in)
    //!            * v2  a +      b
    //!                  a |    + c  u2  <-- lofst = 2
    //!            * v3  b |    | c
    //!                  b |    + c      <-- lsize = 3
    //!                  b +      b  u3
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
            lofst += get_fexp(fnode);
        }

        // overlapped?
        while (true) {
            buf.push_back(fnode);
            if (lofst_e <= lofst) {
                break;
            }
            fnode = fnode->get_next();
            lofst += get_fexp(fnode);
            if (is_head_fnode(fnode)) {  // cycled?
                break;
            }
        }
    }

    //! Return lcover_type(lnode, fnode, lofst, lwght) such that
    //!  - 'lnode' is the input L-node,
    //!  - 'fnode' is the first F-node whose head is covered by 'lnode',
    //!  - 'lofst' means that 'lnode[lofst]' overlaps the back of 'fnode', and
    //!  - 'lwght' means that 'lwght's F-nodes from 'fnode' are covered by 'lnode'.
    //!
    //!     e.g.,
    //!                     [F]      [L]
    //!                  v1  a
    //!                      a        b  u1 (in)
    //!                * v2  a        b
    //!                      a  <---  b           The output is (u1, v2, 2, 3)
    //!                * v3  b        b           since u1[2] overlaps v2.back.
    //!                      b        b
    //!                * v4  c        b
    //!                      c
    //!
    lcover_type get_lcover(lnode_type* lnode) const {
        auto [fnode, lofst_e] = get_first_overlapped_fnode(lnode);
        offset_type lofst_b = lofst_e - get_fexp(fnode) + 1;

        // Is the first overlapped F-node not covered? (i.e., the head is not overlapped)
        if (lofst_b < 0) {
            fnode = fnode->get_next();
            lofst_b = lofst_e + 1;
            lofst_e = lofst_b + get_fexp(fnode) - 1;
        }

        // tfm::printfln("LCOVER: lnode=%s -> fnode=%s, lofst_b=%d, lofst_e=%d",  //
        //               get_pc(lnode), get_pc(fnode), lofst_b, lofst_e);

        lcover_type ret = {lnode, fnode, lofst_e, 0};
        const offset_type lexp = get_lexp(lnode);

        while (lofst_b < lexp) {
            ret.lwght += 1;
            lofst_b += get_fexp(fnode);
            fnode = fnode->get_next();
            if (is_head_fnode(fnode)) {  // cycled?
                break;
            }
        }
        return ret;
    }

    //! Return fcover_type(fnode, lnode, fofst, weight) such that
    //!  - 'fnode' is the input F-node,
    //!  - 'lnode' is the first L-node whose head is covered by 'fnode',
    //!  - 'fofst' means that 'fnode[fofst]' overlaps the back of 'lnode', and
    //!  - 'weight' means that 'weight's L-nodes from 'lnode' are covered by 'fnode'.
    //!
    //!     e.g.,
    //!                     [F]      [L]
    //!                               b  u1
    //!             (in) v1  a        b
    //!                      a        a  u2 *
    //!                      a        a         The output is (v1, u2, 3, 3)
    //!                      a  --->  a           since v1[3] overlaps u2.back.
    //!                      a        b  u3 *
    //!                      a        c  u4 *
    //!                      a        c
    //!
    fcover_type get_fcover(fnode_type* fnode) const {
        auto [lnode, fofst_e] = get_first_overlapped_lnode(fnode);
        offset_type fofst_b = fofst_e - get_lexp(lnode) + 1;

        // Is the first overlapped L-node not covered? (i.e., the head is not overlapped)
        if (fofst_b < 0) {
            lnode = lnode->get_next();
            fofst_b = fofst_e + 1;
            fofst_e = fofst_b + get_lexp(lnode) - 1;
        }

        fcover_type ret = {fnode, lnode, fofst_e, 0};
        const offset_type fexp = get_fexp(fnode);

        while (fofst_b < fexp) {
            ret.fwght += 1;
            fofst_b += get_lexp(lnode);
            lnode = lnode->get_next();
            if (is_dmmy_lnode(lnode)) {  // cycled?
                break;
            }
        }
        return ret;
    }

    /* * * * * * * * * * * * * * * * * * * * * * * *
     *
     *          Updaters for linkages
     *
     * * * * * * * * * * * * * * * * * * * * * * * */

    //! Remove the H/T links associated with the given L-nodes.
    void remove_lnode_links(lbuffer_type& lnodes) {
        for (lnode_type* lnode : lnodes) {
            remove_lnode_link(lnode);
        }
    }

    //! Remove the H/T links associated with the given F-nodes.
    void remove_fnode_links(fbuffer_type& fnodes) {
        for (fnode_type* fnode : fnodes) {
            remove_fnode_link(fnode);
        }
    }

    //! Remove the H/T links of the given L-node and the corresponding F-node.
    void remove_lnode_link(lnode_type* lnode) {
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

    //! Remove the H/T links of the given F-node and the corresponding L-node.
    void remove_fnode_link(fnode_type* fnode) {
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

    //! Reset the H/T links and weights associated with the given L-nodes.
    void reset_lnode_links_and_weights(lbuffer_type& lnodes) {
        for (lnode_type* lnode : lnodes) {
            reset_lnode_link_and_weight(lnode);
        }
    }

    //! Reset the H/T links and weights associated with the given F-nodes.
    void reset_fnode_links_and_weights(fbuffer_type& fnodes) {
        for (fnode_type* fnode : fnodes) {
            reset_fnode_link_and_weight(fnode);
        }
    }

    //! Reset the H/T links and weight associated with the given L-node.
    void reset_lnode_link_and_weight(lnode_type* lnode) {
        reset_lnode_link(lnode);
        reset_lnode_weight(lnode);
    }

    //! Reset the H/T links and weight associated with the given F-node.
    void reset_fnode_link_and_weight(fnode_type* fnode) {
        reset_fnode_link(fnode);
        reset_fnode_weight(fnode);
    }

    //! Reset the H/T links associated with the given L-node.
    void reset_lnode_link(lnode_type* lnode) {
        DEBUG_ABORT_IF(is_dmmy_lnode(lnode));

        if (lnode->get_data().get_hnode()) {
            return;
        }

        // Consider 'lnode' as the pivot, the L-interval is [lofst_lb, lofst_le].
        const offset_type lofst_lb = 0;
        const offset_type lofst_le = get_lexp(lnode) - 1;

        // Consider 'lnode' as the pivot, the F-interval is [lofst_fb, lofst_fe].
        auto [fnode, lofst_fe] = get_first_overlapped_fnode(lnode);
        offset_type lofst_fb = lofst_fe - get_fexp(fnode) + 1;

        if ((lofst_fb <= lofst_lb) and (lofst_lb <= lofst_fe) and (lofst_fe <= lofst_le)) {
            DEBUG_ABORT_IF(fnode->get_data().get_tnode());
            lnode->get_data().set_hlink(fnode, lofst_fe - lofst_lb);
            fnode->get_data().set_tnode(lnode);
        }
    }

    //! Reset the H/T links associated with the given F-node.
    void reset_fnode_link(fnode_type* fnode) {
        if (fnode->get_data().get_tnode()) {
            return;
        }

        // Consider 'fnode' as the pivot, the F-interval is [fofst_fb, fofst_fe].
        const offset_type fofst_fb = 0;
        const offset_type fofst_fe = get_fexp(fnode) - 1;

        // Consider 'fnode' as the pivot, the L-interval is [fofst_lb, fofst_le].
        auto [lnode, fofst_le] = get_first_overlapped_lnode(fnode);
        offset_type fofst_lb = fofst_le - get_lexp(lnode) + 1;

        while (fofst_lb <= fofst_fe) {
            if ((fofst_fb <= fofst_lb) and (fofst_lb <= fofst_fe) and (fofst_fe <= fofst_le)) {
                DEBUG_ABORT_IF(lnode->get_data().get_hnode());
                lnode->get_data().set_hlink(fnode, fofst_fe - fofst_lb);
                fnode->get_data().set_tnode(lnode);
                break;
            }
            lnode = lnode->get_next();
            fofst_lb = fofst_le + 1;
            fofst_le = fofst_le + get_lexp(lnode);
        }
    }

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

    //! Reset the weight associated with the given L-node.
    void reset_lnode_weight(lnode_type* lnode) {
        const lcover_type lcover = get_lcover(lnode);
        lnode->get_data().set_weight(lcover.lwght);
        push_fat_lnode(lnode);
    }

    //! Reset the weight associated with the given F-node.
    void reset_fnode_weight(fnode_type* fnode) {
        const fcover_type fcover = get_fcover(fnode);
        fnode->get_data().set_weight(fcover.fwght);
        push_fat_fnode(fnode);
    }

#ifdef ENABLE_DEBUG_PRINT
    void debug_print() {
        tfm::printfln("");
        tfm::printfln("**** LFIntervalGraph::debug_print() ****");
        tfm::printfln("");

        tfm::printfln("=== $ ===");
        tfm::printfln("m_em_lnode=%s, m_em_lofst=%d", get_pc(m_em_lnode), m_em_lofst);
        tfm::printfln("m_em_fnode=%s, m_em_fofst=%d", get_pc(m_em_fnode), m_em_fofst);

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
        return tfm::format("%s => %s, F=%s, H=%s, O=%d, W=%d, SA=%d",  //
                           get_pc(p),  //
                           make_run(p->get_data().get_chr(), get_lexp(p)),  //
                           get_pc(p->get_data().get_fnode()),  //
                           get_pc(p->get_data().get_hnode()),  //
                           p->get_data().get_hofst(),  //
                           p->get_data().get_weight(),  //
                           p->get_data().get_sae());
    }

    std::string get_fnode_str(const fnode_type* p) {
        return tfm::format("%s => %s, L=%s, T=%s, W=%d",  //
                           get_pc(p),  //
                           make_run(p->get_data().get_chr(), get_fexp(p)),  //
                           get_pc(p->get_data().get_lnode()),  //
                           get_pc(p->get_data().get_tnode()),  //
                           p->get_data().get_weight());
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