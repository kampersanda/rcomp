/**
 * @file Compressor_LFIG.hpp
 * @author Shunsuke Kanda (kampersanda)
 */
#pragma once

#include <array>
#include <map>

#include "basics.hpp"
#include "defs.hpp"
#include "utils.hpp"

#include "StaticQueue.hpp"
#include "StaticVector.hpp"

#ifdef ENABLE_STAT_MONITOR
#include "Timer.hpp"
#endif

namespace rcomp {

/**
 * A class for online RLBWT compressor with the straightforward LF-interval graph.
 *
 * @tparam t_LIndex Index for L-nodes
 * @tparam t_FIndex Index for F-nodes
 * @tparam t_DivBound The upper bound of dividing overlapped nodes
 */
template <class t_LIndex, class t_FIndex, size_type t_DivBound = 7>
class Compressor_LFIG {
  public:
    using this_type = Compressor_LFIG<t_LIndex, t_FIndex, t_DivBound>;

    using lindex_type = t_LIndex;
    using findex_type = t_FIndex;

    using lnode_type = typename lindex_type::lnode_type;
    using fnode_type = typename findex_type::fnode_type;

    // 3 is enough, but 4 is used for faster modulo
    using lqueue_type = StaticQueue<lnode_type*, 4>;
    using fqueue_type = StaticQueue<fnode_type*, 4>;

    using lbuffer_type = StaticVector<lnode_type*, t_DivBound + 2>;
    using fbuffer_type = StaticVector<fnode_type*, t_DivBound + 2>;

    static constexpr size_type DIV_BOUND = t_DivBound;
    static_assert(7 <= DIV_BOUND, "The division bound must be no less than 7");

  private:
    enum class insert_modes : uint8_t {
        // Case C
        MERGE_PREV,
        MERGE_CURR,
        MERGE_NEXT,
        // Case B
        SPLIT_FRONT,
        SPLIT_MIDDLE,
        SPLIT_BACK,
    };

    struct lcover_type {
        lnode_type* lnode;
        fnode_type* fnode;
        offset_type lofst;
        size_type   weight;
    };

    struct fcover_type {
        fnode_type* fnode;
        lnode_type* lnode;
        offset_type fofst;
        size_type   weight;
    };

    // Data Structure
    lindex_type m_lindex;
    findex_type m_findex;

    // i.e., original text size (with $)
    size_type m_num_chars = 1;

    // Result of $-query
    lnode_type* m_em_lnode = nullptr;
    offset_type m_em_lofst = -1;

    // The corresponding $'s position in F
    fnode_type* m_em_fnode = nullptr;
    offset_type m_em_fofst = -1;

    lqueue_type m_fat_lnodes;
    fqueue_type m_fat_fnodes;

    lbuffer_type m_lnodes_buffer;
    fbuffer_type m_fnodes_buffer;

#ifdef ENABLE_DEBUG_PRINT
    std::map<intptr_t, char> m_pcmap = {{intptr_t(nullptr), '?'}};
    char                     m_pcmax = 'A';
#endif

#ifdef ENABLE_STAT_MONITOR
    Timer m_timer_4_case;
    Timer m_timer_4_phase;

    size_type m_num_case_C  = 0;
    size_type m_num_case_B  = 0;
    size_type m_num_case_LF = 0;

    double m_ns_case_C  = 0.0;
    double m_ns_case_B  = 0.0;
    double m_ns_case_LF = 0.0;

    double m_ns_phase_B_in_C = 0.0;
    double m_ns_phase_C_in_C = 0.0;
    double m_ns_phase_D_in_C = 0.0;
    double m_ns_phase_E_in_C = 0.0;
    double m_ns_phase_F_in_C = 0.0;

    double m_ns_phase_B_in_B = 0.0;
    double m_ns_phase_C_in_B = 0.0;
    double m_ns_phase_D_in_B = 0.0;
    double m_ns_phase_E_in_B = 0.0;
    double m_ns_phase_F_in_B = 0.0;

    size_type m_num_merge_curr = 0;
    size_type m_num_merge_prev = 0;
    size_type m_num_merge_next = 0;
#endif

  public:
    //! Default constructor
    Compressor_LFIG() = default;

    //! Default destructor
    virtual ~Compressor_LFIG() = default;

    //! Copy constructor (deleted)
    Compressor_LFIG(const Compressor_LFIG&) = delete;

    //! Copy constructor (deleted)
    Compressor_LFIG& operator=(const Compressor_LFIG&) = delete;

    //! Move constructor
    Compressor_LFIG(Compressor_LFIG&&) noexcept = default;

    //! Move constructor
    Compressor_LFIG& operator=(Compressor_LFIG&&) noexcept = default;

    //! Get if the data structure is empty.
    bool is_empty() const {
        return m_em_lnode == nullptr;
    }

    /**
     * @brief Extend the RLBWT text by appending one character (i.e., \f$ T := c + T \f$).
     * @param[in] new_chr Character to be appended
     */
    void extend(const uchar_type new_chr) {
        ABORT_IF_LE(new_chr, END_MARKER);

        m_num_chars += 1;

        DEBUG_PRINT(tfm::printfln("");)
        DEBUG_PRINT(tfm::printfln("*********************************");)
        DEBUG_PRINT(tfm::printfln("     extend(chr=%c, size=%d)", new_chr, m_num_chars);)
        DEBUG_PRINT(tfm::printfln("*********************************");)
        DEBUG_PRINT(tfm::printfln("");)

        if (is_empty()) {
            extend_init(new_chr);
            DEBUG_PRINT(debug_print();)
            DEBUG_PRINT(test_all();)
            return;
        }

        // Case L or F
        STAT_MONITOR(m_timer_4_case.start();)
        [[maybe_unused]] const size_type num_divided = divide_covers();
        STAT_MONITOR(m_ns_case_LF += m_timer_4_case.stop_and_get_ns();)
        STAT_MONITOR(m_num_case_LF += num_divided;)

        // Phase A
        STAT_MONITOR(m_timer_4_case.start();)
        const insert_modes ins_mode = check_joinable(new_chr);

        // Phase B--F
        if (is_merge_mode(ins_mode)) {
            extend_C(new_chr, ins_mode);
            STAT_MONITOR(m_ns_case_C += m_timer_4_case.stop_and_get_ns(););
            STAT_MONITOR(m_num_case_C += 1;)
        } else {
            extend_B(new_chr, ins_mode);
            STAT_MONITOR(m_ns_case_B += m_timer_4_case.stop_and_get_ns();)
            STAT_MONITOR(m_num_case_B += 1;)
        }

        DEBUG_PRINT(debug_print();)
        DEBUG_PRINT(test_all();)
    }

    /**
     * @brief Output the RLBWT text.
     * @param[in] fn Callback function to get the RLBWT text run by run.
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

        for (const lnode_type* lnode = get_head_lnode();; lnode = lnode->get_next()) {
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

    /**
     * @brief Output the original text.
     * @param[in] fn Callback function to get the original text charcater by charcater.
     */
    void decode(const std::function<void(uchar_type)>& fn) const {
        if (is_empty()) {
            return;
        }

        const lnode_type* lnode = m_em_lnode;
        offset_type       lofst = m_em_lofst;

        while (true) {
            const auto [fnode, fofst] = get_fpartner(lnode, lofst);

            if (is_head_fnode(fnode)) {
                break;
            }
            fn(fnode->get_data().get_chr());

            lnode = fnode->get_data().get_lnode();
            lofst = fofst;

            if (lnode == m_em_lnode and m_em_lofst <= lofst) {
                lofst += 1;
            }
        }
    }

    //! Get the number of characters in the text.
    size_type get_num_chars() const {
        return m_num_chars;
    }

    //! Get the allocated memory in bytes.
    size_type get_memory_in_bytes() const {
        size_type mem = 0;
        mem += m_lindex.get_memory_in_bytes();
        mem += m_findex.get_memory_in_bytes();
        mem += sizeof(m_num_chars);
        mem += sizeof(m_em_lnode);
        mem += sizeof(m_em_lofst);
        mem += sizeof(m_em_fnode);
        mem += sizeof(m_em_fofst);
        mem += sizeof(m_fat_lnodes);
        mem += sizeof(m_fat_fnodes);
        mem += sizeof(m_lnodes_buffer);
        mem += sizeof(m_fnodes_buffer);
        return mem;
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

    //! Print the statistics measured with internal monitors.
    void show_monitored_statistics() const {
#ifdef ENABLE_STAT_MONITOR
        tfm::printfln("[Monitor_Case]");
        tfm::reportfln("num_case_C:\t%d", m_num_case_C);
        tfm::reportfln("num_case_B:\t%d", m_num_case_B);
        tfm::reportfln("num_case_LF:\t%d", m_num_case_LF);
        tfm::reportfln("sec_case_C:\t%g", utils::ns_to_sec(m_ns_case_C));
        tfm::reportfln("sec_case_B:\t%g", utils::ns_to_sec(m_ns_case_B));
        tfm::reportfln("sec_case_LF:\t%g", utils::ns_to_sec(m_ns_case_LF));
        tfm::reportfln("ave_ns_in_case_C:\t%g", m_ns_case_C / m_num_case_C);
        tfm::reportfln("ave_ns_in_case_B:\t%g", m_ns_case_B / m_num_case_B);
        tfm::reportfln("ave_ns_in_case_LF:\t%g", m_ns_case_LF / m_num_case_LF);

        tfm::printfln("[Monitor_PhaseInCaseC]");
        tfm::reportfln("sec_phase_B:\t%g", utils::ns_to_sec(m_ns_phase_B_in_C));
        tfm::reportfln("sec_phase_C:\t%g", utils::ns_to_sec(m_ns_phase_C_in_C));
        tfm::reportfln("sec_phase_D:\t%g", utils::ns_to_sec(m_ns_phase_D_in_C));
        tfm::reportfln("sec_phase_E:\t%g", utils::ns_to_sec(m_ns_phase_E_in_C));
        tfm::reportfln("sec_phase_F:\t%g", utils::ns_to_sec(m_ns_phase_F_in_C));

        tfm::printfln("[Monitor_PhaseInCaseB]");
        tfm::reportfln("sec_phase_B:\t%g", utils::ns_to_sec(m_ns_phase_B_in_B));
        tfm::reportfln("sec_phase_C:\t%g", utils::ns_to_sec(m_ns_phase_C_in_B));
        tfm::reportfln("sec_phase_D:\t%g", utils::ns_to_sec(m_ns_phase_D_in_B));
        tfm::reportfln("sec_phase_E:\t%g", utils::ns_to_sec(m_ns_phase_E_in_B));
        tfm::reportfln("sec_phase_F:\t%g", utils::ns_to_sec(m_ns_phase_F_in_B));

        tfm::printfln("[Monitor_MergeInCaseC]");
        tfm::reportfln("num_merge_curr:\t%d", m_num_merge_curr);
        tfm::reportfln("num_merge_prev:\t%d", m_num_merge_prev);
        tfm::reportfln("num_merge_next:\t%d", m_num_merge_next);
#endif
    }

    //! Test the data structure.
    void test_all() const {
        test_alphabet_order();
        test_list_order();
        test_lf_mapping();
        test_ht_links();
        test_weights();
    }

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
                fnode  = fnode->get_next();
                lnode  = lnode->get_next();
            } else if (fpos_e < lpos_e) {
                ABORT_IF(fnode->get_data().get_tnode());
                fpos_b = fpos_e + 1;
                fnode  = fnode->get_next();
            } else if (fpos_e > lpos_e) {
                ABORT_IF(lnode->get_data().get_hnode());
                lpos_b = lpos_e + 1;
                lnode  = lnode->get_next();
            } else {
                ABORT_IF(true);
            }
        }

        while (!is_head_fnode(fnode)) {
            ABORT_IF(fnode->get_data().get_tnode());
            fpos_b = fpos_b + get_fexp(fnode);
            fnode  = fnode->get_next();
        }

        ABORT_IF_NE(fpos_b, lpos_b);
        ABORT_IF_NE(fpos_b, offset_type(m_num_chars));
    }

    //! Test the weights of L/F-nodes.
    void test_weights() const {
        for (fnode_type* fnode = get_head_fnode();; fnode = fnode->get_next()) {
            const fcover_type fcover = get_fcover(fnode);
            ABORT_IF_NE(fcover.weight, fnode->get_data().get_weight());
            if (fnode->get_data().get_weight() >= DIV_BOUND) {
                ABORT_IF(!m_fat_fnodes.is_member(fnode));
            }
            if (fnode == get_tail_fnode()) {
                break;
            }
        }

        for (lnode_type* lnode = get_head_lnode();; lnode = lnode->get_next()) {
            const lcover_type lcover = get_lcover(lnode);
            ABORT_IF_NE(lcover.weight, lnode->get_data().get_weight());
            if (lnode->get_data().get_weight() >= DIV_BOUND) {
                ABORT_IF(!m_fat_lnodes.is_member(lnode));
            }
            if (lnode == get_tail_lnode()) {
                break;
            }
        }
    }

#ifdef ENABLE_DEBUG_PRINT
    void debug_print() {
        tfm::printfln("");
        tfm::printfln("**** debug_print() ****");
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
#endif

  private:
    bool is_head_fnode(const fnode_type* fnode) const {
        return fnode == m_findex.get_head();
    }
    bool is_dmmy_lnode(const lnode_type* lnode) const {
        return lnode == m_lindex.get_head();
    }

    fnode_type* get_head_fnode() const {
        return m_findex.get_head();
    }
    lnode_type* get_head_lnode() const {
        return m_lindex.get_head()->get_next();
    }
    fnode_type* get_tail_fnode() const {
        return m_findex.get_head()->get_prev();
    }
    lnode_type* get_tail_lnode() const {
        return m_lindex.get_head()->get_prev();
    }

    offset_type get_fexp(const fnode_type* fnode) const {
        return fnode->get_data().get_exp();
    }
    offset_type get_lexp(const lnode_type* lnode) const {
        const offset_type lexp = lnode->get_data().get_exp();
        return lnode == m_em_lnode ? lexp + 1 : lexp;
    }

    run_type get_frun(const fnode_type* fnode) const {
        return run_type{fnode->get_data().get_chr(), size_type(get_fexp(fnode))};
    }
    run_type get_lrun(const lnode_type* lnode) const {
        return run_type{lnode->get_data().get_chr(), size_type(get_lexp(lnode))};
    }

    inline bool is_first_lofst(const lnode_type*, offset_type lofst) const {
        return lofst == 0;
    }
    inline bool is_last_lofst(const lnode_type* lnode, offset_type lofst) const {
        return lnode->get_data().get_exp() == size_type(lofst + 1);
    }
    inline bool is_first_fofst(const fnode_type*, offset_type fofst) const {
        return fofst == 0;
    }
    inline bool is_last_fofst(const fnode_type* fnode, offset_type fofst) const {
        return fnode->get_data().get_exp() == size_type(fofst + 1);
    }

    bool is_merge_mode(insert_modes ins_mode) {
        return ins_mode == insert_modes::MERGE_PREV ||  //
               ins_mode == insert_modes::MERGE_CURR ||  //
               ins_mode == insert_modes::MERGE_NEXT;
    }

    // The distance between the begin of the given F-node and that of the linked L-node.
    offset_type get_begin_offset(const fnode_type* fnode) const {
        DEBUG_ABORT_IF(!fnode->get_data().get_tnode());
        DEBUG_ABORT_IF_LT(0, fnode->get_data().get_tofst() - get_fexp(fnode) + 1);
        return -1 * (fnode->get_data().get_tofst() - get_fexp(fnode) + 1);
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

    //! Check if $'s L-node can be merged with the character.
    insert_modes check_joinable(const uchar_type new_chr) const {
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
        } else if (m_em_lofst == get_lexp(m_em_lnode) - 1) {  // -1 is for $
            // Mergeable with the next node?
            if (new_chr == m_em_lnode->get_next()->get_data().get_chr()) {
                return insert_modes::MERGE_NEXT;  // Case C
            }
        }

        if (m_em_lofst == 0) {
            return insert_modes::SPLIT_FRONT;
        } else if (m_em_lofst == get_lexp(m_em_lnode) - 1) {  // -1 is for $
            ABORT_IF(m_em_lnode != get_tail_lnode());
            return insert_modes::SPLIT_BACK;
        } else {
            return insert_modes::SPLIT_MIDDLE;
        }
    }

    //! Make the data structure for the RLBWT text of T = (new_chr,$).
    void extend_init(uchar_type new_chr) {
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
    }

    //! Extend the RLBWT text in Case C.
    void extend_C(const uchar_type, const insert_modes ins_mode) {
        DEBUG_ABORT_IF(!is_merge_mode(ins_mode));
        DEBUG_PRINT(tfm::printfln("\n==== Case C ====");)

        DEBUG_ABORT_IF(!m_fat_lnodes.is_empty());
        DEBUG_ABORT_IF(!m_fat_fnodes.is_empty());

        lnode_type* ins_lnode = nullptr;
        fnode_type* ins_fnode = nullptr;
        offset_type ins_fofst = 0;

        // NOTE:
        // Depending on 'ins_mode', the L-node with $ can be updated, but Do NOT change the value of 'm_em_lnode' until
        // Phase E. This is because this change will affect 'get_lexp'.

        //
        // Phase B
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase B ==");)
        STAT_MONITOR(m_timer_4_phase.start();)

        if (ins_mode == insert_modes::MERGE_CURR) {
            // Then, the new char will be inserted into 'm_em_lnode'
            DEBUG_ABORT_IF_OUT(m_em_lofst, 0, get_lexp(m_em_lnode) - 1);
            // L (the overlapped F-nodes do not need to be updated)
            ins_lnode = m_em_lnode;
            // F
            ins_fnode = ins_lnode->get_data().get_fnode();
            ins_fofst = m_em_lofst;
            STAT_MONITOR(m_num_merge_curr += 1;)
        } else if (ins_mode == insert_modes::MERGE_PREV) {
            // Then, the new char will be inserted into m_em_lnode->get_prev()
            DEBUG_ABORT_IF_NE(m_em_lofst, 0);
            // L
            ins_lnode = m_em_lnode->get_prev();
            // F
            ins_fnode = ins_lnode->get_data().get_fnode();
            ins_fofst = get_fexp(ins_fnode);
            STAT_MONITOR(m_num_merge_prev += 1;)
        } else {
            ABORT_IF(true);
        }
        STAT_MONITOR(m_ns_phase_B_in_C += m_timer_4_phase.stop_and_get_ns();)

        //
        // Phase C (Remove H/T-links of overlapped nodes)
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase C ==");)
        STAT_MONITOR(m_timer_4_phase.start();)

        // $-replacement
        if (ins_mode == insert_modes::MERGE_PREV) {
            lnode_type* em_tnode = m_em_fnode->get_data().get_tnode();

            const bool is_first_em_fofst = is_first_fofst(m_em_fnode, m_em_fofst);
            const bool is_last_em_fofst  = is_last_fofst(m_em_fnode, m_em_fofst);

            // Update: OV
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

            // Update the linkage
            if (is_last_em_fofst) {
                DEBUG_ABORT_IF(em_tnode != m_em_lnode);
                DEBUG_ABORT_IF_NE(em_tnode->get_data().get_hofst(), 0);
                {
                    fnode_type*       fnode_befo = m_em_fnode;
                    lnode_type*       lnode_befo = m_em_lnode->get_prev();
                    const offset_type new_fexp   = fnode_befo->get_data().get_exp();
                    const offset_type new_lexp   = lnode_befo->get_data().get_exp() + 1;
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
                    fnode_type*       fnode_aftr = m_em_fnode->get_next();
                    lnode_type*       lnode_aftr = m_em_lnode;
                    const offset_type new_fexp   = fnode_aftr->get_data().get_exp();
                    const offset_type new_lexp   = lnode_aftr->get_data().get_exp();
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
        STAT_MONITOR(m_ns_phase_C_in_C += m_timer_4_phase.stop_and_get_ns();)

        //
        // Phase D
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase D ==");)
        STAT_MONITOR(m_timer_4_phase.start();)

        ins_lnode->get_data().add_exp(1);

        STAT_MONITOR(m_ns_phase_D_in_C += m_timer_4_phase.stop_and_get_ns();)
        DEBUG_PRINT(tfm::printfln("ins_lnode=%s", get_lnode_str(ins_lnode));)

        //
        // Phase E
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase E ==");)
        STAT_MONITOR(m_timer_4_phase.start();)

        // Clear since the L-run has already updated.
        m_em_lnode = nullptr;
        m_em_lofst = 0;

        // Cache to the next process
        m_em_fnode = ins_fnode;
        m_em_fofst = ins_fofst;

        std::tie(m_em_lnode, m_em_lofst) = get_lpartner(ins_fnode, ins_fofst);

        STAT_MONITOR(m_ns_phase_E_in_C += m_timer_4_phase.stop_and_get_ns();)
        DEBUG_PRINT(tfm::printfln("m_em_lnode=%s", get_lnode_str(m_em_lnode));)
        DEBUG_PRINT(tfm::printfln("m_em_lofst=%d", m_em_lofst);)

        //
        // Phase F
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase F ==");)
        STAT_MONITOR(m_timer_4_phase.start();)

        // ?-replacement
        {
            const bool is_last_em_fofst  = is_last_fofst(m_em_fnode, m_em_fofst);
            const bool is_first_em_lofst = is_first_lofst(m_em_lnode, m_em_lofst);

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

        STAT_MONITOR(m_ns_phase_F_in_C += m_timer_4_phase.stop_and_get_ns();)
    }

    // Extend the BWT text in Case B.
    void extend_B(const uchar_type new_chr, const insert_modes ins_mode) {
        DEBUG_ABORT_IF(is_merge_mode(ins_mode));
        DEBUG_PRINT(tfm::printfln("\n==== Case B ====");)

        DEBUG_ABORT_IF(!m_fat_lnodes.is_empty());
        DEBUG_ABORT_IF(!m_fat_fnodes.is_empty());

        //
        // Phase B
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase B ==");)
        STAT_MONITOR(m_timer_4_phase.start();)

        lnode_type* lpred = ins_mode == insert_modes::SPLIT_FRONT ? m_em_lnode->get_prev() : m_em_lnode;
        fnode_type* fpred = m_findex.predecessor(new_chr, lpred->get_data().get_fnode());

        STAT_MONITOR(m_ns_phase_B_in_B += m_timer_4_phase.stop_and_get_ns();)
        DEBUG_PRINT(tfm::printfln("lpred=%s, fpred=%s", get_pc(lpred), get_pc(fpred));)

        //
        // Phase C -- F (for division)
        //
        if (ins_mode == insert_modes::SPLIT_MIDDLE) {
            fnode_type* div_fnode = divide_lnode(m_em_lnode, m_em_lofst);
            if (fpred == div_fnode) {
                fpred = div_fnode->get_next();
            }
        }

        //
        // Phase C (for insert)
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase C =="););

        lnode_type* ins_lnode = m_em_lnode;
        fnode_type* ins_fnode = get_tail_fnode() != fpred ? fpred->get_next() : fpred;
        get_overlapped_fnodes(ins_lnode, 0, get_lexp(ins_lnode), m_fnodes_buffer);
        get_overlapped_lnodes(ins_fnode, 0, get_fexp(ins_fnode), m_lnodes_buffer);
        remove_fnode_links(m_fnodes_buffer);
        remove_lnode_links(m_lnodes_buffer);

        DEBUG_PRINT(tfm::printf("m_fnodes_buffer:");)
        DEBUG_PRINT(for (const auto* p : m_fnodes_buffer) tfm::printf(" %s", get_pc(p));)
        DEBUG_PRINT(tfm::printf("\n");)
        DEBUG_PRINT(tfm::printf("m_lnodes_buffer:");)
        DEBUG_PRINT(for (const auto* p : m_lnodes_buffer) tfm::printf(" %s", get_pc(p));)
        DEBUG_PRINT(tfm::printf("\n");)

        //
        // Phase D (for insert)
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase D ==");)
        STAT_MONITOR(m_timer_4_phase.start();)

        lnode_type* new_lnode = m_lindex.insert_after(lpred, new_chr, 1);
        fnode_type* new_fnode = m_findex.insert_after(fpred, new_lnode);

        STAT_MONITOR(m_ns_phase_D_in_B += m_timer_4_phase.stop_and_get_ns();)
        DEBUG_PRINT(tfm::printfln("new_lnode =%s", get_lnode_str(new_lnode));)
        DEBUG_PRINT(tfm::printfln("new_fnode =%s", get_fnode_str(new_fnode));)

        //
        // Phase E
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase E ==");)
        STAT_MONITOR(m_timer_4_phase.start();)

        // Clear since the L-run has already updated.
        m_em_lnode = nullptr;
        m_em_lofst = 0;

        // Cache to the next process
        m_em_fnode = new_fnode;
        m_em_fofst = 0;

        std::tie(m_em_lnode, m_em_lofst) = get_lpartner(new_fnode, 0);

        STAT_MONITOR(m_ns_phase_E_in_B += m_timer_4_phase.stop_and_get_ns();)
        DEBUG_PRINT(tfm::printfln("m_em_lnode=%s", get_lnode_str(m_em_lnode));)
        DEBUG_PRINT(tfm::printfln("m_em_lofst=%d", m_em_lofst);)

        //
        // Phase F
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase F ==");)
        STAT_MONITOR(m_timer_4_phase.start();)

        reset_lnode_links(m_lnodes_buffer);
        reset_fnode_links(m_fnodes_buffer);
        reset_lnode_link(ins_lnode);
        reset_lnode_link(new_lnode);
        reset_fnode_link(ins_fnode);
        reset_fnode_link(new_fnode);

        // Release
        m_fnodes_buffer.clear();
        m_lnodes_buffer.clear();

        STAT_MONITOR(m_ns_phase_F_in_B += m_timer_4_phase.stop_and_get_ns();)
    }

    // cover-queryの結果がdivedeで変わる場合がありので、
    // on-the-flyでやるひつようあり
    size_type divide_covers() {
        size_type num_divided = 0;

        while (!m_fat_lnodes.is_empty()) {
            lnode_type* lnode = m_fat_lnodes.pop();
            DEBUG_ABORT_IF(!lnode);

            if (lnode->get_data().get_weight() >= DIV_BOUND) {
                const lcover_type lcover = get_lcover(lnode);
                DEBUG_ABORT_IF_NE(lcover.weight, lnode->get_data().get_weight());

                divide_lcover(lcover);
                num_divided += 1;
            }
        }

        while (!m_fat_fnodes.is_empty()) {
            fnode_type* fnode = m_fat_fnodes.pop();
            DEBUG_ABORT_IF(!fnode);

            if (fnode->get_data().get_weight() >= DIV_BOUND) {
                const fcover_type fcover = get_fcover(fnode);
                DEBUG_ABORT_IF_NE(fcover.weight, fnode->get_data().get_weight());

                divide_fcover(fcover);
                num_divided += 1;
            }
        }

        return num_divided;
    }

    // Runを割るのはSplitで統一したほうがよさげ
    // Groupはdiv
    void divide_lcover(const lcover_type& lcover) {
        DEBUG_PRINT(tfm::printfln("\n==== Case L ====");)
        DEBUG_PRINT(tfm::printfln("lcover => %d", get_lcover_str(lcover));)

        // Phase B
        lnode_type* div_lnode = lcover.lnode;  // to be split
        offset_type div_lofst = lcover.lofst;
        {
            const size_type   fnum1 = lcover.weight / 2;
            const fnode_type* fnode = lcover.fnode;
            for (size_type i = 1; i < fnum1; i++) {
                fnode = fnode->get_next();
                div_lofst += get_fexp(fnode);
            }
        }
        // Here, div_lnode[div_lofst] is the last character of the front part

        // Phase C--F
        if (div_lnode == m_em_lnode and m_em_lofst <= div_lofst) {
            divide_lnode(div_lnode, div_lofst);
        } else {
            divide_lnode(div_lnode, div_lofst + 1);
        }
    }

    void divide_fcover(const fcover_type& fcover) {
        DEBUG_PRINT(tfm::printfln("\n==== Case F ====");)
        DEBUG_PRINT(tfm::printfln("fcover => %d", get_fcover_str(fcover));)

        // Phase B
        fnode_type* div_fnode = fcover.fnode;  // to be split
        lnode_type* div_lnode = div_fnode->get_data().get_lnode();
        offset_type div_fofst = fcover.fofst;
        {
            const size_type   lnum1 = fcover.weight / 2;
            const lnode_type* lnode = fcover.lnode;
            for (size_type i = 1; i < lnum1; i++) {
                lnode = lnode->get_next();
                div_fofst += get_lexp(lnode);
            }
        }
        // Here, div_fnode[div_fofst] is the last character of the front part

        // Phase C--F
        divide_lnode(div_lnode, div_fofst + 1);
    }

    //! Divide the L-node with $-mark and the corresponding F-node.
    //! As the result, $-mark will be placed at the first of the after L-node.
    fnode_type* divide_lnode(lnode_type* div_lnode, const offset_type new_exp1) {
        DEBUG_PRINT(tfm::printfln("\n** divide_lnode **");)
        DEBUG_PRINT(tfm::printfln("div_lnode=%s, new_exp1=%d", get_pc(div_lnode), new_exp1);)

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

        DEBUG_PRINT(tfm::printf("m_fnodes_buffer:");)
        DEBUG_PRINT(for (const auto* p : m_fnodes_buffer) tfm::printf(" %s", get_pc(p));)
        DEBUG_PRINT(tfm::printf("\n");)
        DEBUG_PRINT(tfm::printf("m_lnodes_buffer:");)
        DEBUG_PRINT(for (const auto* p : m_lnodes_buffer) tfm::printf(" %s", get_pc(p));)
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

        reset_fnode_links(m_fnodes_buffer);
        reset_lnode_links(m_lnodes_buffer);
        reset_lnode_link(div_lnode);
        reset_lnode_link(new_lnode);
        reset_fnode_link(div_fnode);
        reset_fnode_link(new_fnode);

        // Release
        m_fnodes_buffer.clear();
        m_lnodes_buffer.clear();

        return div_fnode;
    }

    // Return (lnode, lofst) such that q_fnode[q_fofst] corresponds to lnode[lofst].
    inline std::pair<lnode_type*, offset_type> get_lpartner(const fnode_type* q_fnode, offset_type q_fofst) const {
        DEBUG_ABORT_IF(!q_fnode);
        DEBUG_ABORT_IF(is_head_fnode(q_fnode));

        // Exception case: \Lambda_t indicates the tail character of F?
        if (q_fnode == get_tail_fnode() and get_fexp(q_fnode) == q_fofst + 1) {
            // Then, 'lofst' to be returned indicates exclusive position (only in this case).
            lnode_type* ltail = get_tail_lnode();
            return {ltail, get_lexp(ltail)};
        }

        auto [lnode, fofst] = get_first_overlapped_lnode(q_fnode);

        while (fofst < q_fofst) {
            lnode = lnode->get_next();
            fofst += get_lexp(lnode);
        }
        return {lnode, (get_lexp(lnode) - 1) - (fofst - q_fofst)};
    }

    // Return (fnode, fofst) such that q_lnode[q_lofst] corresponds to fnode[fofst].
    inline std::pair<fnode_type*, offset_type> get_fpartner(const lnode_type* q_lnode, offset_type q_lofst) const {
        DEBUG_ABORT_IF(!q_lnode);
        DEBUG_ABORT_IF(is_dmmy_lnode(q_lnode));

        auto [fnode, lofst] = get_first_overlapped_fnode(q_lnode);

        while (lofst < q_lofst) {
            fnode = fnode->get_next();
            lofst += get_fexp(fnode);
        }
        return {fnode, (get_fexp(fnode) - 1) - (lofst - q_lofst)};
    }

    //  Return (lnode, fofst) such that fnode[fofst] corresponds to the back of lnode.
    //
    //    e.g.,
    //               [F]      [L]
    //            v1  a
    //                a        b  u1
    //            v2  a        b
    //                a        c  u2 *      The output is (u2, 1)
    //       (in) v3  b        c            since v3[1] corresponds to u2.back().
    //                b  --->  c
    //                b        b  u3
    //                         b
    //
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
        fofst += get_begin_offset(fnode);

        // Move to the first overlapped L-node with the given F-node
        do {
            fofst += get_lexp(lnode);
            lnode = lnode->get_next();
        } while (fofst <= 0);

        return {lnode->get_prev(), fofst - 1};
    }

    //  Return (fnode, lofst) such that lnode[lofst] corresponds to the back of fnode.
    //
    //     e.g.,
    //               [F]      [L]
    //            v1  a
    //                a        b  u1
    //          * v2  a        b
    //                a        c  u2 (in)   The output is (v2, 2)
    //                a        c            since u2[2] corresponds to v2.back.
    //                a  <---  c
    //            v3  b        b  u3
    //                b        b
    //
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
        lofst -= get_begin_offset(fnode);

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

    //  Return the L-nodes overlapped with the F-interval of size 'fsize' starting at the 'fofst'-th of 'fnode'.
    //
    //     e.g.,
    //                     [F]      [L]
    //                  v1  a
    //                      a        b  u1
    //             (in) v2  a        b
    //     fofst = 1 -->    a +    + c  u2 *       The output is [u2, u3].
    //                  v3  b |    | c
    //                      b |    | c
    //     fsize = 4 -->    b +    | b  u3 *
    //                             + b
    //
    inline void get_overlapped_lnodes(const fnode_type* fnode, offset_type fofst, const offset_type fsize,
                                      lbuffer_type& buf) {
        DEBUG_ABORT_IF(is_head_fnode(fnode));

        if (fsize <= 0) {
            return;
        }

        const offset_type fofst_b = fofst;
        const offset_type fofst_e = fofst + fsize - 1;

        lnode_type* lnode      = nullptr;
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

    //  Return the F-nodes overlapped with the L-interval of size 'lsize' starting at the 'lofst'-th of 'lnode'.
    //
    //     e.g.,
    //                 [F]      [L]
    //              v1  a                          The output is [v2, v3].
    //                  a        b  u1 (in)
    //            * v2  a +      b
    //                  a |    + c  u2  <-- lofst = 2
    //            * v3  b |    | c
    //                  b |    + c      <-- lsize = 3
    //                  b +      b  u3
    //                           b
    //
    inline void get_overlapped_fnodes(const lnode_type* lnode, offset_type lofst, const offset_type lsize,
                                      fbuffer_type& buf) {
        DEBUG_ABORT_IF(is_dmmy_lnode(lnode));

        if (lsize <= 0) {
            return;
        }

        const offset_type lofst_b = lofst;
        const offset_type lofst_e = lofst + lsize - 1;

        fnode_type* fnode      = nullptr;
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

    //  Return lcover_type(lnode, fnode, lofst, weight) such that
    //   - 'lnode' is the input L-node,
    //   - 'fnode' is the first F-node whose head is covered by 'lnode',
    //   - 'lofst' means that 'lnode[lofst]' corresponds to the back of 'fnode', and
    //   - 'weight' means that 'weight's F-nodes from 'fnode' are covered by 'lnode'.
    //
    //     e.g.,
    //                     [F]      [L]
    //                  v1  a
    //                      a        b  u1 (in)
    //                * v2  a        b
    //                      a  <---  b           The output is (u1, v2, 2, 3)
    //                * v3  b        b           since u1[2] corresponds to v2.back.
    //                      b        b
    //                * v4  c        b
    //                      c
    //
    lcover_type get_lcover(lnode_type* lnode) const {
        auto [fnode, lofst_e] = get_first_overlapped_fnode(lnode);
        offset_type lofst_b   = lofst_e - get_fexp(fnode) + 1;

        // Is the first overlapped F-node not covered? (i.e., the head is not overlapped)
        if (lofst_b < 0) {
            fnode   = fnode->get_next();
            lofst_b = lofst_e + 1;
            lofst_e = lofst_b + get_fexp(fnode) - 1;
        }

        // tfm::printfln("LCOVER: lnode=%s -> fnode=%s, lofst_b=%d, lofst_e=%d",  //
        //               get_pc(lnode), get_pc(fnode), lofst_b, lofst_e);

        lcover_type       ret  = {lnode, fnode, lofst_e, 0};
        const offset_type lexp = get_lexp(lnode);

        while (lofst_b < lexp) {
            ret.weight += 1;
            lofst_b += get_fexp(fnode);
            fnode = fnode->get_next();
            if (is_head_fnode(fnode)) {  // cycled?
                break;
            }
        }
        return ret;
    }

    //  Return fcover_type(fnode, lnode, fofst, weight) such that
    //   - 'fnode' is the input F-node,
    //   - 'lnode' is the first L-node whose head is covered by 'fnode',
    //   - 'fofst' means that 'fnode[fofst]' corresponds to the back of 'lnode', and
    //   - 'weight' means that 'weight's L-nodes from 'lnode' are covered by 'fnode'.
    //
    //     e.g.,
    //                     [F]      [L]
    //                               b  u1
    //             (in) v1  a        b
    //                      a        a  u2 *
    //                      a        a         The output is (v1, u2, 3, 3)
    //                      a  --->  a           since v1[3] corresponds to u2.back.
    //                      a        b  u3 *
    //                      a        c  u4 *
    //                      a        c
    //
    fcover_type get_fcover(fnode_type* fnode) const {
        auto [lnode, fofst_e] = get_first_overlapped_lnode(fnode);
        offset_type fofst_b   = fofst_e - get_lexp(lnode) + 1;

        // Is the first overlapped L-node not covered? (i.e., the head is not overlapped)
        if (fofst_b < 0) {
            lnode   = lnode->get_next();
            fofst_b = fofst_e + 1;
            fofst_e = fofst_b + get_lexp(lnode) - 1;
        }

        fcover_type       ret  = {fnode, lnode, fofst_e, 0};
        const offset_type fexp = get_fexp(fnode);

        while (fofst_b < fexp) {
            ret.weight += 1;
            fofst_b += get_lexp(lnode);
            lnode = lnode->get_next();
            if (is_dmmy_lnode(lnode)) {  // cycled?
                break;
            }
        }
        return ret;
    }

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

    void remove_lnode_links(lbuffer_type& lnodes) {
        for (lnode_type* lnode : lnodes) {
            remove_lnode_link(lnode);
        }
    }

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

    void remove_fnode_links(fbuffer_type& fnodes) {
        for (fnode_type* fnode : fnodes) {
            remove_fnode_link(fnode);
        }
    }

    void reset_lnode_link(lnode_type* lnode) {
        DEBUG_ABORT_IF(is_dmmy_lnode(lnode));

        if (lnode->get_data().get_hnode()) {
            reset_lweight(lnode);
            return;
        }

        // Consider 'lnode' as the pivot, the L-interval is [lofst_lb, lofst_le].
        const offset_type lofst_lb = 0;
        const offset_type lofst_le = get_lexp(lnode) - 1;

        // Consider 'lnode' as the pivot, the F-interval is [lofst_fb, lofst_fe].
        auto [fnode, lofst_fe] = get_first_overlapped_fnode(lnode);
        offset_type lofst_fb   = lofst_fe - get_fexp(fnode) + 1;

        if ((lofst_fb <= lofst_lb) and (lofst_lb <= lofst_fe) and (lofst_fe <= lofst_le)) {
            DEBUG_ABORT_IF(fnode->get_data().get_tnode());
            lnode->get_data().set_hlink(fnode, lofst_fe - lofst_lb);
            fnode->get_data().set_tnode(lnode);
        }

        reset_lweight(lnode);
    }

    void reset_lnode_links(lbuffer_type& lnodes) {
        for (lnode_type* lnode : lnodes) {
            reset_lnode_link(lnode);
        }
    }

    void reset_fnode_link(fnode_type* fnode) {
        if (fnode->get_data().get_tnode()) {
            reset_fweight(fnode);
            return;
        }

        // Consider 'fnode' as the pivot, the F-interval is [fofst_fb, fofst_fe].
        const offset_type fofst_fb = 0;
        const offset_type fofst_fe = get_fexp(fnode) - 1;

        // Consider 'fnode' as the pivot, the L-interval is [fofst_lb, fofst_le].
        auto [lnode, fofst_le] = get_first_overlapped_lnode(fnode);
        offset_type fofst_lb   = fofst_le - get_lexp(lnode) + 1;

        while (fofst_lb <= fofst_fe) {
            if ((fofst_fb <= fofst_lb) and (fofst_lb <= fofst_fe) and (fofst_fe <= fofst_le)) {
                DEBUG_ABORT_IF(lnode->get_data().get_hnode());
                lnode->get_data().set_hlink(fnode, fofst_fe - fofst_lb);
                fnode->get_data().set_tnode(lnode);
                break;
            }
            lnode    = lnode->get_next();
            fofst_lb = fofst_le + 1;
            fofst_le = fofst_le + get_lexp(lnode);
        }

        reset_fweight(fnode);
    }

    void reset_fnode_links(fbuffer_type& fnodes) {
        for (fnode_type* fnode : fnodes) {
            reset_fnode_link(fnode);
        }
    }

    void reset_lweight(lnode_type* lnode) {
        const lcover_type lcover = get_lcover(lnode);
        lnode->get_data().set_weight(lcover.weight);
        push_fat_lnode(lnode);
    }

    void reset_fweight(fnode_type* fnode) {
        const fcover_type fcover = get_fcover(fnode);
        fnode->get_data().set_weight(fcover.weight);
        push_fat_fnode(fnode);
    }

#ifdef ENABLE_DEBUG_PRINT
    void update_pcmap(const void* p) {
        if (m_pcmap.find(intptr_t(p)) == m_pcmap.end()) {
            m_pcmap.insert(std::make_pair(intptr_t(p), m_pcmax++));
        }
    }
    std::string get_pc(const void* p) {
        update_pcmap(p);
        return tfm::format("[%c]", m_pcmap.find(intptr_t(p))->second);
    }
    std::string get_lnode_str(const lnode_type* p) {
        return tfm::format("%s => %s, F=%s, H=%s, O=%d, W=%d",  //
                           get_pc(p),  //
                           make_run(p->get_data().get_chr(), get_lexp(p)),  //
                           get_pc(p->get_data().get_fnode()),  //
                           get_pc(p->get_data().get_hnode()),  //
                           p->get_data().get_hofst(),  //
                           p->get_data().get_weight());
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
        return tfm::format("L=%s, F=%s, O=%d, N=%d", get_pc(c.lnode), get_pc(c.fnode), c.lofst, c.weight);
    }
    std::string get_fcover_str(const fcover_type& c) {
        return tfm::format("F=%s, L=%s, O=%d, N=%d", get_pc(c.fnode), get_pc(c.lnode), c.fofst, c.weight);
    }
#endif
};

}  // namespace rcomp
