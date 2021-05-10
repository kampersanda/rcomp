/**
 * @file Compressor_GroupLFIG.hpp
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
 * A class for online RLBWT compressor with the grouped LF-interval graph.
 *
 * @tparam t_LIndex Index for grouped L-nodes
 * @tparam t_FIndex Index for grouped F-nodes
 * @tparam t_DivBound The upper bound of dividing overlapped nodes
 */
template <class t_LIndex, class t_FIndex, size_type t_DivBound = 7>
class Compressor_GroupLFIG {
  public:
    using this_type = Compressor_GroupLFIG<t_LIndex, t_FIndex, t_DivBound>;

    using lindex_type = t_LIndex;
    using findex_type = t_FIndex;

    using lnode_type = typename lindex_type::lnode_type;
    using fnode_type = typename findex_type::fnode_type;

    using lcursor_type = typename lnode_type::data_type::lcursor_type;
    using fcursor_type = typename fnode_type::data_type::fcursor_type;

    // 3 is enough, but 4 is used for faster modulo
    using lqueue_type = StaticQueue<lnode_type*, 4>;
    using fqueue_type = StaticQueue<fnode_type*, 4>;

    using lbuffer_type = StaticVector<lnode_type*, t_DivBound + 2>;
    using fbuffer_type = StaticVector<fnode_type*, t_DivBound + 2>;

    static constexpr size_type GROUP_BOUND = lnode_type::data_type::GROUP_BOUND;
    static constexpr size_type DIV_BOUND   = t_DivBound;

    static_assert(1 < GROUP_BOUND);
    static_assert(7 <= DIV_BOUND, "The division bound must be no less than 7");

  private:
    enum class insert_modes : uint8_t {
        // Case C
        MERGE_CURR,
        MERGE_PREV_IN,
        MERGE_PREV_OUT,
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
    size_type   m_em_lrnid = 0;
    offset_type m_em_lofst = -1;

    // The corresponding $'s position in F
    fnode_type* m_em_fnode = nullptr;
    size_type   m_em_frnid = 0;
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
    double m_ns_phase_X_in_B = 0.0;

    size_type m_num_merge_curr     = 0;
    size_type m_num_merge_prev_in  = 0;
    size_type m_num_merge_prev_out = 0;
#endif

  public:
    //! Default constructor
    Compressor_GroupLFIG() = default;

    //! Default destructor
    virtual ~Compressor_GroupLFIG() = default;

    //! Copy constructor (deleted)
    Compressor_GroupLFIG(const Compressor_GroupLFIG&) = delete;

    //! Copy constructor (deleted)
    Compressor_GroupLFIG& operator=(const Compressor_GroupLFIG&) = delete;

    //! Move constructor
    Compressor_GroupLFIG(Compressor_GroupLFIG&&) noexcept = default;

    //! Move constructor
    Compressor_GroupLFIG& operator=(Compressor_GroupLFIG&&) noexcept = default;

    //! Get if the data structure is empty.
    inline bool is_empty() const {
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

        auto [lnode, lcrsr] = get_head_lcursor();

        while (!is_dmmy_lnode(lnode)) {
            if (!is_em_lrun(lnode, lcrsr)) {
                update_run(make_run(lnode->get_data().get_chr(lcrsr), get_lexp(lnode, lcrsr)));
            } else {
                // $-node
                const offset_type this_exp = get_lexp(lnode, lcrsr) - 1;  // without $

                if (m_em_lofst == 0) {
                    // Push front
                    update_run(make_run(END_MARKER, 1));
                    update_run(make_run(lnode->get_data().get_chr(lcrsr), this_exp));
                } else if (m_em_lofst == this_exp) {
                    // Push back
                    update_run(make_run(lnode->get_data().get_chr(lcrsr), this_exp));
                    update_run(make_run(END_MARKER, 1));
                } else {
                    const offset_type new_exp1 = m_em_lofst;
                    const offset_type new_exp2 = this_exp - m_em_lofst;
                    update_run(make_run(lnode->get_data().get_chr(lcrsr), new_exp1));
                    update_run(make_run(END_MARKER, 1));
                    update_run(make_run(lnode->get_data().get_chr(lcrsr), new_exp2));
                }
            }
            std::tie(lnode, lcrsr) = get_next_lcursor(lnode, lcrsr);
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
        size_type         lrnid = m_em_lrnid;
        offset_type       lofst = m_em_lofst;

        while (true) {
            const auto [fnode, frnid, fofst] = get_fpartner(lnode, lrnid, lofst);

            DEBUG_PRINT(tfm::printfln("lnode=%s, lrnid=%d, lofst=%d", get_pc(lnode), lrnid, lofst);)
            DEBUG_PRINT(tfm::printfln("fnode=%s, frnid=%d, fofst=%d", get_pc(fnode), frnid, fofst);)

            if (is_head_fnode(fnode)) {
                break;
            }
            fn(fnode->get_data().get_chr());

            std::tie(lnode, lrnid) = fnode->get_data().get_lrptr(frnid);

            lofst = fofst;
            if (is_em_lrun(lnode, lrnid) and m_em_lofst <= lofst) {
                lofst += 1;
            }
        }
    }

    //! Get the number of characters in the text.
    inline size_type get_num_chars() const {
        return m_num_chars;
    }

    //! Get the allocated memory in bytes.
    size_type get_memory_in_bytes() const {
        size_type mem = 0;
        mem += m_lindex.get_memory_in_bytes();
        mem += m_findex.get_memory_in_bytes();
        mem += sizeof(m_num_chars);
        mem += sizeof(m_em_lnode);
        mem += sizeof(m_em_lrnid);
        mem += sizeof(m_em_lofst);
        mem += sizeof(m_em_fnode);
        mem += sizeof(m_em_frnid);
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
        tfm::reportfln("sec_phase_X:\t%g", utils::ns_to_sec(m_ns_phase_X_in_B));

        tfm::printfln("[Monitor_MergeInCaseC]");
        tfm::reportfln("num_merge_curr:\t%d", m_num_merge_curr);
        tfm::reportfln("num_merge_prev_in:\t%d", m_num_merge_prev_in);
        tfm::reportfln("num_merge_prev_out:\t%d", m_num_merge_prev_out);
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
        auto [curr_fnode, curr_fcrsr] = get_head_fcursor();

        ABORT_IF(!is_dmmy_lnode(curr_fnode->get_data().get_lrptr(curr_fcrsr).first));
        ABORT_IF_NE(curr_fcrsr.get_frnid(), curr_fnode->get_data().get_last_fcursor().get_frnid());
        ABORT_IF_NE(END_MARKER, curr_fnode->get_data().get_chr());

        std::tie(curr_fnode, curr_fcrsr) = get_next_fcursor(curr_fnode, curr_fcrsr);

        while (true) {
            // Check corresponding nodes
            {
                const auto [corr_lnode, corr_lcrsr] = get_corr_lcursor(curr_fnode, curr_fcrsr);
                const auto [corr_fnode, corr_fcrsr] = get_corr_fcursor(corr_lnode, corr_lcrsr);

                ABORT_IF(curr_fnode != corr_fnode);
                ABORT_IF_NE(curr_fcrsr.get_frnid(), corr_fcrsr.get_frnid())
                ABORT_IF_NE(curr_fcrsr.get_frpos(), corr_fcrsr.get_frpos())
            }

            const auto [next_fnode, next_fcrsr] = get_next_fcursor(curr_fnode, curr_fcrsr);
            if (is_head_fnode(next_fnode)) {
                break;
            }

            const uchar_type curr_chr = curr_fnode->get_data().get_chr(curr_fcrsr);
            const uchar_type next_chr = next_fnode->get_data().get_chr(next_fcrsr);

            // F-node characters are not sorted in lex.
            ABORT_IF_LT(next_chr, curr_chr);

            if (curr_chr < next_chr) {
                curr_fnode = next_fnode;
                curr_fcrsr = next_fcrsr;
                continue;
            }

            const auto [curr_lnode, curr_lcrsr] = get_corr_lcursor(curr_fnode, curr_fcrsr);
            const auto [next_lnode, next_lcrsr] = get_corr_lcursor(next_fnode, next_fcrsr);

            const loint_type curr_order = lo_common::get_basic_order(curr_lnode);
            const loint_type next_order = lo_common::get_basic_order(next_lnode);

            if (curr_order != next_order) {
                ABORT_IF_LE(next_order, curr_order);
            } else {
                ABORT_IF_LE(next_lcrsr.get_lrpos(), curr_lcrsr.get_lrpos());
            }

            curr_fnode = next_fnode;
            curr_fcrsr = next_fcrsr;
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
            fpos_b = fpos_b + get_sum_fexps(fnode);
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
#endif

  private:
    //! Get the head F-node.
    inline bool is_head_fnode(const fnode_type* fnode) const {
        return fnode == m_findex.get_head();
    }
    //! Get the dummy L-node.
    inline bool is_dmmy_lnode(const lnode_type* lnode) const {
        return lnode == m_lindex.get_head();
    }
    // Does L-run (lnode, lrnid) have $ mark?
    inline bool is_em_lrun(const lnode_type* lnode, const size_type lrnid) const {
        return (lnode == m_em_lnode) and (lrnid == m_em_lrnid);
    }
    // Does L-run (lnode, lcrsr) have $ mark?
    inline bool is_em_lrun(const lnode_type* lnode, const lcursor_type lcrsr) const {
        return (lnode == m_em_lnode) and (lcrsr.get_lrnid() == m_em_lrnid);
    }

    inline bool is_em_frun(const fnode_type* fnode, const size_type frnid) const {
        return fnode == m_em_fnode and frnid == m_em_frnid;
    }
    inline bool is_em_frun(const fnode_type* fnode, const fcursor_type fcrsr) const {
        return fnode == m_em_fnode and fcrsr.get_frnid() == m_em_frnid;
    }

    inline fnode_type* get_head_fnode() const {
        return m_findex.get_head();
    }
    inline lnode_type* get_head_lnode() const {
        return m_lindex.get_head()->get_next();
    }
    inline fnode_type* get_tail_fnode() const {
        return m_findex.get_head()->get_prev();
    }
    inline lnode_type* get_tail_lnode() const {
        return m_lindex.get_head()->get_prev();
    }

    /* * * * * * * * * * * * * * * * * * * * *
     *
     *     Iterator tools for F/L-runs
     *
     * * * * * * * * * * * * * * * * * * * * */

    inline std::pair<fnode_type*, fcursor_type> get_head_fcursor() {
        fnode_type* fnode = get_head_fnode();
        return {fnode, fnode->get_data().get_first_fcursor()};
    }
    inline std::pair<const fnode_type*, fcursor_type> get_head_fcursor() const {
        const fnode_type* fnode = get_head_fnode();
        return {fnode, fnode->get_data().get_first_fcursor()};
    }
    inline std::pair<fnode_type*, fcursor_type> get_tail_fcursor() {
        fnode_type* fnode = get_tail_fnode();
        return {fnode, fnode->get_data().get_last_fcursor()};
    }
    inline std::pair<const fnode_type*, fcursor_type> get_tail_fcursor() const {
        const fnode_type* fnode = get_tail_fnode();
        return {fnode, fnode->get_data().get_last_fcursor()};
    }
    inline std::pair<lnode_type*, lcursor_type> get_head_lcursor() {
        lnode_type* lnode = get_head_lnode();
        return {lnode, lnode->get_data().get_first_lcursor()};
    }
    inline std::pair<const lnode_type*, lcursor_type> get_head_lcursor() const {
        const lnode_type* lnode = get_head_lnode();
        return {lnode, lnode->get_data().get_first_lcursor()};
    }
    inline std::pair<lnode_type*, lcursor_type> get_tail_lcursor() {
        lnode_type* lnode = get_tail_lnode();
        return {lnode, lnode->get_data().get_last_lcursor()};
    }
    inline std::pair<const lnode_type*, lcursor_type> get_tail_lcursor() const {
        const lnode_type* lnode = get_tail_lnode();
        return {lnode, lnode->get_data().get_last_lcursor()};
    }

    inline std::pair<fnode_type*, fcursor_type> get_next_fcursor(fnode_type* fnode, fcursor_type fcrsr) {
        if (fnode->get_data().is_last_fcursor(fcrsr)) {
            fnode = fnode->get_next();
            fcrsr = fnode->get_data().get_first_fcursor();
        } else {
            fcrsr = fnode->get_data().get_next_fcursor(fcrsr);
        }
        return {fnode, fcrsr};
    }
    inline std::pair<const fnode_type*, fcursor_type> get_next_fcursor(const fnode_type* fnode,
                                                                       fcursor_type      fcrsr) const {
        if (fnode->get_data().is_last_fcursor(fcrsr)) {
            fnode = fnode->get_next();
            fcrsr = fnode->get_data().get_first_fcursor();
        } else {
            fcrsr = fnode->get_data().get_next_fcursor(fcrsr);
        }
        return {fnode, fcrsr};
    }
    inline std::pair<lnode_type*, lcursor_type> get_next_lcursor(lnode_type* lnode, lcursor_type lcrsr) {
        if (lnode->get_data().is_last_lcursor(lcrsr)) {
            lnode = lnode->get_next();
            lcrsr = lnode->get_data().get_first_lcursor();
        } else {
            lcrsr = lnode->get_data().get_next_lcursor(lcrsr);
        }
        return {lnode, lcrsr};
    }
    inline std::pair<const lnode_type*, lcursor_type> get_next_lcursor(const lnode_type* lnode,
                                                                       lcursor_type      lcrsr) const {
        if (lnode->get_data().is_last_lcursor(lcrsr)) {
            lnode = lnode->get_next();
            lcrsr = lnode->get_data().get_first_lcursor();
        } else {
            lcrsr = lnode->get_data().get_next_lcursor(lcrsr);
        }
        return {lnode, lcrsr};
    }

    inline bool is_first_lofst(const lnode_type* lnode, size_type lrnid, offset_type lofst) const {
        return is_first_lofst(lnode, lnode->get_data().get_lcursor(lrnid), lofst);
    }
    inline bool is_first_lofst(const lnode_type* lnode, lcursor_type lcrsr, offset_type lofst) const {
        return lofst == 0 and lnode->get_data().is_first_lcursor(lcrsr);
    }
    inline bool is_last_lofst(const lnode_type* lnode, size_type lrnid, offset_type lofst) const {
        return is_last_lofst(lnode, lnode->get_data().get_lcursor(lrnid), lofst);
    }
    inline bool is_last_lofst(const lnode_type* lnode, lcursor_type lcrsr, offset_type lofst) const {
        const auto& ldata = lnode->get_data();
        return ldata.is_last_lcursor(lcrsr) and ldata.get_exp(lcrsr) == size_type(lofst + 1);
    }
    inline bool is_first_fofst(const fnode_type* fnode, size_type frnid, offset_type fofst) const {
        return is_first_fofst(fnode, fnode->get_data().get_fcursor(frnid), fofst);
    }
    inline bool is_first_fofst(const fnode_type* fnode, fcursor_type fcrsr, offset_type fofst) const {
        return fofst == 0 and fnode->get_data().is_first_fcursor(fcrsr);
    }
    inline bool is_last_fofst(const fnode_type* fnode, size_type frnid, offset_type fofst) const {
        return is_last_fofst(fnode, fnode->get_data().get_fcursor(frnid), fofst);
    }
    inline bool is_last_fofst(const fnode_type* fnode, fcursor_type fcrsr, offset_type fofst) const {
        const auto& fdata = fnode->get_data();
        return fdata.is_last_fcursor(fcrsr) and fdata.get_exp(fcrsr) == size_type(fofst + 1);
    }

    inline std::pair<lnode_type*, lcursor_type> get_corr_lcursor(const fnode_type* fnode, fcursor_type fcrsr) const {
        const auto [lnode, lrnid] = fnode->get_data().get_lrptr(fcrsr);
        return {lnode, lnode->get_data().get_lcursor(lrnid)};
    }
    inline std::pair<fnode_type*, fcursor_type> get_corr_fcursor(const lnode_type* lnode, lcursor_type lcrsr) const {
        const auto [fnode, frnid] = lnode->get_data().get_frptr(lcrsr);
        return {fnode, fnode->get_data().get_fcursor(frnid)};
    }

    // Get the exponent of the run-unit with ID 'frnid' in F-group 'fnode'.
    inline offset_type get_fexp(const fnode_type* fnode, const fcursor_type& fcrsr) const {
        return fnode->get_data().get_exp(fcrsr);
    }
    // Get the exponent of the run-unit with ID 'lrnid' in L-group 'lnode'.
    // Note if the run has $ mark, the retruned value is icremented by 1.
    inline offset_type get_lexp(const lnode_type* lnode, const size_type lrnid) const {
        const offset_type lexp = lnode->get_data().get_exp(lrnid);
        return is_em_lrun(lnode, lrnid) ? lexp + 1 : lexp;
    }
    inline offset_type get_lexp(const lnode_type* lnode, const lcursor_type& lcrsr) const {
        const offset_type lexp = lnode->get_data().get_exp(lcrsr);
        return is_em_lrun(lnode, lcrsr) ? lexp + 1 : lexp;
    }
    // Get the sum of exponents of runs in F-group 'fnode'.
    inline offset_type get_sum_fexps(const fnode_type* fnode) const {
        return fnode->get_data().get_sum_exps();
    }
    // Get the sum of exponents of runs in L-group 'lnode'.
    // Note if the group has $ mark, the retruned value is icremented by 1.
    inline offset_type get_sum_lexps(const lnode_type* lnode) const {
        const offset_type lsize = lnode->get_data().get_sum_exps();
        return lnode == m_em_lnode ? lsize + 1 : lsize;
    }

    // Check whether 'ins_mode' is MERGE_*
    inline bool is_merge_mode(insert_modes ins_mode) {
        return insert_modes::MERGE_CURR <= ins_mode and ins_mode <= insert_modes::MERGE_PREV_OUT;
    }

    // The distance between the begin of the given F-node and that of the linked L-node.
    inline offset_type get_begin_offset(const fnode_type* fnode) const {
        DEBUG_ABORT_IF(!fnode->get_data().get_tnode());
        DEBUG_ABORT_IF_LT(0, fnode->get_data().get_tofst() - get_sum_fexps(fnode) + 1);
        return -1 * (fnode->get_data().get_tofst() - get_sum_fexps(fnode) + 1);
    }
    inline offset_type get_end_offset(const fnode_type* fnode) const {
        DEBUG_ABORT_IF(!fnode->get_data().get_tnode());
        const lnode_type* lnode = fnode->get_data().get_tnode();
        DEBUG_ABORT_IF_LT(0, lnode->get_data().get_hofst() - get_sum_lexps(lnode) + 1);
        return -1 * (lnode->get_data().get_hofst() - get_sum_lexps(lnode) + 1);
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
    inline insert_modes check_joinable(const uchar_type new_chr) const {
        DEBUG_ABORT_IF(is_dmmy_lnode(m_em_lnode));
        DEBUG_ABORT_IF_OUT(m_em_lofst, 0, get_lexp(m_em_lnode, m_em_lrnid) - 1);  // -1 is for $

        const lcursor_type em_lcrsr = m_em_lnode->get_data().get_lcursor(m_em_lrnid);

        // Mergeable with the current node?
        if (new_chr == m_em_lnode->get_data().get_chr(em_lcrsr)) {
            return insert_modes::MERGE_CURR;  // Case C
        }

        if (m_em_lofst == 0) {
            // Try to check the previous run
            if (!m_em_lnode->get_data().is_first_lcursor(em_lcrsr)) {
                // Mergeable with the previous run in the same group?
                const lcursor_type prev_lcurs = m_em_lnode->get_data().get_prev_lcursor(em_lcrsr);
                if (new_chr == m_em_lnode->get_data().get_chr(prev_lcurs)) {
                    return insert_modes::MERGE_PREV_IN;  // Case C
                }
            } else {
                // Mergeable with the previous run in the previous group?
                const lnode_type*  prev_lnode = m_em_lnode->get_prev();
                const lcursor_type prev_lcurs = prev_lnode->get_data().get_last_lcursor();
                if (new_chr == prev_lnode->get_data().get_chr(prev_lcurs)) {
                    return insert_modes::MERGE_PREV_OUT;  // Case C
                }
            }
        }

        // new_chr will not be merged with the next run because $-mark is never placed at the end of some L-run.
        // Only when $-mark is placed at the last of L, it will be placed at the end of the last L-run.
        // But, then new_chr is not be merged with the next run anyway.

        if (m_em_lofst == 0) {
            return insert_modes::SPLIT_FRONT;  // Case B
        } else if (m_em_lofst == get_lexp(m_em_lnode, em_lcrsr) - 1) {  // -1 is for $
            return insert_modes::SPLIT_BACK;  // Case B
        } else {
            return insert_modes::SPLIT_MIDDLE;  // Case B
        }
    }

    //! Make the data structure for the RLBWT text of T = (new_chr,$).
    void extend_init(const uchar_type new_chr) {
        // Dummy node for the partner of $'s F-node (for convenience)
        auto [dlr_lnode, dlr_lcrsr] = m_lindex.clear();
        // $'s node at the top of F
        auto dlr_fnode = m_findex.clear(dlr_lnode, dlr_lcrsr.get_lrnid()).first;

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
    }

    //! Extend the RLBWT text in Case C.
    void extend_C(const uchar_type, const insert_modes ins_mode) {
        DEBUG_ABORT_IF(!is_merge_mode(ins_mode));
        DEBUG_PRINT(tfm::printfln("\n==== Case C ====");)

        DEBUG_ABORT_IF(!m_fat_lnodes.is_empty());
        DEBUG_ABORT_IF(!m_fat_fnodes.is_empty());

        // Variables for the L-intervals to be inserted.
        // i.e., the L-interval is of size 'ins_lsize' starting at the 'ins_lofst'-th position in group 'ins_lnode'.
        // 'ins_lrnid' indicates the L-run to be incremented.
        lnode_type*  ins_lnode = nullptr;
        lcursor_type ins_lcrsr = {};

        // Variables for the F-intervals to be inserted.
        // i.e., the F-interval is of size 'ins_fsize' starting at the 'ins_fofst'-th position in group 'ins_fnode'.
        // 'ins_frnid' indicates the F-run to be incremented.
        fnode_type*  ins_fnode = nullptr;
        fcursor_type ins_fcrsr = {};
        offset_type  ins_fofst = 0;

        // [NOTE]
        // Depending on 'ins_mode', the L-node with $ can be updated.
        // But Do NOT change the value of 'm_em_lnode' until Phase E.
        // This is because this change will affect 'get_sum_lexps' (in Phase C).

        //
        // Phase B
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase B ==");)
        STAT_MONITOR(m_timer_4_phase.start();)

        const lcursor_type em_lcrsr = m_em_lnode->get_data().get_lcursor(m_em_lrnid);

        if (ins_mode == insert_modes::MERGE_CURR) {
            // Then, the new char is inserted into 'm_em_lnode'
            DEBUG_ABORT_IF_OUT(m_em_lofst, 0, get_lexp(m_em_lnode, em_lcrsr) - 1);
            // L (the overlapped F-nodes do not need to be updated)
            ins_lnode = m_em_lnode;
            ins_lcrsr = em_lcrsr;
            // F
            std::tie(ins_fnode, ins_fcrsr) = get_corr_fcursor(ins_lnode, ins_lcrsr);
            ins_fofst                      = m_em_lofst;
            STAT_MONITOR(m_num_merge_curr += 1;)
        } else if (ins_mode == insert_modes::MERGE_PREV_IN) {
            // Then, the new char is inserted into the previous run in the same group.
            DEBUG_ABORT_IF_NE(m_em_lofst, 0);
            DEBUG_ABORT_IF(m_em_lnode->get_data().is_first_lcursor(em_lcrsr));
            // L
            ins_lnode = m_em_lnode;
            ins_lcrsr = m_em_lnode->get_data().get_prev_lcursor(em_lcrsr);
            // F
            std::tie(ins_fnode, ins_fcrsr) = get_corr_fcursor(ins_lnode, ins_lcrsr);
            ins_fofst                      = get_fexp(ins_fnode, ins_fcrsr);
            STAT_MONITOR(m_num_merge_prev_in += 1;)
        } else if (ins_mode == insert_modes::MERGE_PREV_OUT) {
            // Then, the new char will be inserted into the previous run in the previous group.
            DEBUG_ABORT_IF_NE(m_em_lofst, 0);
            DEBUG_ABORT_IF(!m_em_lnode->get_data().is_first_lcursor(em_lcrsr));
            // L
            ins_lnode = m_em_lnode->get_prev();
            ins_lcrsr = ins_lnode->get_data().get_last_lcursor();
            // F
            std::tie(ins_fnode, ins_fcrsr) = get_corr_fcursor(ins_lnode, ins_lcrsr);
            ins_fofst                      = get_fexp(ins_fnode, ins_fcrsr);
            STAT_MONITOR(m_num_merge_prev_out += 1;)
        } else {
            ABORT_IF(true);
        }

        STAT_MONITOR(m_ns_phase_B_in_C += m_timer_4_phase.stop_and_get_ns();)
        DEBUG_PRINT(tfm::printfln("ins_lnode=%s", get_pc(ins_lnode));)
        DEBUG_PRINT(tfm::printfln("ins_fnode=%s, ins_fofst=%d", get_pc(ins_fnode), ins_fofst);)
        DEBUG_PRINT(tfm::printfln("ins_lcrsr=%s, ins_fcrsr=%s", get_lcrsr_str(ins_lcrsr), get_fcrsr_str(ins_fcrsr));)

        //
        // Phase C (Remove H/T-links of overlapped nodes)
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase C ==");)
        STAT_MONITOR(m_timer_4_phase.start();)

        // $-replacement
        if (ins_mode == insert_modes::MERGE_PREV_OUT) {
            lnode_type* em_tnode = m_em_fnode->get_data().get_tnode();

            const bool is_first_em_fofst = is_first_fofst(m_em_fnode, m_em_frnid, m_em_fofst);
            const bool is_last_em_fofst  = is_last_fofst(m_em_fnode, m_em_frnid, m_em_fofst);

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
                    fnode_type*       fnode_befo    = m_em_fnode;
                    lnode_type*       lnode_befo    = m_em_lnode->get_prev();
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
                    fnode_type*       fnode_aftr    = m_em_fnode->get_next();
                    lnode_type*       lnode_aftr    = m_em_lnode;
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
        STAT_MONITOR(m_ns_phase_C_in_C += m_timer_4_phase.stop_and_get_ns();)

        //
        // Phase D
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase D ==");)
        STAT_MONITOR(m_timer_4_phase.start();)

        ins_lnode->get_data().prefetch();
        ins_fnode->get_data().prefetch();

        ins_lnode->get_data().add_exp(ins_lcrsr, 1);
        ins_fnode->get_data().add_exp(ins_fcrsr, 1);

        STAT_MONITOR(m_ns_phase_D_in_C += m_timer_4_phase.stop_and_get_ns();)
        DEBUG_PRINT(tfm::printfln("ins_lnode=%s, ins_lcrsr=%s", get_pc(ins_lnode), get_lcrsr_str(ins_lcrsr));)
        DEBUG_PRINT(tfm::printfln("ins_lrun=%s", get_lrun_str(ins_lnode, ins_lcrsr));)

        //
        // Phase E
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase E ==");)
        STAT_MONITOR(m_timer_4_phase.start();)

        // Clear since the L-run has already updated.
        m_em_lnode = nullptr;
        m_em_lrnid = 0;
        m_em_lofst = 0;

        // Cache to the next process
        m_em_fnode = ins_fnode;
        m_em_frnid = ins_fcrsr.get_frnid();
        m_em_fofst = ins_fofst;

        // $-query
        std::tie(m_em_lnode, m_em_lrnid, m_em_lofst) = get_lpartner(ins_fnode, ins_fcrsr, ins_fofst);

        STAT_MONITOR(m_ns_phase_E_in_C += m_timer_4_phase.stop_and_get_ns();)
        DEBUG_PRINT(tfm::printfln("m_em_lnode=%s, m_em_lrnid=%d, m_em_lofst=%d",  //
                                  get_pc(m_em_lnode), m_em_lrnid, m_em_lofst);)

        //
        // Phase F
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase F ==");)
        STAT_MONITOR(m_timer_4_phase.start();)

        // ?-replacement
        {
            const bool is_last_em_fofst  = is_last_fofst(m_em_fnode, ins_fcrsr, m_em_fofst);
            const bool is_first_em_lofst = is_first_lofst(m_em_lnode, m_em_lrnid, m_em_lofst);

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
        DEBUG_ABORT_IF((ins_mode == insert_modes::SPLIT_FRONT) and (m_em_lofst != 0));
        DEBUG_ABORT_IF((ins_mode == insert_modes::SPLIT_BACK) and (m_em_lofst != get_lexp(m_em_lnode, m_em_lrnid) - 1));
        DEBUG_ABORT_IF((ins_mode == insert_modes::SPLIT_MIDDLE) and  //
                       (m_em_lofst < 1) and (get_lexp(m_em_lnode, m_em_lrnid) - 2 < m_em_lofst));

        DEBUG_ABORT_IF(!m_fat_lnodes.is_empty());
        DEBUG_ABORT_IF(!m_fat_fnodes.is_empty());

        DEBUG_PRINT(tfm::printfln("\n==== Case B ====");)

        // Note:
        // In Case B, L-interval is never changed because the split will be done within the group.

        //
        // Phase B
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase B ==");)
        STAT_MONITOR(m_timer_4_phase.start();)

        // Variables for the predecessor F-run
        fnode_type*  pred_fnode = nullptr;
        fcursor_type pred_fcrsr = {};
        {
            lnode_type*  pred_lnode = m_em_lnode;
            lcursor_type pred_lcrsr = m_em_lnode->get_data().get_lcursor(m_em_lrnid);

            if (ins_mode == insert_modes::SPLIT_FRONT) {
                if (!pred_lnode->get_data().is_first_lcursor(pred_lcrsr)) {
                    // The predecessor L-node is in the group.
                    pred_lcrsr = pred_lnode->get_data().get_prev_lcursor(pred_lcrsr);
                } else {
                    // The predecessor L-node is the tail of the previous group.
                    pred_lnode = pred_lnode->get_prev();
                    pred_lcrsr = pred_lnode->get_data().get_last_lcursor();
                }
            }

            // Predecessor query in O(log(r/B) + B) time
            std::tie(pred_fnode, pred_fcrsr) = predecessor(new_chr, pred_lnode, pred_lcrsr);
        }

        // A new F-run will be inserted before/after F-run (ins_fnode, ins_fcrsr).
        auto [ins_fnode, ins_fcrsr] = get_next_fcursor(pred_fnode, pred_fcrsr);

        // In Case B, the three patterns should be considered.
        //  - INSERT_BEFORE: A new F-run is inserted before F-run (ins_fnode, ins_fcrsr) within the group.
        //  - INSERT_AFTER : A new F-run is inserted after F-run (ins_fnode, ins_fcrsr) within the group.
        //  - NEW_CREATE   : A new F-run is inserted before F-group ins_fnode as a new group.
        enum class fnode_insert_modes : uint8_t { INSERT_BEFORE, INSERT_AFTER, NEW_CREATE };
        fnode_insert_modes fins_mode = fnode_insert_modes::INSERT_BEFORE;

        // Set the pivotal F-run and the insertion mode.
        if (new_chr == ins_fnode->get_data().get_chr()) {
            fins_mode = fnode_insert_modes::INSERT_BEFORE;
        } else if (new_chr == pred_fnode->get_data().get_chr()) {
            DEBUG_ABORT_IF(!ins_fnode->get_data().is_first_fcursor(ins_fcrsr));
            ins_fnode = pred_fnode;
            ins_fcrsr = pred_fcrsr;
            fins_mode = fnode_insert_modes::INSERT_AFTER;
        } else {
            DEBUG_ABORT_IF(!ins_fnode->get_data().is_first_fcursor(ins_fcrsr));
            fins_mode = fnode_insert_modes::NEW_CREATE;
        }

        // Variables for the F-intervals to be inserted.
        // i.e., the F-interval is of size 'ins_fsize' starting at the 'ins_fofst'-th position in group 'ins_fnode'.
        // offset_type ins_fofst = 0;
        // offset_type ins_fsize = get_sum_fexps(ins_fnode) + 1;  // +1 for new char

        STAT_MONITOR(m_ns_phase_B_in_B += m_timer_4_phase.stop_and_get_ns();)
        DEBUG_PRINT(tfm::printfln("ins_fnode=%s, ins_fcrsr=%s", get_pc(ins_fnode), get_fcrsr_str(ins_fcrsr));)
        DEBUG_PRINT(tfm::printfln("pred_fnode=%s, pred_fcrsr=%s", get_pc(pred_fnode), get_fcrsr_str(pred_fcrsr));)

        //
        // Phase C
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase C ==");)

        //
        // Phase D
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase D ==");)
        STAT_MONITOR(m_timer_4_phase.start();)

        // Variables for indicating a new L-run inserted
        lnode_type*  ins_lnode = m_em_lnode;
        lcursor_type ins_lcrsr = m_em_lnode->get_data().get_lcursor(m_em_lrnid);

        fnode_type*  div_fnode = nullptr;
        fcursor_type div_fcrsr = {};

        // Insert a new L-run
        if (ins_mode == insert_modes::SPLIT_FRONT) {
            ins_lcrsr = ins_lnode->get_data().insert_before(ins_lcrsr, new_chr, 1);
        } else if (ins_mode == insert_modes::SPLIT_BACK) {
            ins_lcrsr = ins_lnode->get_data().insert_after(ins_lcrsr, new_chr, 1);
        } else if (ins_mode == insert_modes::SPLIT_MIDDLE) {
            // Get the corresponding F-run to be split.
            std::tie(div_fnode, div_fcrsr) = get_corr_fcursor(ins_lnode, ins_lcrsr);
            // Split and Link
            ins_lcrsr = ins_lnode->get_data().split_after(ins_lcrsr, m_em_lofst);
            div_fcrsr = div_fnode->get_data().split_after(div_fcrsr, ins_lnode, ins_lcrsr.get_lrnid());
            ins_lnode->get_data().set_frptr(ins_lcrsr, div_fnode, div_fcrsr.get_frnid());
            // Insert a new L-run between the split L-runs
            ins_lcrsr = ins_lnode->get_data().insert_before(ins_lcrsr, new_chr, 1);
        } else {
            ABORT_IF(true);
        }

        // Insert a new F-run
        if (fins_mode == fnode_insert_modes::INSERT_BEFORE) {
            ins_fcrsr = ins_fnode->get_data().insert_before(ins_fcrsr, ins_lnode, ins_lcrsr.get_lrnid());
            ins_lnode->get_data().set_frptr(ins_lcrsr, ins_fnode, ins_fcrsr.get_frnid());
        } else if (fins_mode == fnode_insert_modes::INSERT_AFTER) {
            ins_fcrsr = ins_fnode->get_data().insert_after(ins_fcrsr, ins_lnode, ins_lcrsr.get_lrnid());
            ins_lnode->get_data().set_frptr(ins_lcrsr, ins_fnode, ins_fcrsr.get_frnid());
        } else if (fins_mode == fnode_insert_modes::NEW_CREATE) {
            DEBUG_ABORT_IF(ins_fnode != pred_fnode->get_next());
            std::tie(ins_fnode, ins_fcrsr) = m_findex.insert_after(pred_fnode, ins_lnode, ins_lcrsr.get_lrnid());
        } else {
            ABORT_IF(true);
        }

        STAT_MONITOR(m_ns_phase_D_in_B += m_timer_4_phase.stop_and_get_ns();)
        DEBUG_PRINT(tfm::printfln("ins_lnode=%s, ins_lcrsr=%s", get_pc(ins_lnode), get_lcrsr_str(ins_lcrsr));)
        DEBUG_PRINT(tfm::printfln("ins_lrun=%s", get_lrun_str(ins_lnode, ins_lcrsr));)
        DEBUG_PRINT(tfm::printfln("div_fnode=%s, div_fcrsr=%s", get_pc(div_fnode), get_fcrsr_str(div_fcrsr));)
        DEBUG_PRINT(tfm::printfln("div_frun=%s", get_frun_str(div_fnode, div_fcrsr));)

        //
        // Phase E
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase E ==");)
        STAT_MONITOR(m_timer_4_phase.start();)

        // Clear since the L-run has already updated.
        m_em_lnode = nullptr;
        m_em_lrnid = 0;
        m_em_lofst = 0;

        // Cache to the next process
        m_em_fnode = ins_fnode;
        m_em_frnid = ins_fcrsr.get_frnid();
        m_em_fofst = 0;  // always 0 because a new F-run is inserted.

        // $-query
        std::tie(m_em_lnode, m_em_lrnid, m_em_lofst) = get_lpartner(ins_fnode, ins_fcrsr, 0);

        STAT_MONITOR(m_ns_phase_E_in_B += m_timer_4_phase.stop_and_get_ns();)
        DEBUG_PRINT(tfm::printfln("m_em_lnode=%s, m_em_lrnid=%d, m_em_lofst=%d",  //
                                  get_pc(m_em_lnode), m_em_lrnid, m_em_lofst);)

        // Phase F
        DEBUG_PRINT(tfm::printfln("\n== Phase F ==");)
        STAT_MONITOR(m_timer_4_phase.start();)

        // ?-replacement
        {
            const bool is_last_em_fofst  = is_last_fofst(m_em_fnode, ins_fcrsr, m_em_fofst);
            const bool is_first_em_lofst = is_first_lofst(m_em_lnode, m_em_lrnid, m_em_lofst);

            if (fins_mode != fnode_insert_modes::NEW_CREATE) {
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
            } else {  // fins_mode == fnode_insert_modes::NEW_CREATE
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

        STAT_MONITOR(m_ns_phase_F_in_B += m_timer_4_phase.stop_and_get_ns();)

        //
        // Extra Phase: Group split
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase X (Extra Phase) ==");)
        STAT_MONITOR(m_timer_4_phase.start();)

        DEBUG_ABORT_IF(!ins_lnode);
        if (ins_lnode->get_data().get_num_runs() >= GROUP_BOUND) {
            fnode_type* updated_fnode = divide_lnode(ins_lnode);
            push_fat_fnode(updated_fnode);
        }
        DEBUG_ABORT_IF(!ins_fnode);
        if (ins_fnode->get_data().get_num_runs() >= GROUP_BOUND) {
            lnode_type* updated_lnode = divide_fnode(ins_fnode);
            push_fat_lnode(updated_lnode);
        }
        if (div_fnode and div_fnode->get_data().get_num_runs() >= GROUP_BOUND) {
            lnode_type* updated_lnode = divide_fnode(div_fnode);
            push_fat_lnode(updated_lnode);
        }

        STAT_MONITOR(m_ns_phase_X_in_B += m_timer_4_phase.stop_and_get_ns();)
    }

    /**
     * @brief Devide nodes whose weight is no less than DIV_BOUND, stored in m_fat_nodes.
     * @return The number of nodes divided.
     */
    size_type divide_covers() {
        size_type num_divided = 0;

        while (!m_fat_lnodes.is_empty()) {
            lnode_type* lnode = m_fat_lnodes.pop();
            DEBUG_ABORT_IF(!lnode);

            if (lnode->get_data().get_weight() >= DIV_BOUND) {
                const lcover_type lcover = get_lcover(lnode);
                DEBUG_ABORT_IF_NE(lcover.weight, lnode->get_data().get_weight());

                lnode_type* new_updated_lnode = divide_lcover(lcover);
                push_fat_lnode(new_updated_lnode);
                num_divided += 1;
            }
        }

        while (!m_fat_fnodes.is_empty()) {
            fnode_type* fnode = m_fat_fnodes.pop();
            DEBUG_ABORT_IF(!fnode);

            if (fnode->get_data().get_weight() >= DIV_BOUND) {
                const fcover_type fcover = get_fcover(fnode);
                DEBUG_ABORT_IF_NE(fcover.weight, fnode->get_data().get_weight());

                fnode_type* new_updated_fnode = divide_fcover(fcover);
                push_fat_fnode(new_updated_fnode);
                num_divided += 1;
            }
        }

        return num_divided;
    }

    /**
     * @brief Devide the fat L-node whose weight is no less than DIV_BOUND.
     * @param[in] lcover The fat L-node and information of the overlapped F-node.
     * @return The L-node whose weight can be updated.
     * @note To split the F-run corresponding to the split L-run results in change of the weight of some L-node.
     */
    lnode_type* divide_lcover(const lcover_type& lcover) {
        DEBUG_ABORT_IF(!m_lnodes_buffer.is_empty());
        DEBUG_ABORT_IF(!m_fnodes_buffer.is_empty());
        DEBUG_ABORT_IF_LT(lcover.weight, DIV_BOUND);

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

        // If div_lnode == m_em_lnode, may need to update m_em_l*.
        // So, in advance keep the offset of $ position in the L-node.
        offset_type em_lofst_in_node = -1;
        if (div_lnode == m_em_lnode) {
            // Note get_offset does not consider the imaginary $-mark,
            // but the processing does not have a problem.
            em_lofst_in_node = m_em_lnode->get_data().get_offset(m_em_lrnid) + m_em_lofst;
        }

        // Search the L-run to be divided.
        offset_type div_lofst_in_node = lcover.lofst;
        {
            const size_type   fnum_befo = lcover.weight / 2;
            const fnode_type* fnode     = lcover.fnode;

            for (size_type i = 1; i < fnum_befo; i++) {
                fnode = fnode->get_next();
                div_lofst_in_node += get_sum_fexps(fnode);
            }

            if (em_lofst_in_node == div_lofst_in_node) {
                // Then, after split, $-mark will be placed at the end of the front L-run.
                // So, shift the division point.
                fnode = fnode->get_next();
                div_lofst_in_node += get_sum_fexps(fnode);
                ABORT_IF_LE(get_sum_lexps(div_lnode), div_lofst_in_node);
            }
        }
        // Here, div_lnode[div_lofst_in_node] is the last character of the front part.

        DEBUG_PRINT(tfm::printfln("div_lofst_in_node=%d", div_lofst_in_node);)
        DEBUG_PRINT(tfm::printfln("em_lofst_in_node=%d", em_lofst_in_node);)

        //
        // Phase C: Remove head links of divided nodes.
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase C ==");)

        get_overlapped_fnodes(div_lnode, 0, get_sum_lexps(div_lnode), m_fnodes_buffer);
        remove_fnode_links(m_fnodes_buffer);

        //
        // Phase D: Update nodes
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase D ==");)
        ABORT_IF_EQ(em_lofst_in_node, div_lofst_in_node);

        // 1) Split the target L-run in the node.
        lcursor_type div_lcrsr_befo = {};
        offset_type  div_lofst_befo = -1;

        if (div_lnode == m_em_lnode and em_lofst_in_node < div_lofst_in_node) {
            // Then, $-mark will be placed at the before part, and consider the imaginary $-mark in access_with_offset.
            std::tie(div_lcrsr_befo, div_lofst_befo) = div_lnode->get_data().access_with_offset(div_lofst_in_node - 1);
        } else {
            std::tie(div_lcrsr_befo, div_lofst_befo) = div_lnode->get_data().access_with_offset(div_lofst_in_node);
        }
        const size_type div_lexp_befo = div_lofst_befo + 1;

        DEBUG_PRINT(tfm::printfln("div_lcrsr_befo=%s", get_lcrsr_str(div_lcrsr_befo));)
        DEBUG_PRINT(tfm::printfln("div_lofst_befo=%d, div_lexp_befo=%d", div_lofst_befo, div_lexp_befo);)

        // Need to consider OV change depending on L-run split
        fnode_type* div_fnode = nullptr;

        // Do NOT use get_lexp()
        if (div_lnode->get_data().get_exp(div_lcrsr_befo) != div_lexp_befo) {  // the L-run needs to be split?
            std::tie(div_fnode, std::ignore) = split_lrun(div_lnode, div_lcrsr_befo, div_lexp_befo);
        }

        // 2) Divide the L-node
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
        // Release
        m_fnodes_buffer.clear();

        {
            lcover_type lcover = get_lcover(div_lnode);
            div_lnode->get_data().set_weight(lcover.weight);
        }
        {
            lcover_type lcover = get_lcover(new_lnode);
            new_lnode->get_data().set_weight(lcover.weight);
            lcover.fnode->get_data().add_weight(1);  // for new_lnode
        }

        //
        // Extra X: Group split
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase X ==");)

        lnode_type* updated_lnode = nullptr;
        if (div_fnode and div_fnode->get_data().get_num_runs() >= GROUP_BOUND) {
            updated_lnode = divide_fnode(div_fnode);
        }
        return updated_lnode;
    }

    /**
     * @brief Devide the fat F-node whose weight is no less than DIV_BOUND.
     * @param[in] fcover The fat F-node and information of the overlapped L-node.
     * @return The F-node whose weight can be updated.
     * @note To split the L-run corresponding to the split F-run results in change of the weight of some F-node.
     */
    fnode_type* divide_fcover(const fcover_type& fcover) {
        DEBUG_ABORT_IF(!m_lnodes_buffer.is_empty());
        DEBUG_ABORT_IF(!m_fnodes_buffer.is_empty());
        DEBUG_ABORT_IF_LT(fcover.weight, DIV_BOUND);

        if (!fcover.fnode) {
            return nullptr;
        }

        DEBUG_PRINT(tfm::printfln("\n==== Case F ====");)
        DEBUG_PRINT(tfm::printfln("fcover => %d", get_fcover_str(fcover));)

        //
        // Phase B: Search the F-run to be split.
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase B ==");)

        fnode_type* div_fnode         = fcover.fnode;
        offset_type div_fofst_in_node = fcover.fofst;
        {
            const size_type   lnum_befo = fcover.weight / 2;
            const lnode_type* lnode     = fcover.lnode;
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
        const size_type div_fexp_befo         = div_fofst_befo + 1;

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
            fcover_type fcover = get_fcover(div_fnode);
            div_fnode->get_data().set_weight(fcover.weight);
        }
        {
            fcover_type fcover = get_fcover(new_fnode);
            new_fnode->get_data().set_weight(fcover.weight);
            fcover.lnode->get_data().add_weight(1);  // for new_fnode
        }

        //
        // Phase X: Group split
        //
        DEBUG_PRINT(tfm::printfln("\n== Phase X ==");)

        fnode_type* updated_fnode = nullptr;
        if (div_lnode and div_lnode->get_data().get_num_runs() >= GROUP_BOUND) {
            updated_fnode = divide_lnode(div_lnode);
        }
        return updated_fnode;
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
            div_lnode->get_data().set_weight(lcover.weight);
        }
        {
            const lcover_type lcover = get_lcover(new_lnode);
            new_lnode->get_data().set_weight(lcover.weight);
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
            div_fnode->get_data().set_weight(fcover.weight);
        }
        {
            const fcover_type fcover = get_fcover(new_fnode);
            new_fnode->get_data().set_weight(fcover.weight);
        }

        lnode_type* updated_lnode = std::get<0>(get_first_overlapped_lnode(new_fnode));
        updated_lnode->get_data().add_weight(1);  // for new_fnode

        return updated_lnode;
    }

    // Split L-run (lnode, lcrsr_befo) such that the front part has lexp_befo characters.
    // Note that, if L-run (lnode, lcrsr_befo) has $-mark, lexp_befo should NOT consider $.
    std::pair<fnode_type*, fcursor_type> split_lrun(lnode_type* lnode, lcursor_type lcrsr_befo, size_type lexp_befo) {
        DEBUG_ABORT_IF_EQ(lexp_befo, get_lexp(lnode, lcrsr_befo));  // not need to split?

        auto [fnode, fcrsr_befo] = get_corr_fcursor(lnode, lcrsr_befo);

        lcursor_type lcrsr_aftr = lnode->get_data().split_after(lcrsr_befo, lexp_befo);
        fcursor_type fcrsr_aftr = fnode->get_data().split_after(fcrsr_befo, lnode, lcrsr_aftr.get_lrnid());
        lnode->get_data().set_frptr(lcrsr_aftr, fnode, fcrsr_aftr.get_frnid());

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

    // Split F-run (fnode, fcrsr_befo) such that the front part has fexp_befo characters.
    std::pair<lnode_type*, lcursor_type> split_frun(fnode_type* fnode, fcursor_type fcrsr_befo, size_type fexp_befo) {
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

    // Search the most backword F-node and run in the node such that
    //  - F-run's character is no greater than 'chr',
    //  - F-run's order is no greater than 'lo_common::get_basic_order(q_lnode)', and
    //  - F-run's position in the corresponding L-node is no greater than 'q_lnode->get_data().get_lrpos(q_lrnid)'.
    // Return its pointer and run-ID.
    std::pair<fnode_type*, fcursor_type> predecessor(const uchar_type chr, const lnode_type* q_lnode,
                                                     const lcursor_type q_lcurs) {
        const loint_type q_order        = lo_common::get_basic_order(q_lnode);
        const size_type  q_second_order = q_lcurs.get_lrpos();

        // TODO: If the group contains the same character, this can be used to omit binary search.

        return m_findex.predecessor(chr, q_order, q_second_order);
        // return m_findex.predecessor_naive(chr, q_order, q_second_order);
    }

    // Return (lnode, lrnid, lofst) such that the 'em_fofst'-th character of F-run (em_fnode, em_frnid) corresponds to
    // the 'lofst'-th character of L-run (lnode, lrnid).
    inline std::tuple<lnode_type*, size_type, offset_type> get_lpartner(const fnode_type* q_fnode,
                                                                        const size_type   q_frnid,
                                                                        const offset_type q_fofst) const {
        return get_lpartner(q_fnode, q_fnode->get_data().get_fcursor(q_frnid), q_fofst);
    }
    inline std::tuple<lnode_type*, size_type, offset_type> get_lpartner(const fnode_type*  q_fnode,
                                                                        const fcursor_type q_fcrsr,
                                                                        const offset_type  q_fofst) const {
        DEBUG_ABORT_IF(!q_fnode);
        DEBUG_ABORT_IF(is_head_fnode(q_fnode));

        // Exception case: \Lambda_t indicates the tail character of F.
        {
            auto [tail_fnode, tail_fcrsr] = get_tail_fcursor();
            if (q_fnode == tail_fnode and q_fcrsr.get_frnid() == tail_fcrsr.get_frnid() and
                get_fexp(q_fnode, q_fcrsr) == q_fofst + 1) {
                // Then, 'lofst' to be returned indicates exclusive position.
                auto [tail_lnode, tail_lcrsr] = get_tail_lcursor();
                return {const_cast<lnode_type*>(tail_lnode), tail_lcrsr.get_lrnid(), get_lexp(tail_lnode, tail_lcrsr)};
            }
        }

        auto [lnode, fofst]               = get_first_overlapped_lnode(q_fnode);
        const offset_type q_fofst_in_node = q_fnode->get_data().get_offset(q_fcrsr) + q_fofst;

        while (fofst < q_fofst_in_node) {
            lnode = lnode->get_next();
            fofst += get_sum_lexps(lnode);
        }

        lcursor_type lcrsr     = {};
        offset_type  lofst     = (get_sum_lexps(lnode) - 1) - (fofst - q_fofst_in_node);
        std::tie(lcrsr, lofst) = lnode->get_data().access_with_offset(lofst);

        return {lnode, lcrsr.get_lrnid(), lofst};
    }

    inline std::tuple<fnode_type*, size_type, offset_type> get_fpartner(const lnode_type* q_lnode,
                                                                        const size_type   q_lrnid,
                                                                        const offset_type q_lofst) const {
        return get_fpartner(q_lnode, q_lnode->get_data().get_lcursor(q_lrnid), q_lofst);
    }
    inline std::tuple<fnode_type*, size_type, offset_type> get_fpartner(const lnode_type*  q_lnode,
                                                                        const lcursor_type q_lcrsr,
                                                                        const offset_type  q_lofst) const {
        DEBUG_ABORT_IF(!q_lnode);
        DEBUG_ABORT_IF(is_dmmy_lnode(q_lnode));

        auto [fnode, lofst]         = get_first_overlapped_fnode(q_lnode);
        offset_type q_lofst_in_node = q_lnode->get_data().get_offset(q_lcrsr) + q_lofst;

        // Solution for that get_offset() does not consider $-mark
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

        fcursor_type fcrsr     = {};
        offset_type  fofst     = (get_sum_fexps(fnode) - 1) - (lofst - q_lofst_in_node);
        std::tie(fcrsr, fofst) = fnode->get_data().access_with_offset(fofst);

        return {fnode, fcrsr.get_frnid(), fofst};
    }

    //  Return (lnode, fofst) such that fnode[fofst] corresponds to the back of lnode.
    //
    //    e.g.,
    //               [F]      [L]
    //            v1  a
    //                a        b  u1
    //            v2  a        a
    //                a        c  u2 *      The output is (u2, 1)
    //       (in) v3  b        c            since v3[1] corresponds to u2.back().
    //                b  --->  b
    //                b        a  u3
    //                         b
    //
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
        fofst += get_begin_offset(fnode);

        // Move to the first overlapped L-node with the given F-node
        do {
            fofst += get_sum_lexps(lnode);
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
    //          * v2  a        a
    //                a        c  u2 (in)   The output is (v2, 2)
    //                a        c            since u2[2] corresponds to v2.back.
    //                a  <---  b
    //            v3  b        a  u3
    //                b        b
    //
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
        lofst -= get_begin_offset(fnode);

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

    //  Return the L-nodes overlapped with the F-interval of size 'fsize' starting at the 'fofst'-th of 'fnode'.
    //
    //     e.g.,
    //                     [F]      [L]
    //                  v1  a
    //                      a        b  u1
    //             (in) v2  a        a
    //     fofst = 1 -->    a +    + c  u2 *       The output is [u2, u3].
    //                  v3  b |    | c
    //                      b |    | b
    //     fsize = 4 -->    b +    | a  u3 *
    //                             + b
    //
    // If on_cache, then store the result of get_first_overlapped_lnode in dc_*
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

    //  Return the F-nodes overlapped with the L-interval of size 'lsize' starting at the 'lofst'-th of 'lnode'.
    //
    //     e.g.,
    //                 [F]      [L]
    //              v1  a                          The output is [v2, v3].
    //                  a        b  u1 (in)
    //            * v2  a +      a
    //                  a |    + c  u2  <-- lofst = 2
    //            * v3  b |    | c
    //                  b |    + b      <-- lsize = 3
    //                  b +      a  u3
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
    void remove_lnode_links(lbuffer_type& lnodes) {
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
    void remove_fnode_links(fbuffer_type& fnodes) {
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
        offset_type lofst_fb   = lofst_fe - get_sum_fexps(fnode) + 1;

        if ((lofst_fb <= lofst_lb) and (lofst_lb <= lofst_fe) and (lofst_fe <= lofst_le)) {
            DEBUG_ABORT_IF(fnode->get_data().get_tnode());
            lnode->get_data().set_hlink(fnode, lofst_fe - lofst_lb);
            fnode->get_data().set_tnode(lnode);
        }
    }

    void reset_lnode_links(lbuffer_type& lnodes) {
        for (lnode_type* lnode : lnodes) {
            reset_lnode_link(lnode);
        }
    }

    void reset_lnode_links_in_order(const lbuffer_type& lnodes) {
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
        offset_type lofst_lb  = 0;
        offset_type lofst_le  = sum_lexps - 1;

        // Consider 'lnode' as the pivot, the F-interval is [lofst_fb, lofst_fe].
        auto [fnode, lofst_fe] = get_first_overlapped_fnode(lnode);
        offset_type sum_fexps  = get_sum_fexps(fnode);
        offset_type lofst_fb   = lofst_fe - sum_fexps + 1;

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
                lofst_lb  = 0;
                lofst_le  = sum_lexps - 1;
            }

            if (fskip) {
                fnode = fnode->get_next();

                sum_fexps = get_sum_fexps(fnode);
                lofst_fe  = lofst_fe + sum_fexps;
                lofst_fb  = lofst_fe - sum_fexps + 1;
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
        offset_type fofst_lb   = fofst_le - get_sum_lexps(lnode) + 1;

        while (fofst_lb <= fofst_fe) {
            if ((fofst_fb <= fofst_lb) and (fofst_lb <= fofst_fe) and (fofst_fe <= fofst_le)) {
                DEBUG_ABORT_IF(lnode->get_data().get_hnode());
                lnode->get_data().set_hlink(fnode, fofst_fe - fofst_lb);
                fnode->get_data().set_tnode(lnode);
                break;
            }
            lnode    = lnode->get_next();
            fofst_lb = fofst_le + 1;
            fofst_le = fofst_le + get_sum_lexps(lnode);
        }
    }

    void reset_fnode_links(fbuffer_type& fnodes) {
        for (fnode_type* fnode : fnodes) {
            reset_fnode_link(fnode);
        }
    }

    void reset_fnode_links_in_order(const fbuffer_type& fnodes) {
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
        offset_type fofst_fb  = 0;
        offset_type fofst_fe  = sum_fexps - 1;

        // Consider 'lnode' as the pivot, the F-interval is [lofst_fb, lofst_fe].
        auto [lnode, fofst_le] = get_first_overlapped_lnode(fnode);
        offset_type sum_lexps  = get_sum_lexps(lnode);
        offset_type fofst_lb   = fofst_le - sum_lexps + 1;

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
                fofst_fb  = 0;
                fofst_fe  = sum_fexps - 1;
            }

            if (lskip) {
                lnode = lnode->get_next();

                sum_lexps = get_sum_lexps(lnode);
                fofst_le  = fofst_le + sum_lexps;
                fofst_lb  = fofst_le - sum_lexps + 1;
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
    inline lcover_type get_lcover(const lnode_type* lnode) const {
        auto [fnode, lofst_e] = get_first_overlapped_fnode(lnode);
        offset_type lofst_b   = lofst_e - get_sum_fexps(fnode) + 1;

        // Is the first overlapped F-node not covered? (i.e., the head is not overlapped)
        if (lofst_b < 0) {
            fnode   = fnode->get_next();
            lofst_b = lofst_e + 1;
            lofst_e = lofst_b + get_sum_fexps(fnode) - 1;
        }

        lcover_type       ret  = {const_cast<lnode_type*>(lnode), fnode, lofst_e, 0};
        const offset_type lexp = get_sum_lexps(lnode);

        while (lofst_b < lexp) {
            ret.weight += 1;
            lofst_b += get_sum_fexps(fnode);
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
    inline fcover_type get_fcover(const fnode_type* fnode) const {
        auto [lnode, fofst_e] = get_first_overlapped_lnode(fnode);
        offset_type fofst_b   = fofst_e - get_sum_lexps(lnode) + 1;

        // Is the first overlapped L-node not covered? (i.e., the head is not overlapped)
        if (fofst_b < 0) {
            lnode   = lnode->get_next();
            fofst_b = fofst_e + 1;
            fofst_e = fofst_b + get_sum_lexps(lnode) - 1;
        }

        fcover_type       ret  = {const_cast<fnode_type*>(fnode), lnode, fofst_e, 0};
        const offset_type fexp = get_sum_fexps(fnode);

        while (fofst_b < fexp) {
            ret.weight += 1;
            fofst_b += get_sum_lexps(lnode);
            lnode = lnode->get_next();
            if (is_dmmy_lnode(lnode)) {  // cycled?
                break;
            }
        }
        return ret;
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
    std::string get_pc(const void* p) const {
        return tfm::format("[%c]", m_pcmap.find(intptr_t(p))->second);
    }

    std::string get_lnode_str(const lnode_type* p) {
        if (!p) return tfm::format("(%s)", get_pc(p));
        const auto& data = p->get_data();
        std::string str  = tfm::format("%s => H=%s, O=%d, W=%d",  //
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
        std::string str  = tfm::format("%s => T=%s, W=%d", get_pc(p), get_pc(data.get_tnode()), data.get_weight());
        for (size_type i = 0; i < data.get_num_runs(); i++) {
            const fcursor_type crsr = data.get_fcursor_from_position(i);
            str += tfm::format("\n\t%s", get_frun_str(p, crsr));
        }
        return str;
    }

    std::string get_lrun_str(const lnode_type* p, const lcursor_type crsr) {
        if (!p) return "()";
        const auto& data = p->get_data();
        return tfm::format("(R=%s, F=(%s,%d), L=%d)", make_run(data.get_chr(crsr), get_lexp(p, crsr)),
                           get_pc(data.get_frptr(crsr).first), data.get_frptr(crsr).second, crsr.get_lrnid());
    }
    std::string get_frun_str(const fnode_type* p, const fcursor_type crsr) {
        if (!p) return "()";
        const auto& data = p->get_data();
        return tfm::format("(R=%s, L=(%s,%d), F=%d)", make_run(data.get_chr(crsr), get_fexp(p, crsr)),
                           get_pc(data.get_lrptr(crsr).first), data.get_lrptr(crsr).second, crsr.get_frnid());
    }

    std::string get_lcrsr_str(const lcursor_type crsr) {
        return tfm::format("(%d,%d)", crsr.get_lrnid(), crsr.get_lrpos());
    }
    std::string get_fcrsr_str(const fcursor_type crsr) {
        return tfm::format("(%d,%d)", crsr.get_frnid(), crsr.get_frpos());
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