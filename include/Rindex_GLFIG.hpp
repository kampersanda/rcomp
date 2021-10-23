/**
 * @file Rindex_GLFIG.hpp
 */
#pragma once

#include "GroupedLFIntervalGraph.hpp"
#include "Successor.hpp"
#include "utils.hpp"

#include "Rindex_Interface.hpp"

#ifdef ENABLE_STAT_MONITOR
#include "Timer.hpp"
#endif

namespace rcomp {

/**
 * An r-index implementation with the grouped LF-interval graph.
 *
 * @tparam t_Graph The grouped LF-interval graph.
 */
template <class t_Graph>
class Rindex_GLFIG : public Rindex_Interface {
  public:
    using this_type = Rindex_GLFIG<t_Graph>;
    using graph_type = t_Graph;
    using lnode_type = typename graph_type::lnode_type;
    using fnode_type = typename graph_type::fnode_type;
    using lcursor_type = typename lnode_type::data_type::lcursor_type;
    using fcursor_type = typename fnode_type::data_type::fcursor_type;
    using successor_type = Successor<size_type, size_type>;
    using literator_type = typename graph_type::literator_type;
    using fiterator_type = typename graph_type::fiterator_type;

    using insert_etype = typename graph_type::insert_etype;
    using bwres_etype = typename graph_type::bwres_etype;

    static_assert(t_Graph::ldata_type::WITH_SAE);

  private:
    // Grouped LF-interval graph
    graph_type m_graph;

    // Successor index storing SA-entries on the L-run boundary.
    // This data structure is equivalent to pred_B in paper "Refining the r-index".
    // In our implementation, it's "Successor" since the SA-entries are reversed.
    successor_type m_sae_succ;

    // Temporal SA-entries around $-marker in L.
    size_type m_prev_sae = MAX_SIZE_INT;  // SA[i-1] for $'s position i
    size_type m_curr_sae = MAX_SIZE_INT;  // SA[i]   for $'s position i
    size_type m_next_sae = MAX_SIZE_INT;  // SA[i+1] for $'s position i

#ifdef ENABLE_STAT_MONITOR
    Timer m_timer_4_case;
    Timer m_timer_4_phase;

    size_type m_num_case_C = 0;
    size_type m_num_case_B = 0;
    size_type m_num_case_LF = 0;

    double m_ns_case_C = 0.0;
    double m_ns_case_B = 0.0;
    double m_ns_case_LF = 0.0;
#endif

  public:
    //! Default constructor
    Rindex_GLFIG() = default;

    //! Default destructor
    virtual ~Rindex_GLFIG() = default;

    //! Copy constructor (deleted)
    Rindex_GLFIG(const Rindex_GLFIG&) = delete;

    //! Copy constructor (deleted)
    Rindex_GLFIG& operator=(const Rindex_GLFIG&) = delete;

    //! Move constructor
    Rindex_GLFIG(Rindex_GLFIG&&) noexcept = default;

    //! Move constructor
    Rindex_GLFIG& operator=(Rindex_GLFIG&&) noexcept = default;

    //! Check if the data structure is empty.
    inline bool is_empty() const override {
        return m_graph.is_empty();
    }

    //! Get the number of stored characters.
    inline size_type get_num_chars() const override {
        return m_graph.get_num_chars();
    }

    //! Get the allocated memory in bytes.
    size_type get_memory_in_bytes(bool include_this = true) const override {
        const size_type this_bytes = sizeof(*this) * include_this;
        return this_bytes + m_graph.get_memory_in_bytes(false) + m_sae_succ.get_memory_in_bytes(false);
    }

    //! Print the statistics related to memory (via stdout).
    void show_memory_statistics() const override {
        m_graph.show_memory_statistics();
    }

    //! Print the detailed statistics (via stdout).
    void show_detailed_statistics() const override {
        m_graph.show_detailed_statistics();
    }

    //! Print the statistics measured with internal monitors (via stdout).
    void show_monitored_statistics() const override {
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
#endif
    }

    //! Extend the RLBWT text by appending the given character (i.e., T := c + T).
    void extend(const uchar_type new_chr) override {
        ABORT_IF_LE(new_chr, END_MARKER);

        DEBUG_PRINT(tfm::printfln("");)
        DEBUG_PRINT(tfm::printfln("*********************************");)
        DEBUG_PRINT(tfm::printfln("     extend(chr=%c, size=%d)", new_chr, get_num_chars());)
        DEBUG_PRINT(tfm::printfln("*********************************");)
        DEBUG_PRINT(tfm::printfln("");)

        if (is_empty()) {
            m_graph.extend_init(new_chr);

            m_graph.set_sae(m_graph.get_em_literator(), 0);
            m_sae_succ.insert(0, 1);
            m_sae_succ.insert(1, MAX_SIZE_INT);
            m_prev_sae = 0;
            m_curr_sae = 1;
            m_next_sae = MAX_SIZE_INT;

            DEBUG_PRINT(debug_print());
            return;
        }

        // Case L or F
        STAT_MONITOR(m_timer_4_case.start();)
        [[maybe_unused]] size_type num_divided = 0;
        while (divide_fat_node()) {
            num_divided += 1;
        }
        STAT_MONITOR(m_ns_case_LF += m_timer_4_case.stop_and_get_ns();)
        STAT_MONITOR(m_num_case_LF += num_divided;)

        // Phase A
        STAT_MONITOR(m_timer_4_case.start();)
        auto ins_mode = m_graph.check_mergable(new_chr);

        // Remove successor entries that are no longer needed due to the merge.
        {
            auto em_liter = m_graph.get_em_literator();
            if (ins_mode == insert_etype::MERGE_PREV_IN or ins_mode == insert_etype::MERGE_PREV_OUT) {
                m_sae_succ.remove(m_prev_sae);
            } else if (ins_mode == insert_etype::MERGE_CURR) {
                if (m_graph.is_first_char_in_lrun(em_liter)) {
                    m_sae_succ.remove(m_curr_sae);
                } else if (m_graph.is_last_char_in_lrun(em_liter)) {
                    m_sae_succ.remove(m_prev_sae);
                } else {
                    m_sae_succ.remove(m_prev_sae);
                    m_sae_succ.remove(m_curr_sae);
                }
            } else {
                DEBUG_ABORT_IF(m_graph.is_mergable_mode(ins_mode));
            }
        }

        // Phase B--F
        if (m_graph.is_mergable_mode(ins_mode)) {
            m_graph.extend_with_merge(new_chr, ins_mode);

            // If the new character is the last one of the merged L-node, need to update the SA-entry.
            auto em_fiter = m_graph.get_em_fiterator();
            if (m_graph.is_last_char_in_frun(em_fiter)) {
                m_graph.set_sae(m_graph.get_corr_lcursor(em_fiter), m_curr_sae);
            }

            STAT_MONITOR(m_ns_case_C += m_timer_4_case.stop_and_get_ns(););
            STAT_MONITOR(m_num_case_C += 1;)
        } else {
            // Temporarily keep the SA-entry of $'s L-node.
            const uint64_t em_sae = m_graph.get_sae(m_graph.get_em_literator());

            // Divide
            m_graph.extend_with_split(new_chr, ins_mode);

            // Set the SA-entry for the new inserted L-node
            auto em_fiter = m_graph.get_em_fiterator();
            auto ins_liter = m_graph.get_corr_lcursor(em_fiter);

            // Set the old $'s L-run.
            m_graph.set_sae(ins_liter, m_curr_sae);

            if (ins_mode == insert_etype::SPLIT_FRONT or ins_mode == insert_etype::SPLIT_MIDDLE) {
                auto next_liter = ins_liter;
                m_graph.set_next_lcursor(next_liter);
                m_graph.set_sae(next_liter, em_sae);
            }

            if (ins_mode == insert_etype::SPLIT_BACK or ins_mode == insert_etype::SPLIT_MIDDLE) {
                auto prev_liter = ins_liter;
                m_graph.set_prev_lcursor(prev_liter);
                m_graph.set_sae(prev_liter, m_prev_sae);
            }

            STAT_MONITOR(m_ns_case_B += m_timer_4_case.stop_and_get_ns();)
            STAT_MONITOR(m_num_case_B += 1;)
        }

        {
            // Increment the SA-entry for $'s character in L
            m_curr_sae += 1;

            auto em_fiter = m_graph.get_em_fiterator();

            // Update the previous SA-entry
            if (m_graph.is_first_char_in_frun(em_fiter)) {
                auto prev_fiter = em_fiter;
                m_graph.set_prev_fcursor(prev_fiter);
                if (m_graph.is_head_fcursor(prev_fiter)) {
                    m_prev_sae = 0;
                } else {
                    m_prev_sae = get_last_sae(prev_fiter) + 1;
                }
            } else {
                m_prev_sae += 1;
            }

            // Since the result may be used in get_first_sae() of the next step,
            // we need to do the update here
            m_sae_succ.insert(m_prev_sae, m_curr_sae);

            // Update the next SA-entry
            if (m_graph.is_last_char_in_frun(em_fiter)) {
                auto next_fiter = em_fiter;
                m_graph.set_next_fcursor(next_fiter);
                if (!next_fiter.is_valid()) {
                    m_next_sae = MAX_SIZE_INT;
                } else {
                    m_next_sae = get_first_sae(next_fiter) + 1;
                }
            } else {
                m_next_sae += 1;
            }

            m_sae_succ.insert(m_curr_sae, m_next_sae);
        }

        DEBUG_PRINT(debug_print());
    }

    //! Count the number of occurrences for the query pattern.
    size_type count(range_type<const uchar_type*> pat) const override {
        auto [lo_liter, hi_liter] = bwsearch(pat);
        if (!lo_liter.is_valid()) {
            return 0;
        }
        size_type occ = 1;
        while (!m_graph.is_equal_literator(lo_liter, hi_liter)) {
            occ += 1;
            m_graph.set_next_literator(lo_liter);
        }
        return occ;
    }

    //! Report the positions of occurrence for the query pattern via the callback function.
    void locate(range_type<const uchar_type*> pat, const std::function<void(size_type)>& fn) const override {
        auto [lo_liter, hi_liter, sae] = bwsearch_toehold(pat);
        if (!lo_liter.is_valid()) {
            return;
        }
        fn(sae);
        while (!m_graph.is_equal_literator(lo_liter, hi_liter)) {
            if (!set_next_sae_and_iterator(sae, lo_liter)) {
                break;
            }
            fn(sae);
        }
    }

    //! Extract SA-entries from the head of the BWT-text.
    void extract_sa_entries(const std::function<void(size_type)>& fn) const override {
        auto liter = m_graph.get_head_literator();
        for (size_type sae = 0;;) {
            fn(sae);
            if (!set_next_sae_and_iterator(sae, liter)) {
                break;
            }
        }
    }

    //! Output the original text (in the input order) via the callback function.
    void decode_text(const std::function<void(uchar_type)>& fn) const override {
        m_graph.decode_text(fn);
    }

    //! Output the RLBWT text via the callback function.
    void output_runs(const std::function<void(const run_type&)>& fn) const override {
        m_graph.output_runs(fn);
    }

    //! Test the data structure.
    void test() const override {
        m_graph.test_alphabet_order();
        m_graph.test_list_order();
        m_graph.test_lf_mapping();
        m_graph.test_ht_links();
        m_graph.test_weights();
    }

#ifdef ENABLE_DEBUG_PRINT
    void debug_print() {
        m_graph.debug_print();

        tfm::printfln("");
        tfm::printfln("**** Rindex_GLFIG::debug_print() ****");
        tfm::printfln("");

        tfm::printfln("m_prev_sae=%d", m_prev_sae);
        tfm::printfln("m_curr_sae=%d", m_curr_sae);
        tfm::printfln("m_next_sae=%d", m_next_sae);

        tfm::printfln("m_sae_succ=");
        for (const auto& kv : m_sae_succ) {
            tfm::printfln(" - %d => %d", kv.first, kv.second);
        }
    }
#endif

  private:
    //! Get the SA value of the first character in the given L-run (in log time).
    inline size_type get_first_sae(const literator_type& liter) const {
        if (m_graph.is_head_lcursor(liter)) {
            return 0;  // SA'[0] is always 0.
        }
        if (m_graph.get_lexp(liter) == 1) {
            return m_graph.get_sae(liter);
        }

        auto prev_liter = liter;
        m_graph.set_prev_lcursor(prev_liter);

        auto res = m_sae_succ.find(m_graph.get_sae(prev_liter), MAX_SIZE_INT);
        ABORT_IF_EQ(res, MAX_SIZE_INT);
        return res;
    }

    //! Get the SA value of the first character in the corresponding L-run (in log time).
    inline size_type get_first_sae(const fiterator_type& fiter) const {
        return get_first_sae(m_graph.get_corr_lcursor(fiter));
    }

    //! Get the SA value of the last character in the given L-run (in constant time).
    inline size_type get_last_sae(const literator_type& liter) const {
        return m_graph.get_sae(liter);
    }

    //! Get the SA value of the last character in the corresponding L-run (in constant time).
    inline size_type get_last_sae(const fiterator_type& fiter) const {
        return get_last_sae(m_graph.get_corr_lcursor(fiter));
    }

    //! Get SA[i+1] where SA[i] is the given value.
    inline size_type get_next_sae(size_type sae) const {
        const auto [k, v] = m_sae_succ.search(sae, MAX_SIZE_INT);
        if (k == MAX_SIZE_INT) {
            return MAX_SIZE_INT;
        }
        return v - (k - sae);
    }

    inline bool set_next_sae_and_iterator(size_type& sae, literator_type& liter) const {
        sae = get_next_sae(sae);
        if (sae == MAX_SIZE_INT) {
            return false;
        }
        m_graph.set_next_literator(liter);
        return true;
    }

    // Backward search returns (lo_lnode, lo_lcrsr, lo_lofst, hi_lnode, hi_lcrsr, hi_lofst)
    // LOWER and UPPER BOUND
    std::tuple<literator_type, literator_type> bwsearch(range_type<const uchar_type*> pat) const {
        DEBUG_ABORT_IF(pat.beg >= pat.end);

        // Search range
        literator_type lo_liter = m_graph.get_head_literator();
        literator_type hi_liter = m_graph.get_tail_literator();
        fiterator_type lo_fiter, hi_fiter;

        bwres_etype res;
        for (auto itr = pat.beg; itr < pat.end; ++itr) {
            const uchar_type chr = *itr;

            // Update the F lower bound
            res = m_graph.successor_in_bwsearch(chr, lo_liter, lo_fiter);
            if (res == bwres_etype::FAILED) {
                return {};
            }

            // Update the upper bound
            res = m_graph.predecessor_in_bwsearch(chr, hi_liter, hi_fiter);
            if (res == bwres_etype::FAILED) {
                return {};
            }

            lo_liter = m_graph.get_overlapped_literator(lo_fiter);
            hi_liter = m_graph.get_overlapped_literator(hi_fiter);

            // Check if the search range is valid
            if (m_graph.compare_order(hi_liter, lo_liter)) {
                return {};
            }
        }

        return {lo_liter, hi_liter};
    }

    std::tuple<literator_type, literator_type, size_type> bwsearch_toehold(range_type<const uchar_type*> pat) const {
        DEBUG_ABORT_IF(pat.beg >= pat.end);

        DEBUG_PRINT(tfm::printfln("****** bwsearch_toehold ******"));

        // Search range
        literator_type lo_liter = m_graph.get_head_literator();
        literator_type hi_liter = m_graph.get_tail_literator();
        fiterator_type lo_fiter, hi_fiter;

        bwres_etype res;
        size_type toehold_sae = 0;  // SA entry of the head of L-list

        for (auto itr = pat.beg; itr < pat.end; ++itr) {
            const uchar_type chr = *itr;

            DEBUG_PRINT(tfm::printfln("* chr=%c *", chr));
            DEBUG_PRINT(tfm::printfln("lo_liter=(%s)", m_graph.get_litr_str(lo_liter)));
            DEBUG_PRINT(tfm::printfln("hi_liter=(%s)", m_graph.get_litr_str(hi_liter)));

            // Update the F lower bound
            res = m_graph.successor_in_bwsearch(chr, lo_liter, lo_fiter);
            if (res == bwres_etype::FAILED) {
                return {};
            }

            // DEBUG_PRINT(tfm::printfln("lo_fiter=(%s)", m_graph.get_fitr_str(lo_fiter)));
            DEBUG_PRINT(tfm::printfln("lo_fiter.fnode=(%s)", m_graph.get_pc(lo_fiter.fnode)));
            DEBUG_PRINT(tfm::printfln("lo_fiter.fcrsr.frnid=(%d)", lo_fiter.fcrsr.get_frnid()));
            DEBUG_PRINT(tfm::printfln("lo_fiter.fcrsr.frpos=(%d)", lo_fiter.fcrsr.get_frpos()));

            if (res == bwres_etype::MATCHED) {
                toehold_sae += 1;
            } else if (res == bwres_etype::NOT_MATCHED_BY_EM) {
                // SPECIAL CASE:
                // Since the itr indicates $-marker, the successor is the one simply shifted in the L-node.
                toehold_sae = m_next_sae + 1;
            } else {
                ABORT_IF(!m_graph.is_first_char_in_frun(lo_fiter));

                literator_type co_liter = m_graph.get_corr_lcursor(lo_fiter);
                co_liter.lofst = 0;

                if (m_graph.is_head_lcursor(co_liter)) {
                    // Then, the SA-entry is always zero (for $)
                    toehold_sae = 1;  // 0+1
                } else if (m_graph.is_em_literator(co_liter)) {
                    // Then, the corr L character is the next of $, and we have the SA-entry in m_next_sae.
                    toehold_sae = m_next_sae + 1;
                } else if (m_graph.get_lexp(co_liter) == 1) {
                    // Then, the SA-entry is the last one in the L-node
                    toehold_sae = get_last_sae(co_liter) + 1;
                } else {
                    // We extract corr_lnode[0].sae using Successor query for the last SA of the previous L-node
                    // auto corr_lprev = corr_lnode->get_prev();
                    m_graph.set_prev_lcursor(co_liter);
                    const auto [k, v] = m_sae_succ.search(get_last_sae(co_liter), MAX_SIZE_INT);
                    DEBUG_ABORT_IF_NE(k, get_last_sae(co_liter));
                    toehold_sae = v + 1;
                }
            }

            // Update the upper bound
            res = m_graph.predecessor_in_bwsearch(chr, hi_liter, hi_fiter);
            if (res == bwres_etype::FAILED) {
                return {};
            }

            DEBUG_PRINT(tfm::printfln("hi_fiter=(%s)", m_graph.get_fitr_str(hi_fiter)));

            lo_liter = m_graph.get_overlapped_literator(lo_fiter);
            hi_liter = m_graph.get_overlapped_literator(hi_fiter);

            // Check if the search range is valid
            if (m_graph.compare_order(hi_liter, lo_liter)) {
                return {};
            }
        }

        DEBUG_PRINT(tfm::printfln("lo_liter=(%s)", m_graph.get_litr_str(lo_liter)));
        DEBUG_PRINT(tfm::printfln("hi_liter=(%s)", m_graph.get_litr_str(hi_liter)));

        return {lo_liter, hi_liter, toehold_sae};
    }

    //! Divide the fat nodes in Case L/F, while updating SA-entries.
    inline bool divide_fat_node() {
        lnode_type* div_lnode = nullptr;
        if (m_graph.divide_fat_lnode(div_lnode)) {
            update_sa_entries_in_L(div_lnode);
            return true;
        }
        fnode_type* div_fnode = nullptr;
        if (m_graph.divide_fat_fnode(div_fnode)) {
            update_sa_entries_in_F(div_fnode);
            return true;
        }
        return false;
    }

    /**
     * @brief Update the SA-entries of L-runs split in Case L.
     * @param[in] div_lnode The before part of the divided L-nodes in Case L.
     * @note Assume that
     *   - The boundary of the divided L-nodes does not overlap a F-node.
     *   - The after part of the split L-runs is the one newly created and does not have an SA-entry.
     *   - The before part of the split L-runs has the SA-entry of the after part.
     */
    inline void update_sa_entries_in_L(lnode_type* div_lnode) {
        if (div_lnode == nullptr) {
            return;
        }

        // div_lnode and new_lnode are the divided L-nodes in Case L.
        auto new_lnode = div_lnode->get_next();

        // 1) Get the overlapped L-nodes.
        DEBUG_ABORT_IF(!new_lnode->get_data().get_hnode());
        auto new_ovlp_fnode = std::get<0>(m_graph.get_first_overlapped_fnode(new_lnode));
        auto div_ovlp_fnode = new_ovlp_fnode->get_prev();

        // 2) Get the split L-runs.
        auto div_lcrsr = div_lnode->get_data().get_last_lcursor();
        auto new_lcrsr = new_lnode->get_data().get_first_lcursor();

        // 3) Set the SA-entry of the after L-run.
        const auto new_sae = get_last_sae({div_lnode, div_lcrsr});
        new_lnode->get_data().set_sae(new_lcrsr, new_sae);

        // 4) Set the SA-entry of the before L-run via LF-mapping.
        auto div_ovlp_fcrsr = div_ovlp_fnode->get_data().get_last_fcursor();
        const auto div_sae = get_last_sae({div_ovlp_fnode, div_ovlp_fcrsr}) + 1;
        div_lnode->get_data().set_sae(div_lcrsr, div_sae);

        // 4) Update successor
        const auto next_sae = get_next_sae(div_sae);
        m_sae_succ.insert(div_sae, next_sae);
    }

    /**
     * @brief Update the SA-entries of L-runs split in Case F.
     * @param[in] div_fnode The before part of the divided F-nodes in Case F.
     * @note Assume that
     *   - The boundary of the divided F-nodes does not overlap a L-node.
     *   - The split L-runs are in the same L-node.
     *   - The after part of the split L-runs is the one newly created and does not have an SA-entry.
     *   - The before part of the split L-runs has the SA-entry of the after part.
     */
    inline void update_sa_entries_in_F(fnode_type* div_fnode) {
        if (div_fnode == nullptr) {
            return;
        }

        // div_fnode and new_fnode are the divided F-nodes in Case F.
        auto new_fnode = div_fnode->get_next();

        // 1) Get the overlapped F-nodes.
        DEBUG_ABORT_IF(!new_fnode->get_data().get_tnode());
        auto new_ovlp_lnode = std::get<0>(m_graph.get_first_overlapped_lnode(new_fnode));
        auto div_ovlp_lnode = new_ovlp_lnode->get_prev();

        // 2) Get the split L-runs.
        auto div_fiter = m_graph.get_last_fcursor(div_fnode);
        auto div_corr_liter = m_graph.get_corr_lcursor(div_fiter);
        auto new_corr_liter = m_graph.get_next_lcursor(div_corr_liter);
        DEBUG_ABORT_IF(div_corr_liter.lnode != new_corr_liter.lnode);

        // 3) Set the SA-entry of the after L-run.
        const auto new_sae = get_last_sae(div_corr_liter);
        m_graph.set_sae(new_corr_liter, new_sae);

        // 4) Set the SA-entry of the before L-run via LF-mapping.
        const auto div_ovlp_liter = m_graph.get_last_lcursor(div_ovlp_lnode);
        const auto div_sae = get_last_sae(div_ovlp_liter) - 1;
        m_graph.set_sae(div_corr_liter, div_sae);

        // 5) Update successor
        const auto next_sae = get_next_sae(div_sae);
        m_sae_succ.insert(div_sae, next_sae);
    }
};

}  // namespace rcomp