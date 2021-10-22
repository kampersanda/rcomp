/**
 * @file Rindex_LFIG.hpp
 */
#pragma once

#include "LFIntervalGraph.hpp"
#include "Successor.hpp"
#include "utils.hpp"

#include "Rindex_Interface.hpp"

#ifdef ENABLE_STAT_MONITOR
#include "Timer.hpp"
#endif

namespace rcomp {

/**
 * An r-index class with the straightforward LF-interval graph.
 *
 * @tparam t_Graph The LF-interval graph.
 */
template <class t_Graph>
class Rindex_LFIG : public Rindex_Interface {
  public:
    using this_type = Rindex_LFIG<t_Graph>;
    using graph_type = t_Graph;
    using lnode_type = typename graph_type::lnode_type;
    using fnode_type = typename graph_type::fnode_type;
    using insert_modes = typename graph_type::insert_modes;
    using successor_type = Successor<size_type, size_type>;

    static_assert(t_Graph::ldata_type::WITH_SAE);

  private:
    // LF-interval graph
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
    Rindex_LFIG() = default;

    //! Default destructor
    virtual ~Rindex_LFIG() = default;

    //! Copy constructor (deleted)
    Rindex_LFIG(const Rindex_LFIG&) = delete;

    //! Copy constructor (deleted)
    Rindex_LFIG& operator=(const Rindex_LFIG&) = delete;

    //! Move constructor
    Rindex_LFIG(Rindex_LFIG&&) noexcept = default;

    //! Move constructor
    Rindex_LFIG& operator=(Rindex_LFIG&&) noexcept = default;

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
            auto [lnode, ldmmy] = m_graph.extend_init(new_chr);
            lnode->get_data().set_sae(0);
            ldmmy->get_data().set_sae(MAX_SIZE_INT);
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
            auto [em_lnode, em_lofst] = m_graph.get_em_laddress();

            if (ins_mode == insert_modes::MERGE_PREV) {
                m_sae_succ.remove(m_prev_sae);
            } else if (ins_mode == insert_modes::MERGE_CURR) {
                if (m_graph.is_first_laddress(em_lnode, em_lofst)) {
                    m_sae_succ.remove(m_curr_sae);
                } else if (m_graph.is_last_laddress(em_lnode, em_lofst)) {
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
            auto [em_fnode, em_fofst] = m_graph.get_em_faddress();
            if (m_graph.is_last_faddress(em_fnode, em_fofst)) {
                auto ins_lnode = em_fnode->get_data().get_lnode();
                ins_lnode->get_data().set_sae(m_curr_sae);
            }

            STAT_MONITOR(m_ns_case_C += m_timer_4_case.stop_and_get_ns(););
            STAT_MONITOR(m_num_case_C += 1;)
        } else {
            // In the case, the after part of the divided L-nodes is the one newly created.
            // The new L-node will not have SA-entry and needs to reset it.
            if (ins_mode == insert_modes::SPLIT_MIDDLE) {
                // First, we temporarily store the SA-entry of $'s L-node.
                auto em_lnode = std::get<0>(m_graph.get_em_laddress());
                auto em_sae = em_lnode->get_data().get_sae();
                // Second, we devide $'s L-node into new two L-nodes.
                // Then, em_lnode will be the one newly created (having $ in the first).
                ins_mode = m_graph.divide_if_need(ins_mode);
                em_lnode = std::get<0>(m_graph.get_em_laddress());
                // Third, we set the appropriate SA-entry.
                // Note that, then, the before part has an invalid SA-entry.
                em_lnode->get_data().set_sae(em_sae);
            }

            m_graph.extend_with_split(new_chr, ins_mode);

            // Set the SA-entry for the new inserted L-node
            auto em_fnode = std::get<0>(m_graph.get_em_faddress());
            auto ins_lnode = em_fnode->get_data().get_lnode();
            ins_lnode->get_data().set_sae(m_curr_sae);

            // The previous L-node that may have an invalid SA-entry due to SPLIT_MIDDLE
            auto ins_lprev = ins_lnode->get_prev();
            ins_lprev->get_data().set_sae(m_prev_sae);  // may not have an invalid SA-entry

            STAT_MONITOR(m_ns_case_B += m_timer_4_case.stop_and_get_ns();)
            STAT_MONITOR(m_num_case_B += 1;)
        }

        {
            auto [em_fnode, em_fofst] = m_graph.get_em_faddress();

            // Increment the SA-entry for $'s character in L
            m_curr_sae += 1;

            // Update the previous SA-entry
            if (m_graph.is_first_faddress(em_fnode, em_fofst)) {
                auto fprev = em_fnode->get_prev();
                if (m_graph.is_head_fnode(fprev)) {
                    m_prev_sae = 0;
                } else {
                    m_prev_sae = get_last_sae(fprev) + 1;
                }
            } else {
                m_prev_sae += 1;
            }

            // Since the result may be used in get_first_sae() of the next step,
            // we need to do the update here
            m_sae_succ.insert(m_prev_sae, m_curr_sae);

            // Update the next SA-entry
            if (m_graph.is_last_faddress(em_fnode, em_fofst)) {
                auto fnext = em_fnode->get_next();
                if (m_graph.is_head_fnode(fnext)) {  // cycled?
                    m_next_sae = MAX_SIZE_INT;
                } else {
                    m_next_sae = get_first_sae(fnext) + 1;
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
        auto [lo_lnode, lo_lofst, hi_lnode, hi_lofst] = bwsearch(pat);
        if (lo_lnode == nullptr) {
            return 0;
        }

        DEBUG_ABORT_IF(lo_common::get_basic_order(hi_lnode) < lo_common::get_basic_order(lo_lnode));
        DEBUG_ABORT_IF(hi_lnode == lo_lnode and hi_lofst < lo_lofst);

        size_type occ = 1;
        while (hi_lnode != lo_lnode or hi_lofst != lo_lofst) {
            occ += 1;
            std::tie(lo_lnode, lo_lofst) = get_next(lo_lnode, lo_lofst);
            DEBUG_ABORT_IF(lo_lnode == nullptr);
        }
        return occ;
    }

    //! Report the positions of occurrence for the query pattern via the callback function.
    void locate(range_type<const uchar_type*> pat, const std::function<void(size_type)>& fn) const override {
        auto [lo_lnode, lo_lofst, hi_lnode, hi_lofst, sae] = bwsearch_toehold(pat);
        if (lo_lnode == nullptr) {
            return;
        }

        DEBUG_ABORT_IF(lo_common::get_basic_order(hi_lnode) < lo_common::get_basic_order(lo_lnode));
        DEBUG_ABORT_IF(hi_lnode == lo_lnode and hi_lofst < lo_lofst);

        while (hi_lnode != lo_lnode or hi_lofst != lo_lofst) {
            fn(sae);
            std::tie(lo_lnode, lo_lofst, sae) = get_next_sae(lo_lnode, lo_lofst, sae);
            DEBUG_ABORT_IF(lo_lnode == nullptr);
        }
        fn(sae);
    }

    //! Extract SA-entries from the head of the BWT-text.
    void extract_sa_entries(const std::function<void(size_type)>& fn) const override {
        const lnode_type* lnode = nullptr;
        offset_type lofst = 0;
        size_type sae = 0;

        std::tie(lnode, lofst) = m_graph.get_head_laddress();

        while (lnode != nullptr and sae != MAX_SIZE_INT) {
            fn(sae);
            std::tie(lnode, lofst, sae) = get_next_sae(lnode, lofst, sae);
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
        tfm::printfln("**** Rindex_LFIG.debug_print() ****");
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
    inline size_type get_first_sae(const lnode_type* lnode) const {
        if (lnode == m_graph.get_head_lnode()) {
            return 0;  // SA'[0] is always 0.
        }
        if (m_graph.get_lexp(lnode) == 1) {
            return lnode->get_data().get_sae();
        }
        auto lprev = lnode->get_prev();
        auto result = m_sae_succ.find(lprev->get_data().get_sae(), MAX_SIZE_INT);
        ABORT_IF_EQ(result, MAX_SIZE_INT);
        return result;
    }

    //! Get the SA value of the first character in the corresponding L-run (in log time).
    inline size_type get_first_sae(const fnode_type* fnode) const {
        return get_first_sae(fnode->get_data().get_lnode());
    }

    //! Get the SA value of the last character in the given L-run (in constant time).
    inline size_type get_last_sae(const lnode_type* lnode) const {
        return lnode->get_data().get_sae();
    }

    //! Get the SA value of the last character in the corresponding L-run (in constant time).
    inline size_type get_last_sae(const fnode_type* fnode) const {
        return get_last_sae(fnode->get_data().get_lnode());
    }

    // Backward search returns (lo_lnode, lo_lofst, hi_lnode, hi_lofst)
    inline std::tuple<const lnode_type*, offset_type, const lnode_type*, offset_type>  //
    bwsearch(range_type<const uchar_type*> pat) const {
        DEBUG_ABORT_IF(pat.beg >= pat.end);

        // Search range
        auto [lo_lnode, lo_lofst] = m_graph.get_head_laddress();
        auto [hi_lnode, hi_lofst] = m_graph.get_tail_laddress();

        for (auto itr = pat.beg; itr < pat.end; ++itr) {
            DEBUG_PRINT(tfm::printfln("** bwsearch for %c ** ", *itr));
            DEBUG_PRINT(tfm::printfln(" - lo_lnode=%s, lo_lofst=%d", m_graph.get_pc(lo_lnode), lo_lofst));
            DEBUG_PRINT(tfm::printfln(" - hi_lnode=%s, hi_lofst=%d", m_graph.get_pc(hi_lnode), hi_lofst));

            const uchar_type chr = *itr;

            // Update the lower bound
            auto [lo_fnode, lo_fofst, dummy1] = m_graph.successor_in_bwsearch(chr, lo_lnode, lo_lofst);
            if (lo_fnode == nullptr) {
                return {nullptr, -1, nullptr, -1};
            }

            // Update the upper bound
            auto [hi_fnode, hi_fofst, dummy2] = m_graph.predecessor_in_bwsearch(chr, hi_lnode, hi_lofst);
            if (hi_fnode == nullptr) {
                return {nullptr, -1, nullptr, -1};
            }

            std::tie(lo_lnode, lo_lofst) = m_graph.get_overlapped_laddress(lo_fnode, lo_fofst);
            std::tie(hi_lnode, hi_lofst) = m_graph.get_overlapped_laddress(hi_fnode, hi_fofst);

            DEBUG_PRINT(tfm::printfln(" - lo_fnode=%s, lo_fofst=%d", m_graph.get_pc(lo_fnode), lo_fofst));
            DEBUG_PRINT(tfm::printfln(" - hi_fnode=%s, hi_fofst=%d", m_graph.get_pc(hi_fnode), hi_fofst));

            // Check if the search range is valid
            if (lo_common::get_basic_order(hi_lnode) < lo_common::get_basic_order(lo_lnode)) {
                return {nullptr, -1, nullptr, -1};
            } else if (hi_lnode == lo_lnode and hi_lofst < lo_lofst) {
                return {nullptr, -1, nullptr, -1};
            }
        }

        DEBUG_PRINT(tfm::printfln("** bwsearch is done **"));
        DEBUG_PRINT(tfm::printfln(" - lo_lnode=%s, lo_lofst=%d", m_graph.get_pc(lo_lnode), lo_lofst));
        DEBUG_PRINT(tfm::printfln(" - hi_lnode=%s, hi_lofst=%d", m_graph.get_pc(hi_lnode), hi_lofst));

        return {lo_lnode, lo_lofst, hi_lnode, hi_lofst};
    }

    // Backward search with toehold lemma
    // Return (lo_lnode, lo_lofst, hi_lnode, hi_lofst, toehold_sae)
    // such that toehold_sae = lo_lnode[lo_lofst].sae
    inline std::tuple<const lnode_type*, offset_type, const lnode_type*, offset_type, size_type>  //
    bwsearch_toehold(range_type<const uchar_type*> pat) const {
        DEBUG_ABORT_IF(pat.beg >= pat.end);

        const auto [em_lnode, em_lofst] = m_graph.get_em_laddress();

        // Search range
        auto [lo_lnode, lo_lofst] = m_graph.get_head_laddress();
        auto [hi_lnode, hi_lofst] = m_graph.get_tail_laddress();

        // We start the SA entry of the head of L-list
        size_type toehold_sae = 0;

        for (auto itr = pat.beg; itr < pat.end; ++itr) {
            DEBUG_PRINT(tfm::printfln("** bwsearch for %c ** ", *itr));
            DEBUG_PRINT(tfm::printfln(" - lo_lnode=%s, lo_lofst=%d", m_graph.get_pc(lo_lnode), lo_lofst));
            DEBUG_PRINT(tfm::printfln(" - hi_lnode=%s, hi_lofst=%d", m_graph.get_pc(hi_lnode), hi_lofst));
            DEBUG_PRINT(tfm::printfln(" - toehold_sae=%d", toehold_sae));

            const uchar_type chr = *itr;

            // Update the lower bound
            auto [lo_fnode, lo_fofst, matched] = m_graph.successor_in_bwsearch(chr, lo_lnode, lo_lofst);
            if (lo_fnode == nullptr) {
                return {nullptr, -1, nullptr, -1, MAX_SIZE_INT};
            }

            if (matched) {
                toehold_sae += 1;
            } else if (lo_lnode == lo_fnode->get_data().get_lnode()) {
                // SPECIAL CASE: Since the cursor (lo_lnode, lo_lofst) indicates $-marker,
                //               the successor is the one simply shifted in the L-node.
                toehold_sae = m_next_sae + 1;
            } else {
                ABORT_IF(!m_graph.is_first_faddress(lo_fnode, lo_fofst));

                auto corr_lnode = lo_fnode->get_data().get_lnode();

                if (m_graph.get_head_lnode() == corr_lnode) {
                    // Then, the SA-entry is always zero (for $)
                    toehold_sae = 1;  // 0+1
                } else if (em_lofst == 0 and corr_lnode == em_lnode) {
                    // Then, the corresponding L character is the next of $.
                    // So, we have the SA-entry in m_next_sae.
                    toehold_sae = m_next_sae + 1;
                } else if (m_graph.get_lexp(corr_lnode) == 1) {
                    // Then, the SA-entry is the last one in the L-node
                    toehold_sae = get_last_sae(corr_lnode) + 1;
                } else {
                    // We extract corr_lnode[0].sae using Successor query for the last SA of the previous L-node
                    auto corr_lprev = corr_lnode->get_prev();
                    const auto [k, v] = m_sae_succ.search(get_last_sae(corr_lprev), MAX_SIZE_INT);
                    ABORT_IF_NE(k, get_last_sae(corr_lprev));
                    toehold_sae = v + 1;
                }
            }

            std::tie(lo_lnode, lo_lofst) = m_graph.get_overlapped_laddress(lo_fnode, lo_fofst);

            // Update the upper bound
            // auto [hi_fnode, hi_fofst, matched] = m_graph.predecessor_in_bwsearch(hi_lnode, hi_lofst, chr);
            auto [hi_fnode, hi_fofst, _] = m_graph.predecessor_in_bwsearch(chr, hi_lnode, hi_lofst);
            if (hi_fnode == nullptr) {
                return {nullptr, -1, nullptr, -1, MAX_SIZE_INT};
            }
            std::tie(hi_lnode, hi_lofst) = m_graph.get_overlapped_laddress(hi_fnode, hi_fofst);

            DEBUG_PRINT(tfm::printfln(" - lo_fnode=%s, lo_fofst=%d", m_graph.get_pc(lo_fnode), lo_fofst));
            DEBUG_PRINT(tfm::printfln(" - hi_fnode=%s, hi_fofst=%d", m_graph.get_pc(hi_fnode), hi_fofst));

            // Check if the search range is valid
            if (lo_common::get_basic_order(hi_lnode) < lo_common::get_basic_order(lo_lnode)) {
                return {nullptr, -1, nullptr, -1, MAX_SIZE_INT};
            } else if (hi_lnode == lo_lnode and hi_lofst < lo_lofst) {
                return {nullptr, -1, nullptr, -1, MAX_SIZE_INT};
            }
        }

        DEBUG_PRINT(tfm::printfln("** bwsearch is done **"));
        DEBUG_PRINT(tfm::printfln(" - lo_lnode=%s, lo_lofst=%d", m_graph.get_pc(lo_lnode), lo_lofst));
        DEBUG_PRINT(tfm::printfln(" - hi_lnode=%s, hi_lofst=%d", m_graph.get_pc(hi_lnode), hi_lofst));
        DEBUG_PRINT(tfm::printfln(" - toehold_sae=%d", toehold_sae));

        return {lo_lnode, lo_lofst, hi_lnode, hi_lofst, toehold_sae};
    }

    //! Get the next L-cursor.
    inline auto get_next(const lnode_type* lnode, offset_type lofst) const
        -> std::tuple<const lnode_type*, offset_type> {
        if (m_graph.is_last_laddress(lnode, lofst)) {
            lnode = lnode->get_next();
            lofst = 0;
        } else {
            lofst += 1;
        }
        return {lnode, lofst};
    }

    //! Get SA[i+1] where SA[i] is the given value.
    inline size_type get_next_sae(size_type sae) const {
        const auto [k, v] = m_sae_succ.search(sae, MAX_SIZE_INT);
        if (k == MAX_SIZE_INT) {
            return MAX_SIZE_INT;
        }
        return v - (k - sae);
    }

    //! Get SA[i+1], where SA[i] is the given value, and the corresponding L-cursor.
    inline auto get_next_sae(const lnode_type* lnode, offset_type lofst, size_type sae) const
        -> std::tuple<const lnode_type*, offset_type, size_type> {
        const auto [k, v] = m_sae_succ.search(sae, MAX_SIZE_INT);
        if (k == MAX_SIZE_INT) {
            return {nullptr, 0, 0};
        }

        const size_type next_sae = v - (k - sae);

        // DEBUG_PRINT(tfm::printfln("** get_next_sae **"));
        // DEBUG_PRINT(tfm::printfln(" - succ(%d)=(%d,%d)", sae, k, v));
        // DEBUG_PRINT(tfm::printfln(" - next_sae=%d", next_sae));

        if (m_graph.is_last_laddress(lnode, lofst)) {
            lnode = lnode->get_next();
            lofst = 0;
        } else {
            lofst += 1;
        }

        return {lnode, lofst, next_sae};
    }

    //! Divide the fat nodes in Case L/F, while updating SA-entries.
    inline bool divide_fat_node() {
        auto [div_lnode, lfm_lnode] = m_graph.divide_fat_lnode();
        if (div_lnode) {
            update_sa_entries_in_LF(div_lnode, lfm_lnode, true);
            return true;
        }
        std::tie(div_lnode, lfm_lnode) = m_graph.divide_fat_fnode();
        if (div_lnode) {
            update_sa_entries_in_LF(div_lnode, lfm_lnode, false);
            return true;
        }
        return false;
    }

    //! Update SA-entries in Step L/F.
    //!  - div_lnode: The before part of the divided L-nodes
    //!  - lfm_lnode: The L-node whose last character is the next/previous one
    //!               of the last one of div_lnode, in the original text.
    inline void update_sa_entries_in_LF(lnode_type* div_lnode, lnode_type* lfm_lnode, bool in_L) {
        // New $'s L-node
        auto [em_lnode, em_lofst] = m_graph.get_em_laddress();

        // The SA-entry of the after part is stored in div_lnode.
        auto new_lnode = div_lnode->get_next();
        new_lnode->get_data().set_sae(div_lnode->get_data().get_sae());

        if (em_lnode == new_lnode and em_lofst == 0) {
            // If the head of the latter divided L-node is the new $'s position,
            // then we can update div_lnode.sae by the tmp SA entries (and lfm_lnode can be invalid).
            div_lnode->get_data().set_sae(m_prev_sae);
        } else {
            div_lnode->get_data().set_sae(lfm_lnode->get_data().get_sae() + (in_L ? 1 : -1));
            const size_type prev_sae = div_lnode->get_data().get_sae();
            const size_type next_sae = get_next_sae(prev_sae);
            m_sae_succ.insert(prev_sae, next_sae);
        }

        // DEBUG_PRINT(tfm::printfln("\n** divide_fat_node **"));
        // DEBUG_PRINT(tfm::printfln(" - div_lnode=%s", m_graph.get_lnode_str(div_lnode)));
        // DEBUG_PRINT(tfm::printfln(" - new_lnode=%s", m_graph.get_lnode_str(new_lnode)));
        // DEBUG_PRINT(tfm::printfln(" - lfm_lnode=%s", m_graph.get_lnode_str(lfm_lnode)));
    }
};

}  // namespace rcomp