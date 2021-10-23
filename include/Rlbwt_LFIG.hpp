/**
 * @file Rlbwt_LFIG.hpp
 */
#pragma once

#include "LFIntervalGraph.hpp"

#include "Rlbwt_Interface.hpp"

#ifdef ENABLE_STAT_MONITOR
#include "Timer.hpp"
#endif

namespace rcomp {

/**
 * A RLBWT class with the straightforward LF-interval graph.
 *
 * @tparam t_Graph The LF-interval graph.
 */
template <class t_Graph>
class Rlbwt_LFIG : public Rlbwt_Interface {
  public:
    using this_type = Rlbwt_LFIG<t_Graph>;
    using graph_type = t_Graph;

    static_assert(!t_Graph::ldata_type::WITH_SAE);

  private:
    graph_type m_graph;

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
    Rlbwt_LFIG() = default;

    //! Default destructor
    virtual ~Rlbwt_LFIG() = default;

    //! Copy constructor (deleted)
    Rlbwt_LFIG(const Rlbwt_LFIG&) = delete;

    //! Copy constructor (deleted)
    Rlbwt_LFIG& operator=(const Rlbwt_LFIG&) = delete;

    //! Move constructor
    Rlbwt_LFIG(Rlbwt_LFIG&&) noexcept = default;

    //! Move constructor
    Rlbwt_LFIG& operator=(Rlbwt_LFIG&&) noexcept = default;

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
        return this_bytes + m_graph.get_memory_in_bytes(false);
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
            return;
        }

        // Case L or F
        STAT_MONITOR(m_timer_4_case.start();)
        [[maybe_unused]] size_type num_divided = 0;
        while (m_graph.divide_fat_node()) {
            num_divided += 1;
        }
        STAT_MONITOR(m_ns_case_LF += m_timer_4_case.stop_and_get_ns();)
        STAT_MONITOR(m_num_case_LF += num_divided;)

        // Phase A
        STAT_MONITOR(m_timer_4_case.start();)
        auto ins_mode = m_graph.check_mergable(new_chr);

        // Phase B--F
        if (m_graph.is_mergable_mode(ins_mode)) {
            m_graph.extend_with_merge(new_chr, ins_mode);
            STAT_MONITOR(m_ns_case_C += m_timer_4_case.stop_and_get_ns(););
            STAT_MONITOR(m_num_case_C += 1;)
        } else {
            ins_mode = m_graph.divide_if_need(ins_mode);
            m_graph.extend_with_split(new_chr, ins_mode);
            STAT_MONITOR(m_ns_case_B += m_timer_4_case.stop_and_get_ns();)
            STAT_MONITOR(m_num_case_B += 1;)
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
};

}  // namespace rcomp