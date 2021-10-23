/**
 * @file Bwt_Naive.hpp
 */
#pragma once

#include <array>
#include <string_view>

#include "basics.hpp"
#include "utils.hpp"

#include "Rindex_Interface.hpp"

namespace rcomp {

/**
 * A BWT class in a naive manner (for test).
 *
 * @note Do NOT use this for a large text.
 */
class Bwt_Naive : public Rindex_Interface {
  public:
    struct node_type {
        uchar_type chr;
        size_type sae;  // SA entry
    };

  private:
    std::array<size_type, 256> m_occ;
    std::vector<node_type> m_nodes;
    size_type m_empos = 0;

  public:
    //! Default constructor
    Bwt_Naive() = default;

    //! Default destructor
    virtual ~Bwt_Naive() = default;

    //! Copy constructor (deleted)
    Bwt_Naive(const Bwt_Naive&) = delete;

    //! Copy constructor (deleted)
    Bwt_Naive& operator=(const Bwt_Naive&) = delete;

    //! Move constructor
    Bwt_Naive(Bwt_Naive&&) noexcept = default;

    //! Move constructor
    Bwt_Naive& operator=(Bwt_Naive&&) noexcept = default;

    //! Check if the data structure is empty.
    inline bool is_empty() const override {
        return m_nodes.empty();
    }

    //! Get the number of stored characters.
    inline size_type get_num_chars() const override {
        return m_nodes.size();
    }

    //! Get the allocated memory in bytes.
    inline size_type get_memory_in_bytes(bool include_this = true) const override {
        const size_type this_bytes = sizeof(*this) * include_this;
        return this_bytes + m_nodes.capacity() * sizeof(node_type);
    }

    //! Print the statistics related to memory.
    void show_memory_statistics() const override {}

    //! Print the detailed statistics.
    void show_detailed_statistics() const override {}

    //! Print the statistics measured with internal monitors.
    void show_monitored_statistics() const override {}

    //! Extend the RLBWT text by appending the given character (i.e., T := c + T).
    void extend(const uchar_type new_chr) override {
        if (m_nodes.empty()) {  // init insertion?
            m_occ.fill(0);
            m_occ[END_MARKER] += 1;
            m_occ[new_chr] += 1;
            m_nodes.push_back(node_type{new_chr, 0});
            m_nodes.push_back(node_type{END_MARKER, 1});
            m_empos = 1;
            return;
        }

        const size_type new_sae = m_nodes[m_empos].sae + 1;

        m_occ[new_chr] += 1;
        m_nodes[m_empos].chr = new_chr;

        // rank+occ should be the result before inserting 'chr' to the BWT.
        // But, since get_rank and get_occ are exclusive to the arguments, the result will be the same.
        m_empos = get_rank(new_chr, m_empos) + get_occ(new_chr);
        m_nodes.insert(m_nodes.begin() + m_empos, node_type{END_MARKER, new_sae});
    }

    //! Count the number of occurrences for the query pattern.
    size_type count(range_type<const uchar_type*> pat) const override {
        size_type n = 0;
        locate(pat, [&](size_type) { n++; });
        return n;
    }

    //! Report the positions of occurrence for the query pattern via the callback function.
    void locate(range_type<const uchar_type*> pat, const std::function<void(size_type)>& fn) const override {
        auto rng = bw_search(pat);
        for (size_type i = rng.beg; i < rng.end; i++) {
            fn(m_nodes[i].sae);
        }
    }

    //! Extract SA-entries from the head of the BWT-text.
    void extract_sa_entries(const std::function<void(size_type)>& fn) const override {
        for (auto& n : m_nodes) {
            fn(n.sae);
        }
    }

    //! Output the original text (in the input order) via the callback function.
    void decode_text(const std::function<void(uchar_type)>& fn) const override {
        if (m_nodes.empty()) {
            return;
        }
        for (size_type pos = 0; pos != m_empos;) {
            fn(m_nodes[pos].chr);
            pos = get_lfpos(pos);
        }
    }

    //! Output the RLBWT text via the callback function.
    void output_runs(const std::function<void(const run_type&)>& fn) const override {
        ABORT_IF_EQ(m_empos, 0);

        uchar_type chr = m_nodes[0].chr;
        size_type exp = 1;

        for (size_type i = 1; i < m_nodes.size(); i++) {
            if (chr == m_nodes[i].chr) {
                exp += 1;
            } else {
                fn(make_run(chr, exp));
                chr = m_nodes[i].chr;
                exp = 1;
            }
        }
        fn(make_run(chr, exp));
    }

    //! Test the data structure.
    void test() const override {}

  private:
    // exclusive pos
    size_type get_rank(uchar_type chr, size_type pos) const {
        size_type cnt = 0;
        for (size_type i = 0; i < pos; i++) {
            if (m_nodes[i].chr == chr) {
                cnt += 1;
            }
        }
        return cnt;
    }

    // get_occ for c<chr (assuming end_marker is contained)
    size_type get_occ(uchar_type chr) const {
        return std::accumulate(m_occ.data(), m_occ.data() + chr, size_type(0));
    }

    size_type get_lfpos(size_type pos) const {
        const uchar_type chr = m_nodes[pos].chr;
        return get_occ(chr) + get_rank(chr, pos);
    }

    range_type<size_type> bw_search(range_type<const uchar_type*> pat) const {
        if (m_occ[*pat.beg] == 0) {
            return {MAX_SIZE_INT, MAX_SIZE_INT};
        }

        size_type bpos = get_occ(*pat.beg);
        size_type epos = get_occ(*pat.beg + 1);

        for (auto it = pat.beg + 1; it < pat.end and bpos < epos; ++it) {
            const uchar_type chr = *it;
            bpos = get_occ(chr) + get_rank(chr, bpos);
            epos = get_occ(chr) + get_rank(chr, epos);
        }

        return {bpos, epos};
    }
};

}  // namespace rcomp