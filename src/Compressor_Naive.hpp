/**
 * @file Compressor_Naive.hpp
 * @author Shunsuke Kanda (kampersanda)
 */
#pragma once

#include "basics.hpp"

namespace rcomp {

// Online BWT compression in a naive manner.
// As this is just for testing the compression mechanism, do NOT use it for large texts.
class Compressor_Naive {
  private:
    std::vector<uchar_type> m_bwt;
    size_type               m_em_pos = 0;

  public:
    Compressor_Naive() = default;

    virtual ~Compressor_Naive() = default;

    //! Copy constructor (deleted)
    Compressor_Naive(const Compressor_Naive&) = delete;

    //! Copy constructor (deleted)
    Compressor_Naive& operator=(const Compressor_Naive&) = delete;

    //! Move constructor
    Compressor_Naive(Compressor_Naive&&) noexcept = default;

    //! Move constructor
    Compressor_Naive& operator=(Compressor_Naive&&) noexcept = default;

    void extend(uchar_type chr) {
        if (m_bwt.empty()) {  // init insertion?
            m_bwt.push_back(chr);
            m_em_pos = 1;
            return;
        }

        if (m_em_pos == m_bwt.size()) {
            m_bwt.push_back(chr);
        } else {
            m_bwt.insert(m_bwt.begin() + m_em_pos, chr);
        }

        // The rank+occ should be the result before inserting 'chr' to the BWT.
        // But, since rank and occ are exclusive to the arguments, the result will be the same.
        m_em_pos = rank(chr, m_em_pos) + occ(chr);
    }

    void print(std::ostream& os) const {
        os << "[L]";
        for (size_type i = 0; i < m_bwt.size(); i++) {
            os << " " << to_print(m_bwt[i]);
        }
        os << std::endl;
    }

    void output_runs(const std::function<void(const run_type&)>& fn) const {
        ABORT_IF_EQ(m_em_pos, 0);

        uchar_type chr = m_bwt[0];
        size_type  exp = 1;

        for (size_type i = 1; i < m_bwt.size(); i++) {
            if (i == m_em_pos) {
                fn(make_run(chr, exp));
                chr = END_MARKER;
                exp = 1;
            }
            if (chr == m_bwt[i]) {
                exp += 1;
            } else {
                fn(make_run(chr, exp));
                chr = m_bwt[i];
                exp = 1;
            }
        }
        fn(make_run(chr, exp));

        if (m_em_pos == m_bwt.size()) {
            fn(make_run(END_MARKER, 1));
        }
    }

  private:
    // exclusive pos
    size_type rank(uchar_type chr, size_type pos) const {
        size_type cnt = 0;
        for (size_type i = 0; i < pos; i++) {
            if (m_bwt[i] == chr) {
                cnt += 1;
            }
        }
        return cnt;
    }

    // occ for c<chr (assuming end_marker is contained)
    size_type occ(uchar_type chr) const {
        size_type cnt = 1;  // for $
        for (size_type i = 0; i < m_bwt.size(); i++) {
            if (m_bwt[i] < chr) {
                cnt += 1;
            }
        }
        return cnt;
    }
};

}  // namespace rcomp
