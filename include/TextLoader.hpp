/**
 * @file TextLoader.hpp
 */
#pragma once

#include <stddef.h>
#include <stdio.h>

#include "abort_if.hpp"
#include "basics.hpp"

namespace rcomp {

class TextLoader {
  public:
    static constexpr size_t BUF_SIZE = 2 << 10;

  private:
    const bool m_reverse;

    FILE* m_fp = nullptr;
    char m_buf[BUF_SIZE];

    size_t m_i = 0;  // of file
    size_t m_j = 0;  // of buf
    size_t m_filesize = 0;

    bool m_eof = false;

  public:
    TextLoader(const std::string& filepath, bool reverse = false) : m_reverse(reverse) {
        m_fp = fopen(filepath.c_str(), "rb");
        ABORT_IF(m_fp == nullptr);

        ABORT_IF(fseek(m_fp, 0L, SEEK_END) != 0);

        auto ret = ftell(m_fp);
        ABORT_IF_EQ(ret, -1L);
        m_filesize = static_cast<size_t>(ret);

        if (!m_reverse) {
            ABORT_IF(fseek(m_fp, 0L, SEEK_SET) != 0);
        }
    }

    virtual ~TextLoader() {
        if (m_fp) {
            fclose(m_fp);
        }
    }

    inline bool eof() const {
        return m_eof;
    }

    inline uchar_type next() {
        ABORT_IF(eof());
        if (!m_reverse) {
            return next_in_forward();
        } else {
            return next_in_backward();
        }
    }

  private:
    inline uchar_type next_in_forward() {
        if (m_j == 0) {
            const auto readsize = std::min<size_t>(BUF_SIZE, m_filesize - m_i);
            ABORT_IF_EQ(fread(m_buf, 1, readsize, m_fp), 0);
        }
        if ((++m_i) == m_filesize) {
            m_eof = true;
        }
        const auto c = static_cast<uchar_type>(m_buf[m_j]);
        if ((++m_j) == BUF_SIZE) {
            m_j = 0;
        }
        return c;
    }

    inline uchar_type next_in_backward() {
        if (m_j == 0) {
            const auto readsize = std::min<int64_t>(BUF_SIZE, m_filesize - m_i);
            fseek(m_fp, -readsize, SEEK_CUR);
            ABORT_IF_EQ(fread(m_buf + BUF_SIZE - readsize, 1, readsize, m_fp), 0);
            fseek(m_fp, -readsize, SEEK_CUR);
        }
        if ((++m_i) == m_filesize) {
            m_eof = true;
        }
        const auto c = static_cast<uchar_type>(m_buf[BUF_SIZE - m_j - 1]);
        if ((++m_j) == BUF_SIZE) {
            m_j = 0;
        }
        return c;
    }
};

}  // namespace rcomp