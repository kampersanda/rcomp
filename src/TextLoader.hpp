/**
 * @file TextLoader.hpp
 * @author Shunsuke Kanda (kampersanda)
 */
#pragma once

#include "utils.hpp"

namespace rcomp {

/**
 * A class for text loader.
 *
 * @tparam t_EnableOnline Load the text character by character?
 */
template <bool t_EnableOnline>
class TextLoader;

template <>
class TextLoader<true> {
  public:
    static constexpr bool IN_ONLINE = true;

  private:
    std::ifstream m_ifs;

  public:
    TextLoader(const std::string& input_path) : m_ifs(utils::make_ifstream(input_path)) {}

    virtual ~TextLoader() = default;

    inline bool eof() {
        return m_ifs.peek() == std::ios::traits_type::eof();
    }
    inline uchar_type get() {
        char c = {};
        m_ifs.get(c);
        return static_cast<uchar_type>(c);
    }
};

template <>
class TextLoader<false> {
  public:
    static constexpr bool IN_ONLINE = false;

  private:
    std::vector<uchar_type> m_text;
    size_type               m_pos = 0;

  public:
    TextLoader(const std::string& input_path) : m_text(utils::load_text(input_path)) {}

    virtual ~TextLoader() = default;

    inline bool eof() {
        return m_pos == m_text.size();
    }
    inline uchar_type get() {
        return m_text[m_pos++];
    }
};

}  // namespace rcomp